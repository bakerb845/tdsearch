#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
//#include <mpi.h>
#include "tdsearch.h"
#include "prepmt/prepmt_pickFile.h"
#include "iscl/iscl/iscl.h"
#include "iscl/time/time.h"

#define DMIN_DIST 30.0
#define DMAX_DIST 95.0
#define PROGRAM_NAME "tdsearch"
static int parseCommands(int argc, char *argv[], char iniFile[PATH_MAX]);

int main(int argc, char *argv[])
{
    struct tdSearchData_struct data;
    struct tdSearch_struct tds;
    struct tdSearchEventParms_struct event;
    struct tdSearchHudson_struct ffGrns;
    struct tdSearchGreens_struct grns;
    struct sacData_struct synth;
    double cutEnd, cutStart, targetDt;
    char iniFile[PATH_MAX], synthName[PATH_MAX], heatMap[PATH_MAX],
         pickFile[PATH_MAX];
    int ierr, iobs;
    bool lsetNewPicks, lusePickFile;

    // Fire up MPI 
    // int provided;
    //MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    // Initialize
    memset(iniFile, 0, PATH_MAX*sizeof(char));
    memset(&event, 0, sizeof(struct tdSearchEventParms_struct));
    memset(&data, 0, sizeof(struct tdSearchData_struct));
    memset(&ffGrns, 0, sizeof(struct tdSearchHudson_struct));
    memset(&grns, 0, sizeof(struct tdSearchGreens_struct));
    memset(&tds, 0, sizeof(struct tdSearch_struct));
    iscl_init();
    // Parse the input commands
    parseCommands(argc, argv, iniFile);
    // Initialize from the ini file - event information
    printf("%s: Initializing from %s...\n", PROGRAM_NAME, iniFile);
    ierr = tdsearch_event_initializeFromIniFile(iniFile, &event);
    if (ierr != 0)
    {
        printf("%s: Failed to read event information\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Initialize grid search parameters from ini file
    ierr = tdsearch_gridSearch_initializeFromFile(iniFile, &tds);
    if (ierr != 0)
    {
        printf("%s: Failed to read grid search parameters\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Initialize data to process from ini file
    ierr = tdsearch_data_initializeFromFile(iniFile,
                                            &data);
    if (ierr != 0)
    {
        printf("%s: Failed to read the data\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    ierr = tdsearch_data_getDefaultDTAndWindowFromIniFile(iniFile,
                                                          &targetDt,
                                                          &cutStart,
                                                          &cutEnd);
    if (ierr != 0)
    {   
        printf("%s: Failed to read the processing defaults\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    ierr = tdsearch_data_getPickStrategy(iniFile,
                                         &lsetNewPicks, &lusePickFile,
                                         pickFile);
    // Initialize the forward modeling parameters for hudson and the
    // source time function to be used in hpulse
    ierr =  tdsearch_hudson_initializeParametersFromIniFile(iniFile, &ffGrns);
    if (ierr != 0)
    {
        printf("%s: Failed to read hudson/hpulse parameters\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Get the Green's functions pre-processing commands
    ierr = tdsearch_greens_setPreprocessingCommandsFromIniFile(iniFile,
                                                               data.nobs,  
                                                               &grns);
    if (ierr != 0)
    {
        printf("%s: Failed to get greens preprocessing commands\n",
               PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Set the event information on the data
    ierr = tdsearch_data_setEventInformation(event.latitude,
                                             event.longitude,
                                             event.depth,
                                             event.time, &data);
    if (ierr != 0)
    {
        printf("%s: Failed to set event info\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Now that I have the event information verify the data are in bounds
    ierr = tdsearch_data_verifyDistances(DMIN_DIST, DMAX_DIST, &data);
    if (ierr != 0)
    {
        printf("%s: Error checking distances\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Set the theoretical pick times on KA and A
    ierr = tdsearch_data_setPicks(NULL, "ak135",
                                  lsetNewPicks, lusePickFile, pickFile,
                                  &data);
    if (ierr != 0)
    {
        printf("%s: Error setting picks\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Set the forward modeling grid on ffGrns
    ierr =  tdsearch_hudson_setGrid(tds.ntstar, tds.tstar,
                                    tds.ndepth, tds.depths,
                                    &ffGrns);
    if (ierr != 0)
    {
        printf("%s: Failed to set forward modeling grid on ffGrns\n",
               PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    //------------------------------------------------------------------------//
    // JEFF - here you might consider breaking this and letting the user work //
    // interactively process the data and make picks.  But I'm going to       //
    // automatically process the data.  To do that I first need to modify the //
    // generic processing commands to conform with the data.                  //
    //------------------------------------------------------------------------//
    printf("%s: Processing data...\n", PROGRAM_NAME);
    time_tic();
 tdsearch_data_modifyProcessingCommands(ffGrns.modelingParms.hpulse96Parms.iodva,
                                        cutStart, cutEnd, targetDt, &data);
    ierr = tdsearch_data_process(&data);
    if (ierr != 0)
    {
        printf("%s: Failed to process data!\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    printf("%s: Processing time %4.2f (s)\n", PROGRAM_NAME, time_toc());
tdsearch_data_writeFiles("prepData", NULL, data);
    // Compute the fundamental fault Green's functions (e.g., ZDD, ZDS, etc.).
    // JEFF - This would be like what you read from the h5 Green's functions
    // archive in that they have no context to the data other than they are
    // at an appropriate distance.
    printf("%s: Computing fundamental-fault Green's functions...\n",
           PROGRAM_NAME);
    time_tic();
    ierr = tdsearch_hudson_computeGreensFF(data, &ffGrns);
    if (ierr != 0)
    {
        printf("%s: Error computing ff Green's functions\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    printf("%s: Green's functions computation time %5.2f (seconds)\n",
           PROGRAM_NAME, time_toc());
    // Rotate the fundamental faults into the observation frame and give
    // these Green's functions some context w.r.t. to the data
    printf("%s: Putting Green's functions in data context...\n",
           PROGRAM_NAME);
    ierr = tdsearch_greens_ffGreensToGreens(data, ffGrns, &grns);
    if (ierr != 0)
    {
        printf("%s: Error manipulating Green's functions\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // I'm done with the fundamental faults. 
    ierr = tdsearch_hudson_free(&ffGrns);
    if (ierr != 0)
    {
        printf("%s: Failed to free ffGrns\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    // Fix the Green's functions pre-processing commands
    printf("%s: Modifying Green's functions processing commands...\n",
           PROGRAM_NAME);
    ierr = tdsearch_greens_modifyProcessingCommands(
                 ffGrns.modelingParms.hpulse96Parms.iodva,
                 cutStart, cutEnd, targetDt, &grns);
    if (ierr != 0)
    {
        printf("%s: Failed to modify Green's pre-processing commands\n",
               PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    ierr = tdsearch_greens_process(&grns);
/*
tdsearch_greens_writeSelectGreensFunctions("prepData",
                                           0, 1, 3, grns);
*/
    // At this point the data and Green's functions defined on the (t*,depth)
    // grid have been pre-processed.  Now the estimation begins.
    for (iobs=0; iobs<data.nobs; iobs++)
    {
        printf("%s: Setting forward modeling matrices...\n", PROGRAM_NAME);
        time_tic();
        // Set the forward modeling matrices so that I can compute
        // d*u/(|d|_2 |u|_2) = d*(G m)/(|d|_2 |u|_2) 
        // expediently for any given moment tensor
        ierr = tdsearch_gridSearch_setForwardModelingMatrices(iobs, data,
                                                              grns, &tds);
        if (ierr != 0)
        {
            printf("%s: Failed to set forward modeling matrices\n",
                   PROGRAM_NAME);
            return EXIT_FAILURE;
        }
        printf("%s: Forward modeling matrix computation time: %e (s)\n",
               PROGRAM_NAME, time_toc());
        //--------------------------------------------------------------------//
        // JEFF - here you could put an interactive loop.  In this instance   //
        // the user could continually set the moment tensor and perpetually   //
        // re-run.                                                            //
        //--------------------------------------------------------------------//
        // Set the moment tensor
        ierr = tdSearch_gridSearch_setMomentTensorFromEvent(event, &tds);
        if (ierr != 0)
        {
            printf("%s: Failed to set event\n", PROGRAM_NAME);
            return EXIT_FAILURE;
        }
        // Run the t*/depth grid search
        printf("%s: Performing grid search...\n", PROGRAM_NAME);
        ierr = tdSearch_gridSearch_performGridSearch(&tds); 
        if (ierr != 0)
        {
            printf("%s: Failed to compute grid search\n", PROGRAM_NAME);
            return EXIT_FAILURE;
        }
        // Write the heatmap for gnuplot plotting
        // JEFF - you would probably grab tds.xcorr, tds.tstar, and tds.depth
        // and plot a 2D image (maybe qwt_plot_spectrogram).  Note that depth
        // increases down so you'd have to flip the y (depth) axis.
        memset(heatMap, 0, PATH_MAX*sizeof(char));
        sprintf(heatMap, "%s.%s.%s.%s",
                data.obs[iobs].header.knetwk, data.obs[iobs].header.kstnm,
                data.obs[iobs].header.kcmpnm, data.obs[iobs].header.khole);
        ierr = tdsearch_gridSearch_writeHeatMap("prepData", heatMap, tds);
        // Write the optimal synthetic
        // JEFF - you would plot this synthetic next to the data.
        memset(&synth, 0, sizeof(struct sacData_struct));
        ierr = tdSearch_gridSearch_makeSACSynthetic(iobs, tds.itopt, tds.idopt,
                                                    data, grns, tds, &synth);
        memset(synthName, 0, PATH_MAX*sizeof(char));
        sprintf(synthName, "prepData/%s.%s.%s.%s.EST.SAC",
                synth.header.knetwk, synth.header.kstnm, 
                synth.header.kcmpnm, synth.header.khole);
        sacio_writeTimeSeriesFile(synthName, synth); 
        sacio_free(&synth);
    }
    // Free space
    ierr = tdsearch_gridSearch_free(&tds);
    ierr = tdsearch_data_free(&data);
    ierr = tdsearch_greens_free(&grns);
    iscl_finalize();
    // Done with MPI
    //MPI_Finalize();
    return 0;
}

static int parseCommands(int argc, char *argv[], char iniFile[PATH_MAX])
{
    memset(iniFile, 0, PATH_MAX*sizeof(char));
strcpy(iniFile, "tdsearch.ini");
    return 0;
}
