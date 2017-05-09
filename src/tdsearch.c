#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <mpi.h>
#include "tdsearch.h"
#include "iscl/iscl/iscl.h"

#define PROGRAM_NAME "tdsearch"
static int parseCommands(int argc, char *argv[], char iniFile[PATH_MAX]);

int main(int argc, char *argv[])
{
    struct tdSearchData_struct data;
    struct tdSearch_struct tds;
    struct tdSearchEventParms_struct event;
    char iniFile[PATH_MAX];
    int ierr, provided;
    // Fire up MPI 
    //MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    // Initialize
    memset(iniFile, 0, PATH_MAX*sizeof(char));
    memset(&event, 0, sizeof(struct tdSearchEventParms_struct));
    memset(&data, 0, sizeof(struct tdSearchEventParms_struct));
    memset(&tds, 0, sizeof(struct tdSearch_struct));
    iscl_init();
    // Parse the input commands
    parseCommands(argc, argv, iniFile);
    // Initialize from the ini file - event information
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
    ierr = tdsearch_data_initializeFromFile(iniFile, &data);
    if (ierr != 0)
    {
        printf("%s: Failed to read the data\n", PROGRAM_NAME);
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
    // Set the theoretical pick times on KA and A
    ierr = tdsearch_data_setPPickTimeFromTheoreticalTime(NULL, "ak135",
                                                         SAC_FLOAT_A,
                                                         SAC_CHAR_KA,
                                                         &data);
    if (ierr != 0)
    {
        printf("%s: Failed to set theoretical pick times\n", PROGRAM_NAME);
        return EXIT_FAILURE;
    }
    //------------------------------------------------------------------------//
    // JEFF - here you might consider breaking this and letting the user work //
    // interactively process the data and make picks.  But I'm going to       //
    // automatically process the data.                                        //
    //------------------------------------------------------------------------//

    // Generate Green's functions

    // Free space
    ierr = tdsearch_gridSearch_free(&tds);
    ierr = tdsearch_data_free(&data);
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
