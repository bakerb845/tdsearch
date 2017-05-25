#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <iniparser.h>
#include "tdsearch_data.h"
#include "tdsearch_commands.h"
#include "ispl/process.h"
#include "sacio.h"
#include "ttimes.h"
#include "iscl/array/array.h"
#include "iscl/geodetic/geodetic.h"
#include "iscl/log/log.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"
#include "iscl/string/string.h"

/*!
 * @brief Reads the data and pole-zero files from the list specified in the
 *        an ini file.
 *
 * @param[in] iniFile     Name of ini file.
 *
 * @param[out] data       On successful exit contains the input data.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_data_initializeFromFile(const char *iniFile,
                                     struct tdSearchData_struct *data)
{
    const char *fcnm = "tdsearch_data_initializeFromFile\0";
    FILE *flist;
    const char *s;
    char **sacFiles, **sacpzFiles, **csplit, **cmds,
         varname[128], cline[2*PATH_MAX];
    dictionary *ini;
    size_t lenos;
    int ierr, item, k, ncmds, ncmdsWork, nitems, nfiles, nobs;
    bool luseDataList, luseProcessingList;
    nobs = 0;
    nfiles = 0;
    ierr = 0;
    sacFiles = NULL;
    sacpzFiles = NULL;
    memset(data, 0, sizeof(struct tdSearchData_struct));
    if (!os_path_isfile(iniFile))
    {
        log_errorF("%s: Error ini file %s doesn't exist\n", fcnm, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile); 
    luseDataList = iniparser_getboolean(ini, "tdSearch:data:luseDataList\0",
                                        false);
    // Read from external file list
    if (luseDataList)
    {
        nfiles = 0;
        s = iniparser_getstring(ini, "tdSearch:data:dataList\0", NULL);
        if (!os_path_isfile(s))
        {
            log_errorF("%s: Error file list %s does not exist\n", fcnm, s);
            return -1;
        }
        flist = fopen(s, "r");
        while(fgets(cline, 128, flist) != NULL){nfiles = nfiles + 1;}
        if (nfiles < 1)
        {
            log_errorF("%s: Error no files in file list\n", fcnm);
            fclose(flist);
        }
        sacFiles = (char **) calloc((size_t) nfiles, sizeof(char *));
        sacpzFiles = (char **) calloc((size_t) nfiles, sizeof(char *));
        rewind(flist);
        for (k=0; k<nfiles; k++)
        {
            memset(cline, 0, 2*PATH_MAX*sizeof(char));
            fgets(cline, 2*PATH_MAX, flist);
            csplit = string_rsplit(NULL, cline, &nitems);
            if (csplit != NULL)
            {
                if (nitems > 0)
                {
                    if (!os_path_isfile(csplit[0]))
                    {
                        log_errorF("%s: Error data file %s doesn't exist\n",
                                   fcnm, csplit[0]); 
                        goto NEXT_LINE; 
                    }
                    else
                    {
                        sacFiles[nobs] = (char *)
                                         calloc(PATH_MAX, sizeof(char));
                        strcpy(sacFiles[nobs], csplit[0]); 
                        nobs = nobs + 1;
                    }
                }
                if (nitems == 2)
                {
                    if (!os_path_isfile(csplit[1]))
                    {
                        log_errorF("%s: Warning pz file %s doesn't exist\n",
                                   fcnm, csplit[1]);
                    }
                    else
                    {
                        sacpzFiles[nobs-1] = (char *)
                                             calloc(PATH_MAX, sizeof(char));
                        strcpy(sacpzFiles[nobs-1], csplit[1]);
                    }
                }
                NEXT_LINE:;
                for (item=0; item<nitems; item++){free(csplit[item]);}
                free(csplit);
            }
        } // Loop on files in file list
        fclose(flist);
    }
    // Read from internal file list
    else
    {
        nfiles = iniparser_getint(ini, "tdSearch:data:nobs\0", 0); 
        if (nfiles < 1)
        {
            log_errorF("%s: No files to read!\n", fcnm);
            ierr = 1;
            goto ERROR;
        }
        sacFiles = (char **) calloc((size_t) nfiles, sizeof(char *));
        sacpzFiles = (char **) calloc((size_t) nfiles, sizeof(char *));
        for (k=0; k<nfiles; k++)
        {
            memset(varname, 0, 128*sizeof(char));
            sprintf(varname, "tdSearch:data:sacFile_%d", k+1);
            sacFiles[k] = (char *) calloc(PATH_MAX, sizeof(char));
            sacpzFiles[k] = (char *) calloc(PATH_MAX, sizeof(char));
            s = iniparser_getstring(ini, varname, NULL);
            if (!os_path_isfile(s))
            {
                log_errorF("%s: Error sac file %s doesn't exist\n", fcnm, s);
                continue;
            }
            strcpy(sacFiles[nobs], s);
            memset(varname, 0, 128*sizeof(char));
            sprintf(varname, "tdSearch:data:sacpzFile_%d", k+1); 
            s = iniparser_getstring(ini, varname, NULL);
            if (os_path_isfile(s))
            {
                strcpy(sacpzFiles[nobs], s);
            }
            nobs = nobs + 1;
        }
    }
    // Read the data
    ierr = tdsearch_data_readFiles(nfiles,
                                   (const char **) sacFiles,
                                   (const char **) sacpzFiles, data);
    if (ierr != 0)
    {
        log_errorF("%s: Error reading data files\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    // Read the processing commands
    data->cmds = (struct tdSearchDataProcessingCommands_struct *)
                 calloc((size_t) data->maxobs,
                        sizeof(struct tdSearchDataProcessingCommands_struct));
    ncmds = iniparser_getint(ini, "tdSearch:data:nCommands\0", 0);
    luseProcessingList
        = iniparser_getboolean(ini, "tdSearch:data:useProcessingList\0", false);
    if (!luseProcessingList)
    {
        if (ncmds > 0)
        {
            ncmdsWork = ncmds;
            ncmds = 0;
            cmds = (char **) calloc((size_t) ncmdsWork, sizeof(char *));
            for (k=0; k<ncmdsWork; k++)
            {
                memset(varname, 0, 128*sizeof(char));
                sprintf(varname, "tdSearch:data:command_%d", k+1);
                s = iniparser_getstring(ini, varname, NULL);
                if (s == NULL){continue;}
                lenos = strlen(s);
                cmds[ncmds] = (char *) calloc(lenos+1, sizeof(char));
                strcpy(cmds[ncmds], s);
                //printf("%s\n", cmds[ncmds]);
                ncmds = ncmds + 1;
            }
            // Attach these to the data processing structure
            for (k=0; k<data->nobs; k++)
            {
                ierr = tdsearch_data_attachCommandsToObservation(
                           k, ncmds, (const char **) cmds, data);
            }
            // Free space
            for (k=0; k<ncmdsWork; k++)
            {
                if (cmds[k] != NULL){free(cmds[k]);}
            }
            free(cmds);
        }
    }
    else
    {
        log_errorF("%s: processing commands list not yet programmed\n", fcnm);
        ierr = 1;
    }
ERROR:;
    if (nfiles > 0)
    {
        if (sacFiles != NULL)
        {
            for (k=0; k<nfiles; k++){free(sacFiles[k]);}
            free(sacFiles);
        }
        if (sacpzFiles != NULL)
        {
            for (k=0; k<nfiles; k++){free(sacpzFiles[k]);}
            free(sacpzFiles);
        }
    }
    iniparser_freedict(ini);
    return ierr;
}
//============================================================================//
/*!
 * @brief Utility function for attaching processing commands to the observation
 *        structure.
 *
 * @param[in] iobs      C indexed observation number.
 * @param[in] ncmds     Number of processing commands to attach to data
 *                      structure.
 * @param[in] cmds      Commands to attach to data structure [ncmds].
 * 
 * @param[in,out] data  On input contains the maximum number of observations
 *                      allowed.
 *                      On exit data->cmds[iobs] contains the corresponding
 *                      data commands.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_data_attachCommandsToObservation(const int iobs,
                                              const int ncmds,
                                              const char **cmds,
                                              struct tdSearchData_struct *data)
{
    const char *fcnm = "tdsearch_data_attachCommandsToObservation\0";
    int i;
    size_t lenos;
    // Make sure iobs is in bounds 
    if (iobs < 0 || iobs >= data->maxobs)
    {
        log_errorF("%s: Error iobs is out of bounds [0,%d]\n",
                   fcnm, iobs, data->maxobs);
        return -1;
    }
    // Try to handle space allocation if not already done
    if (data->cmds == NULL && data->maxobs > 0)
    {
        data->cmds = (struct tdSearchDataProcessingCommands_struct *)
                     calloc((size_t) data->maxobs,
                        sizeof(struct tdSearchDataProcessingCommands_struct));
    }
    if (data->cmds == NULL)
    {
        log_errorF("%s: Error data->cmds is NULL\n", fcnm);
        return -1; 
    }
    if (ncmds == 0){return 0;}
    data->cmds[iobs].ncmds = ncmds;
    data->cmds[iobs].cmds = (char **) calloc((size_t) ncmds, sizeof(char *));
    for (i=0; i<ncmds; i++)
    {
        lenos = strlen(cmds[i]);
        data->cmds[iobs].cmds[i] = (char *) calloc(lenos+1, sizeof(char));
        strcpy(data->cmds[iobs].cmds[i], cmds[i]); 
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Releases memory on the data structure.
 *
 * @param[out] data    On exit all memory has been freed and all variables on
 *                     data have been set to 0 or NULL.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_data_free(struct tdSearchData_struct *data)
{
    const char *fcnm = "tdsearch_data_free\0";
    int i, ierr, k;
    ierr = 0;
    if (data->maxobs > 0 && data->obs != NULL)
    {
        for (k=0; k<data->nobs; k++)
        {
            ierr += sacio_free(&data->obs[k]);
            if (ierr != 0)
            {
                log_errorF("%s: Error freeing SAC data struct %d\n",
                           fcnm, k+1);
            }
        }
        free(data->obs);
    }
    if (data->cmds != NULL)
    {
        for (k=0; k<data->nobs; k++)
        {
            if (data->cmds[k].cmds != NULL)
            {
                for (i=0; i<data->cmds[k].ncmds; i++)
                {
                    if (data->cmds[k].cmds[i] != NULL)
                    {
                        free(data->cmds[k].cmds[i]);
                    }
                }
                free(data->cmds[k].cmds);
            }
        }
        free(data->cmds);
    }
    memory_free8l(&data->lskip);
    memset(data, 0, sizeof(struct tdSearchData_struct));
    return ierr;
}
//============================================================================//
/*!
 * @brief Reads all data and meta data onto the data structure.
 *
 * @param[in] nfiles    Number of SAC files to read.
 * @param[in] sacfls    List of SAC data files [nfiles].
 * @param[in] pzfls     List of corresponding SAC pole-zero files [nfiles].
 *                      If this is NULL then it will be ignored.
 *
 * @param[out] data     On successful contains the data and corresponding
 *                      pole-zero data.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_data_readFiles(const int nfiles,
                            const char **sacfls,
                            const char **pzfls,
                            struct tdSearchData_struct *data)
{
    const char *fcnm = "tdsearch_data_readFiles\0";
    int ierr, ierrAll, k;
    bool lreadpz;
    // Make sure there is something to read
    ierrAll = 0;
    if (nfiles < 1)
    {
        log_errorF("%s: No data files to read\n", fcnm);
        return -1;
    }
    lreadpz = true;
    if (pzfls == NULL){lreadpz = false;}
    // Ensure that sufficient space exists
    if (nfiles > data->maxobs)
    {
        tdsearch_data_free(data);
        data->obs = (struct sacData_struct *)
                    calloc((size_t) nfiles, sizeof(struct sacData_struct));
        data->lskip = memory_calloc8l(nfiles);
        data->maxobs = nfiles;
    }
    // Load the data files and pole-zero files
    data->nobs = 0; 
    for (k=0; k<nfiles; k++)
    {
        ierr = sacio_readTimeSeriesFile(sacfls[k], &data->obs[data->nobs]);
        if (ierr != 0)
        {
            log_errorF("%s: Failed to read file %s\n", fcnm, sacfls[k]);
            ierrAll = ierrAll + 1;
            continue;
        }
        if (lreadpz)
        {
            if (!os_path_isfile(pzfls[k]))
            {
                log_warnF("%s: No pole-zero file corresponding to %s\n",
                          fcnm, sacfls[k]);
                goto JUST_DATA;
            }
        }
        ierr = sacio_readPoleZeroFile(pzfls[k], &data->obs[data->nobs].pz);
        if (ierr != 0)
        {
            log_errorF("%s: Couldn't read %s pole-zero file\n", fcnm, pzfls[k]);
        }
        JUST_DATA:;
        data->nobs = data->nobs + 1;
    }
    return ierrAll;
}
//============================================================================//
/*!
 * @brief Sets the event latitude, longitude, depth, and event time on the
 *        SAC header.
 *
 * @param[in] evla         Event latitude (degrees).
 * @param[in] evlo         Event longitude (degrees).
 * @param[in] evdp         Event depth (km).
 * @param[in] evtime       Event time (UTC-seconds).
 * @param[in,out] data     On input contains the SAC data.
 *                         On output contains the event information.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_data_setEventInformation(const double evla,
                                      const double evlo,
                                      const double evdp,
                                      const double evtime,
                                      struct tdSearchData_struct *data)
{
    const char *fcnm = "tdsearch_data_setEventInformation\0";
    double az, baz, dist, epoch, gcarc, o, stla, stlo;
    int ierr, k;
    ierr = 0;
    for (k=0; k<data->nobs; k++)
    {
        ierr = sacio_getEpochalStartTime(data->obs[k].header, &epoch);
        o = evtime - epoch;
        if (ierr != 0)
        {
            log_errorF("%s: Failed to get start time\n", fcnm);
            break;
        }
        // Get the station location
        ierr = sacio_getFloatHeader(SAC_FLOAT_STLA, data->obs[k].header, &stla);
        if (ierr != 0)
        {
            log_errorF("%s: Failed to get station latitude\n", fcnm);
            break;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_STLO, data->obs[k].header, &stlo);
        if (ierr != 0)
        {
            log_errorF("%s: Failed to get station longitude\n", fcnm);
            break;
        }
        // Compute the azimuth, back-azimuth, distances 
        geodetic_gps2distanceAzimuth(evla, evlo, stla, stlo,
                                     &dist, &gcarc, &az, &baz);
        sacio_setFloatHeader(SAC_FLOAT_EVLA,  evla,  &data->obs[k].header);
        sacio_setFloatHeader(SAC_FLOAT_EVLO,  evlo,  &data->obs[k].header);
        sacio_setFloatHeader(SAC_FLOAT_EVDP,  evdp,  &data->obs[k].header);
        sacio_setFloatHeader(SAC_FLOAT_GCARC, gcarc, &data->obs[k].header);
        sacio_setFloatHeader(SAC_FLOAT_AZ,    az,    &data->obs[k].header);
        sacio_setFloatHeader(SAC_FLOAT_BAZ,   baz,   &data->obs[k].header);
        sacio_setFloatHeader(SAC_FLOAT_O,     o,     &data->obs[k].header);
    }
    return ierr;
}
//============================================================================//
/*!
 * @brief Utility routine for verifying all data is in an appropriate
 *        modeling distance range.
 *
 * @param[in] dmin      Closest source/station distance (degrees).  This is
 *                      likely 30.
 * @param[in] dmax      Furthest source/station distance (degrees).  This is
 *                      likely 95.
 * @param[in,out] data  On input contains the data and their great circle
 *                      distances defined in the SAC headers.
 *                      On output lskip will be updated if the observation
 *                      is out of the distance bounds.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_data_verifyDistances(const double dmin, const double dmax,
                                  struct tdSearchData_struct *data)
{
    const char *fcnm = "tdsearch_data_verifyDistances\0";
    double gcarc;
    int ierr, k;
    if (dmin > dmax || dmin < 0.0 || dmax > 180.0)
    {
        if (dmin > dmax)
        {
            log_errorF("%s: Error dmin %f > dmax %f\n", fcnm, dmin, dmax);
        }
        if (dmin < 0.0)
        {
            log_errorF("%s: Error dmin %f must be non-negative\n", fcnm, dmin);
        }
        if (dmax > 180.0)
        {
            log_errorF("%s: Error dmax %f must be < 180\n", fcnm, dmax);
        }
        return -1;
    }
    for (k=0; k<data->nobs; k++)
    {
        ierr = sacio_getFloatHeader(SAC_FLOAT_GCARC,
                                    data->obs[k].header, &gcarc); 
        if (ierr != 0)
        {
            log_errorF("%s: Error gcarc not set on header\n", fcnm);
        }
        if (gcarc < dmin || gcarc > dmax)
        {
            log_errorF("%s: Observations %d is out of bounds [%f %f]\n",
                        fcnm, dmin, dmax);
            data->lskip[k] = true;
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the theoretical P pick time.
 *
 * @param[in] dirnm     Directory containing the ttimes models.  If NULL this
 *                      will use the default as specified in the ttimes config
 *                      header file.
 * @param[in] model     ttimes model.  Likely "ak135" or "iasp91".
 * @param[in] data      Data with event depth, event origin time, and the
 *                      source-to-receiver great circle distance for each 
 *                      observation.
 * @param[in] nwork     Max allotted workspace to ptimes.  This must be 
 *                      >= data.nobs;
 *
 * @param[out] ptimes   Epochal theoretical primary arrival times (UTC-seconds)
 *                      for each observation.
 *
 * @result 0 indicates success.
 * 
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_data_computeTheoreticalPPickTimes(
    const char *dirnm, const char *model,
    const struct tdSearchData_struct data,
    const int nwork, double *__restrict__ ptimes)
{
    const char *fcnm = "tdsearch_data_computeTheoreticalPPickTimes\0";
    struct ttimesTravelTime_struct ppick;
    double delta, depth, epoch, o;
    int ierr, k;
    if (nwork < data.nobs || ptimes == NULL)
    {
        log_errorF("%s: Insufficient space for output ptimes\n", fcnm);
        return -1;
    }
    memset(&ppick, 0, sizeof(struct ttimesTravelTime_struct));
    for (k=0; k<data.nobs; k++)
    {
        ierr = sacio_getEpochalStartTime(data.obs[k].header, &epoch);
        if (ierr != 0)
        {
            log_errorF("%s: Failed to get start time\n", fcnm);
            return -1;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_O, data.obs[k].header, &o);
        if (ierr != 0)
        {
            log_errorF("%s: Origin time not set\n", fcnm);
            return -1;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_GCARC,
                                    data.obs[k].header, &delta); 
        if (ierr != 0)
        {   
            log_errorF("%s: Error event distance not set\n", fcnm);
            return -1; 
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_EVDP,
                                    data.obs[k].header, &depth);
        if (ierr != 0)
        {
            log_errorF("%s: Error event depth not set\n", fcnm);
            return -1;
        }
        if (ierr != 0)
        {   
            log_errorF("%s: Error getting origin time\n", fcnm);
            return -1;
        }
        ierr = ttimes_getFirstPPhase(delta, depth, dirnm, model, &ppick);
        if (ierr != 0)
        {
            log_errorF("%s: Error computing P pick time\n", fcnm);
            return -1;
        }
        ptimes[k] = epoch + o + ppick.tt;
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the P pick time on the header from the theoretical travel time.
 *
 * @param[in] dirnm           Directory name where ttimes models reside.
 *                            If NULL then the default directory will be used.
 * @param[in] model           Earth model to be used in travel time computation.
 *                            If NULL then ak135 will be used.
 * @param[in] pickHeaderTime  Header variable identifier to hold the pick
 *                            time.  Likely SAC_FLOAT_A for a primary arrival,
 *                            but also may be SAC_FLOAT_T0-SAC_FLOAT_T9
 *                            for other phases.
 * @param[in] pickHeaderName  Header variable identifer to hold the phase name.
 *                            Likely SAC_CHAR_KA for a primary arrival,
 *                            but also may be SAC_CHAR_KT0-SAC_CHAR_KT9
 *                            for other phases.
 *
 * @param[in,out] data        On input holds the data.
 *                            On output holds the theoretical primary pick time.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_data_setPPickTimeFromTheoreticalTime(
    const char *dirnm, const char *model,
    const enum sacHeader_enum pickHeaderTime,
    const enum sacHeader_enum pickHeaderName,
    struct tdSearchData_struct *data)
{
    const char *fcnm = "tdsearch_data_setPPickTimeFromTheoreticalTime\0";
    double *ptimes;
    char phaseName[8];
    int ierr, k, nobs;
    //------------------------------------------------------------------------//
    //
    // Initialize and get number of observations
    ierr = 0;
    ptimes = NULL;
    memset(phaseName, 0, 8*sizeof(char));
    phaseName[0] = 'P';
    nobs = data->nobs;
    if (nobs < 1)
    {
        log_errorF("%s: Error - no data\n", fcnm);
        return -1;
    }
    // Compute the theoretical primary arrival times
    ptimes = memory_calloc64f(nobs);
    ierr = tdsearch_data_computeTheoreticalPPickTimes(dirnm, model, *data,
                                                      nobs, ptimes); 
    if (ierr != 0)
    {
        log_errorF("%s: Error computing theorietical P times\n", fcnm);
        goto ERROR;
    }
    // Set the theoretical arrival times on the header
    for (k=0; k<data->nobs; k++)
    {
        ierr = sacio_setEpochalPickOnHeader(ptimes[k], phaseName,
                                            pickHeaderTime,
                                            pickHeaderName,
                                            &data->obs[k].header);
        if (ierr != 0)
        {
            log_errorF("%s: Failed to set pick time on header\n", fcnm);
            goto ERROR;
        }
    }
ERROR:;
    memory_free64f(&ptimes);
    return ierr;
}
//============================================================================//
/*!
 * @brief Some ad-hoc rules for fixing the data processing commands. 
 *
 * @param[in] cut0        Time relative to pick to begin cut (seconds).
 * @param[in] cut1        Time relative to pick to end cut (seconds).
 * @param[in] targetDt    Target sampling period (seconds).
 * @param[in,out] data    On input contains the data and processing commands.
 *                        On output the generic processing commands have
 *                        been modified so that they are parsable by ispl.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_data_modifyProcessingCommands(
    const double cut0, const double cut1, const double targetDt,
    struct tdSearchData_struct *data)
{
    const char *fcnm = "tdsearch_data_modifyProcessingCommands\0";
    struct tdSearchModifyCommands_struct options;
    const char **cmds;
    char **newCmds;
    int i, ierr, k, ncmds;
    //------------------------------------------------------------------------//
    ierr = 0;
    memset(&options, 0, sizeof(struct tdSearchModifyCommands_struct));
    for (k=0; k<data->nobs; k++)
    {
        newCmds = NULL;
        ncmds = data->cmds[k].ncmds;
        if (ncmds < 1){continue;}
        options.cut0 = cut0;
        options.cut1 = cut1;
        options.targetDt = targetDt;
        options.ldeconvolution = true;
        options.iodva = 1; // TODO change me
        cmds = (const char **) data->cmds[k].cmds;
        newCmds = tdsearch_commands_modifyCommands(ncmds, (const char **) cmds,
                                                   options,
                                                   data->obs[k], &ierr);
        if (ierr != 0)
        {
            log_errorF("%s: Error modifying commands for obs: %d\n",
                       fcnm, k+1); 
            break;
        }
        // Reset the command structure and free newCmds
        for (i=0; i<ncmds; i++)
        { 
            if (data->cmds[k].cmds[i] != NULL){free(data->cmds[k].cmds[i]);}
            data->cmds[k].cmds[i] = (char *) calloc(MAX_CMD_LEN, sizeof(char));
            strcpy(data->cmds[k].cmds[i], newCmds[i]);
            free(newCmds[i]);
        }
        free(newCmds);
    }
    return ierr;
}
//============================================================================//
int tdsearch_data_writeFiles(const char *outdir,
                             const char *suffix,
                             const struct tdSearchData_struct data)
{
    const char *fcnm = "tdsearch_data_writeFiles\0";
    char fname[PATH_MAX], root[PATH_MAX], kcmpnm[8],
         knetwk[8], khole[8], kstnm[8];
    int ierr, k;
    size_t lendir;
    bool lsuffix;
    // Set the root output directory
    lendir = 0;
    memset(root, 0, PATH_MAX*sizeof(char));
    if (outdir == NULL)
    {
        strcpy(root, "./\0");
    }
    else
    {
        lendir = strlen(outdir);
        if (lendir > 0)
        {
            strcpy(root, outdir);
            if (root[lendir-1] != '/'){root[lendir] = '/';}
        }
        else
        {
            strcpy(root, "./\0");
        }
    }
    ierr = os_makedirs(root);
    if (ierr != 0)
    {
        log_errorF("%s: Error making output directory: %s\n", fcnm, root);
        return -1;
    }
    lsuffix = false;
    if (suffix != NULL)
    {
        if (strlen(suffix) > 0){lsuffix = true;}
    }
    for (k=0; k<data.nobs; k++)
    {
        memset(fname, 0, PATH_MAX*sizeof(char));
        sacio_getCharacterHeader(SAC_CHAR_KNETWK, data.obs[k].header, knetwk);
        sacio_getCharacterHeader(SAC_CHAR_KSTNM,  data.obs[k].header, kstnm);
        sacio_getCharacterHeader(SAC_CHAR_KCMPNM, data.obs[k].header, kcmpnm);
        sacio_getCharacterHeader(SAC_CHAR_KHOLE,  data.obs[k].header, khole);
        if (lsuffix)
        {
            sprintf(fname, "%s%s.%s.%s.%s.%s.SAC",
                    root, knetwk, kstnm, kcmpnm, khole, suffix);
        }
        else
        {
            sprintf(fname, "%s%s.%s.%s.%s.SAC",
                    root, knetwk, kstnm, kcmpnm, khole); 
        }
        sacio_writeTimeSeriesFile(fname, data.obs[k]);
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Processes the data according to the commands on the data data
 *        structure.
 */
int tdsearch_data_process(struct tdSearchData_struct *data)
{
    const char *fcnm = "tdsearch_data_process\0";
    struct serialCommands_struct *commands;
    double dt, dt0, epoch, epoch0, time, *ycopy;
    int *nyAll, i, i0, ierr, k, npts0, nq, nwork;
    bool lnewDt, lnewStartTime;
    const int nTimeVars = 12;
    const enum sacHeader_enum timeVars[12]
       = {SAC_FLOAT_A, SAC_FLOAT_O, 
          SAC_FLOAT_T0, SAC_FLOAT_T1, SAC_FLOAT_T2, SAC_FLOAT_T3,
          SAC_FLOAT_T4, SAC_FLOAT_T5, SAC_FLOAT_T6, SAC_FLOAT_T7,
          SAC_FLOAT_T8, SAC_FLOAT_T9};
    // No data
    ycopy = NULL;
    nyAll = NULL;
    if (data->nobs < 1){return 0;}
    commands = (struct serialCommands_struct *)
               calloc((size_t)data->nobs, sizeof(struct serialCommands_struct));
    // Set the processing commands and the data
    for (k=0; k<data->nobs; k++)
    {
        ierr = process_stringsToSerialCommandsOptions(data->cmds[k].ncmds,
                                      (const char **) data->cmds[k].cmds,
                                      &commands[k]);
        if (ierr != 0)
        {
            log_errorF("%s: Error setting serial command string\n", fcnm);
            goto ERROR;
        }
        process_setSerialCommandsData64f(data->obs[k].npts,
                                         data->obs[k].data, 
                                         &commands[k]);
    }
    // Process the data
    ierr = process_applyMultipleSerialCommands(data->nobs, commands);
    if (ierr != 0)
    {
        log_errorF("%s: Error applying serial commands chains to data\n", fcnm);
        goto ERROR;
    }
    // Get the workspace
    nwork =-1;
    nyAll = memory_calloc32i(data->nobs);
    for (k=0; k<data->nobs; k++)
    {
        process_getSerialCommandsData64f(commands[k], -1, &nyAll[k], NULL); 
        nwork = MAX(nwork, nyAll[k]);
    }
    ycopy = memory_calloc64f(nwork);
    // Move the results back onto data
    for (k=0; k<data->nobs; k++)
    {
        lnewDt = false;
        lnewStartTime = false;
        // Loop through commands and check if dt or npts changed
        ierr = sacio_getEpochalStartTime(data->obs[k].header, &epoch0);
        if (ierr != 0)
        {
            log_errorF("%s: Failed to get start time of trace\n", fcnm);
            goto ERROR;
        }
        epoch = epoch0;
        sacio_getFloatHeader(SAC_FLOAT_DELTA, data->obs[k].header, &dt0);
        sacio_getIntegerHeader(SAC_INT_NPTS, data->obs[k].header, &npts0);
        dt = dt0;
        for (i=0; i<commands->ncmds; i++)
        {
            if (commands[k].commands[i].type == CUT_COMMAND)
            {
                i0 = commands[k].commands[i].cut.i0;
                epoch = epoch + (double) i0*dt;
                lnewStartTime = true;
            }
            if (commands[k].commands[i].type == DOWNSAMPLE_COMMAND)
            {
                nq = commands[k].commands[i].downsample.nq;
                dt = dt*(double) nq;
                lnewDt = true;
            }
            if (commands[k].commands[i].type == DECIMATE_COMMAND)
            {
                nq = commands[k].commands[i].decimate.nqAll;
                dt = dt*(double) nq;
                lnewDt = true;
            }
        }
        // Extract the data - here a resize is required
        if (nyAll[k] != npts0)
        {
            ierr = process_getSerialCommandsData64f(commands[k], nwork,
                                                    &nyAll[k], ycopy);
            if (ierr != 0)
            {
                log_errorF("%s: Error extracting data onto ycopy\n", fcnm);
                goto ERROR;
            }
            sacio_freeData(&data->obs[k]);
            if (nyAll[k] > 0)
            {
                sacio_freeData(&data->obs[k]);
                data->obs[k].data = sacio_malloc64f(nyAll[k]);
                data->obs[k].npts = nyAll[k];
                sacio_setIntegerHeader(SAC_INT_NPTS, nyAll[k],
                                       &data->obs[k].header);
                ierr = array_copy64f_work(nyAll[k], ycopy, data->obs[k].data);
            }
        }
        // Otherwise just copy back onto the data structure
        else
        {
            ierr = process_getSerialCommandsData64f(commands[k], nwork,
                                                    &nyAll[k],
                                                    data->obs[k].data);
            if (ierr != 0)
            {
                log_errorF("%s: Error extracting data\n", fcnm);
                goto ERROR;
            }
        }
        // Update the times
        if (lnewStartTime)
        {
            // Update the picks
            for (i=0; i<nTimeVars; i++)
            {
                ierr = sacio_getFloatHeader(timeVars[i], data->obs[k].header,
                                            &time);
                if (ierr == 0)
                {
                    time = time + epoch0; // Turn to real time
                    time = time - epoch;  // Make relative to new time 
                    sacio_setFloatHeader(timeVars[i], time,
                                         &data->obs[k].header);
                }
                ierr = 0;
            }
            sacio_setEpochalStartTime(epoch, &data->obs[k].header);
        }
        // Update the sampling period
        if (lnewDt)
        {
            sacio_setFloatHeader(SAC_FLOAT_DELTA, dt, &data->obs[k].header);
        }
    }
    // Release memory
ERROR:;
    for (k=0; k<data->nobs; k++)
    {
        process_freeSerialCommands(&commands[k]);
    }
    memory_free64f(&ycopy);
    memory_free32i(&nyAll);
    free(commands);
    return ierr;
}
