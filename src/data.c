#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iniparser.h>
#include "tdsearch_data.h"
#include "ispl/process.h"
#include "sacio.h"
#include "ttimes.h"
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
    int ierr, item, k, ncmds, ncmdsWork, nitems, nfiles, nlines, nobs;
    bool luseDataList;
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
    int i, ierr, k;
    if (data->maxobs > 0 && data->obs != NULL)
    {
        for (k=0; k<data->nobs; k++)
        {
            ierr = sacio_free(&data->obs[k]);
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
                    free(data->cmds[k].cmds[i]);
                }
                free(data->cmds[k].cmds);
            }
        }
        free(data->cmds);
    }
    memset(data, 0, sizeof(struct tdSearchData_struct));
    return 0;
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
 * @param[out] ttimes   Epochal theoretical primary arrival times (UTC-seconds)
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
 */
char **tdsearch_data_modifyProcessingCommands(
    const int ncmds, const char **cmds,
    const double cut0, const double cut1, const double targetDt,
    const struct sacData_struct data,
    int *ierr)
{
    const char *fcnm = "tdsearch_data_modifyProcessingCommands\0";
    char **newCmds, cwork[MAX_CMD_LEN], cmd1[64], cmd2[64];
    double dt0;
    size_t lenos;
    int i;
    const bool laaFilter = true; // anti-alias filter in decimation
    const bool lfixPhase = true; // don't let anti-alias filter mess up picks
    *ierr = 0;
    newCmds = NULL;
    if (ncmds < 1){return newCmds;}
    // Modify the problematic commands
    *ierr = sacio_getFloatHeader(SAC_FLOAT_DELTA, data.header, &dt0);
    newCmds = (char **) calloc((size_t) ncmds, sizeof(char *));
    for (i=0; i<ncmds; i++)
    {
        lenos = strlen(cmds[i]);
        memset(cwork, 0, MAX_CMD_LEN*sizeof(char));
        memset(cmd1, 0, 64*sizeof(char));
        memset(cmd2, 0, 64*sizeof(char));
        if (strcasecmp(cmds[i], "transfer\0") == 0)
        {

        }
        if (strcasecmp(cmds[i], "decimate\0") == 0)
        {
            *ierr = decimate_createDesignCommandsWithDT(dt0, targetDt,
                                                        laaFilter, lfixPhase,
                                                        cwork);
            if (*ierr != 0)
            {
                log_errorF("%s: Couldn't modify the decimate command\n", fcnm);
                break;
            }
            dt0 = targetDt;
        }
        if (strcasecmp(cmds[i], "downsample\0") == 0)
        {
            *ierr = downsample_downsampleTargetSamplingPeriodToString(
                       dt0, targetDt, cwork);
            if (*ierr != 0)
            {
                log_errorF("%s: Couldn't modify downsample command\n", fcnm);
                break;
            }
        }
        else
        {
            newCmds[i] = (char *) calloc(lenos+1, sizeof(char));
            strcpy(newCmds[i], cmds[i]);
        }
    }
    return newCmds;
}
//============================================================================//
/*!
 * @brief Sets the pick time on the SAC header.
 *
 * @param[in] pickTime        Epochal time of phase arrival (UTC-seconds).
 * @param[in] phaseName       Name of phase.
 * @param[in] pickHeaderTime  Header variable identifier to hold the pick
 *                            time.  Likely SAC_float_A for a primary arrival,
 *                            but also may be SAC_FLOAT_T0-SAC_FLOAT_T9
 *                            for other phases.
 * @param[in] pickHeaderName  Header variable identifer to hold the phase name.
 *                            Likely SAC_CHAR_KA for a primary arrival,
 *                            but also may be SAC_CHAR_KT0-SAC_CHAR_KT9
 *                            for other phases.
 * 
 * @param[out] data           SAC data structure on which the phase pick was
 *                            set.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 */
/*
int tdsearch_data_setPickTime(const double pickTime,
                              const char *phaseName,
                              const enum sacHeader_enum pickHeaderTime,
                              const enum sacHeader_enum pickHeaderName,
                              struct sacData_struct *data)
{
    const char *fcnm = "tdsearch_data_setPickTime\0";
    double dt, epoch, p;
    int ierr, npts;
    // Get the trace epochal start time 
    ierr = sacio_getEpochalStartTime(data->header, &epoch);
    if (ierr != 0)
    {
        log_errorF("%s: Failed to get start time\n", fcnm);
        return -1; 
    }
    sacio_getIntegerHeader(SAC_INT_NPTS, data->header, &npts);
    sacio_getFloatHeader(SAC_FLOAT_DELTA, data->header, &dt);
    if (pickTime < epoch)
    {
        log_warnF("%s: Pick time is before trace start time\n", fcnm);
    }
    if (pickTime > epoch + (double) (npts - 1)*dt)
    {
        log_warnF("%s: Pick time comes in after trace end time\n", fcnm);
    }
    // Make it a relative time
    p = pickTime - epoch; 
    sacio_setFloatHeader(pickHeaderTime, p, &data->header); 
    sacio_setCharacterHeader(pickHeaderName, phaseName, &data->header);
    return 0;
}
*/
//============================================================================//
/*
int tdsearch_data_setObservation(const int iobs,
                                 const struct sacData_struct obs,
                                 struct tdSearchData_struct *data) 
{
    const char *fcnm = "tdsearch_data_setObservation\0";
    if (iobs > data->maxobs)
    {
        log_errorF("%s: Insufficient space; iobs > data->maxobs %d %d\n",
                  fcnm, iobs, data->maxobs);
        return -1;
    }
    sacio_copy(obs, &data->obs[iobs]);
    return 0;
}
*/
