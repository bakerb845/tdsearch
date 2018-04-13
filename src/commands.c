#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "tdsearch_commands.h"
#include "sacio.h"
#include "ispl/process.h"
#include "iscl/array/array.h"
#include "iscl/memory/memory.h"
#include "iscl/string/string.h"

char **tdsearch_commands_modifyCommands(
    const int ncmds, const char **cmds,
    const struct tdSearchModifyCommands_struct options,
    const struct sacData_struct data,
    int *ierr)
{
    char **newCmds, **cmdSplit, cwork[MAX_CMD_LEN], c64[64], cmd1[64], cmd2[64];
    double *freqs, cut0, cut1, dt0, epoch, ptime, t0, t1, targetDt;
    struct signalZPK_struct zpkFrom, zpkTo;
    int i, iodva, l, npolesAdd, nsplit, nzerosAdd;
    size_t lenos;
    bool ldeconvolution, oneCorner;
    const bool lisDigital = false; // SAC-PZ's aren't digital
    const bool laaFilter = true; // anti-alias filter in decimation
    const bool lfixPhase = true; // don't let anti-alias filter mess up picks
    const bool lpad2 = false; // don't pad signals
    *ierr = 0;
    freqs = NULL;
    cut0 = options.cut0;
    cut1 = options.cut1;
    targetDt = options.targetDt;
    ldeconvolution = options.ldeconvolution;
    iodva = options.iodva;
    if (iodva < 0 || iodva > 3)
    {
        fprintf(stdout, "%s: No idea what iodva=%d means\n", __func__, iodva);
    }
    newCmds = NULL;
    if (ncmds < 1){return 0;} // Nothing to do
    // Modify the problematic commands
    *ierr = sacio_getFloatHeader(SAC_FLOAT_DELTA, data.header, &dt0);
    newCmds = (char **) calloc((size_t) ncmds, sizeof(char *));
    for (i=0; i<ncmds; i++)
    {
        newCmds[i] = (char *) calloc(MAX_CMD_LEN, sizeof(char));
        lenos = strlen(cmds[i]);
        memset(cwork, 0, MAX_CMD_LEN*sizeof(char));
        memset(cmd1, 0, 64*sizeof(char));
        memset(cmd2, 0, 64*sizeof(char));
        strncpy(cmd1, cmds[i], MIN(3, lenos));
        strncpy(cmd2, cmds[i], MIN(3, lenos));
        npolesAdd = 0;
        nzerosAdd = 0;
        if (strcasecmp(cmds[i], "transfer\0") == 0)
        {
            memset(&zpkTo, 0, sizeof(struct signalZPK_struct));
            memset(&zpkFrom, 0, sizeof(struct signalZPK_struct));
            if (data.pz.lhavePZ)
            {
                // 1/From(w) -> need to eliminate a zero from To(w)
                // by adding a zero
                if (data.pz.inputUnits == SAC_METERS ||
                    data.pz.inputUnits == SAC_NANOMETERS)
                {
                    // Displacement to displacement
                    if (iodva == 0)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 0;
                    }
                    // Velocity to displacement
                    else if (iodva == 1)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 1;
                    }
                    // Acceleration to displacement
                    else if (iodva == 2)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 2;
                    }
                }
                // 1/From(w) = 1 -> do nothing
                else if (data.pz.inputUnits == SAC_METERS_SECOND ||
                         data.pz.inputUnits == SAC_NANOMETERS_SECOND)
                {
                    // Displacement to velocity
                    if (iodva == 0)
                    {
                        npolesAdd = 1;
                        nzerosAdd = 0;
                    }
                    // Velocity to velocity
                    else if (iodva == 1)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 0;
                    }
                    // Acceleration to velocity
                    else if (iodva == 2)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 1;
                    }
                }
                // 1/From(w) -> add a zero to To(w) by adding a pole
                else if (data.pz.inputUnits == SAC_METERS_SECOND_SECOND ||
                         data.pz.inputUnits == SAC_NANOMETERS_SECOND_SECOND)
                {
                    // Displacement to acceleration
                    if (iodva == 0)
                    {
                        npolesAdd = 2;
                        nzerosAdd = 0;
                    }
                    // Velocity to acceleration
                    else if (iodva == 1)
                    {
                        npolesAdd = 1;
                        nzerosAdd = 0;
                    }
                    // Acceleration to acceleration
                    else if (iodva == 2)
                    {
                        npolesAdd = 0;
                        nzerosAdd = 0;
                    }
                }
                else if (data.pz.inputUnits == SAC_UNKNOWN_UNITS)
                {
                    fprintf(stderr, "%s: Unknown input units\n", __func__);
                    *ierr = 1;
                    break;
                }
            }
            // Set the poles and zeros
            zpkFrom.npoles = npolesAdd;
            zpkFrom.nzeros = nzerosAdd; 
            zpkFrom.k = 1.0;
            if (zpkFrom.npoles > 0)
            {
                zpkFrom.p = memory_calloc64z(zpkFrom.npoles);
            }
            if (zpkFrom.nzeros > 0)
            {
                zpkFrom.z = memory_calloc64z(zpkFrom.nzeros);
            }
            zpkTo.npoles = data.pz.npoles;
            zpkTo.nzeros = data.pz.nzeros;
            zpkTo.k = data.pz.constant;
            if (zpkTo.npoles > 0)
            {
                zpkTo.p = array_copy64z(zpkTo.npoles, data.pz.poles, ierr);
            }
            if (zpkTo.nzeros > 0)
            {
                zpkTo.z = array_copy64z(zpkTo.nzeros, data.pz.zeros, ierr);
            }
            if (!ldeconvolution)
            {
                *ierr = transfer_poleZeroInstrumentResponsesToString(dt0,
                                                             freqs,
                                                             lisDigital,
                                                             true, zpkFrom,
                                                             true, zpkTo,
                                                             lpad2, cwork);
            }
            else
            {
                *ierr = transfer_poleZeroInstrumentResponsesToString(dt0,
                                                             freqs,
                                                             lisDigital,
                                                             true, zpkTo,
                                                             true, zpkFrom,
                                                             lpad2, cwork);
            }
            memory_free64z(&zpkTo.p);
            memory_free64z(&zpkTo.z);
            memory_free64z(&zpkFrom.p);
            memory_free64z(&zpkFrom.z);
            if (*ierr != 0)
            {
                fprintf(stderr, "%s: Error setting transfer command\n",
                        __func__);
                goto ERROR;
            }
/*
transfer_poleZeroInstrumentResponsesToString(
    dt0,
    const double *freqs,
    const bool lisDigital,
    const bool lhaveFrom, const struct signalZPK_struct zpkFrom,
    const bool lhaveTo, const  struct signalZPK_struct zpkTo,
    cwork);
*/
        }
        else if (strcasecmp(cmd1, "lowpass\0") == 0 ||
                 strcasecmp(cmd1, "highpas\0") == 0)
        {
            cmdSplit = string_rsplit(NULL, cmds[i], &nsplit);
            strcpy(cwork, cmdSplit[0]);
            strcat(cwork, " "); 
            for (l=1; l<nsplit; l++)
            {
                if (strcasecmp(cmdSplit[l], "corner\0") == 0)
                {
                    memset(c64, 0, 64*sizeof(char));
                    sprintf(c64, "dt %f corner %s", dt0, cmdSplit[l+1]);
                    strcat(cwork, c64);
                    l = l + 1;
                }
                else
                {
                    strcat(cwork, cmdSplit[l]);
                }
                if (l < nsplit - 1){strcat(cwork, " ");}
            }
            for (l=0; l<nsplit; l++){free(cmdSplit[l]);}
            free(cmdSplit);
        }
        else if (strcasecmp(cmd1, "bandpas\0") == 0 ||
                 strcasecmp(cmd1, "bandrej\0") == 0)
        {
            cmdSplit = string_rsplit(NULL, cmds[i], &nsplit);
            strcpy(cwork, cmdSplit[0]);
            strcat(cwork, " ");
            for (l=1; l<nsplit; l++)
            {
                if (strcasecmp(cmdSplit[l], "corners\0") == 0)
                {
                    memset(c64, 0, 64*sizeof(char));
                    sprintf(c64, "dt %f corners %s %s",
                            dt0, cmdSplit[l+1], cmdSplit[l+2]);
                    strcat(cwork, c64);
                    l = l + 2;
                }
                else
                {
                    strcat(cwork, cmdSplit[l]);
                }
                if (l < nsplit - 1){strcat(cwork, " ");}
            }
            for (l=0; l<nsplit; l++){free(cmdSplit[l]);}
            free(cmdSplit);
        }
        else if (strcasecmp(cmd2, "sos\0") == 0)
        {
            oneCorner = false;
            cmdSplit = string_rsplit(NULL, cmds[i], &nsplit);
            strcpy(cwork, cmdSplit[0]);
            strcat(cwork, " ");
            for (l=1; l<nsplit; l++)
            {
                if (strcasecmp(cmdSplit[l], "lowpass\0") == 0 ||
                    strcasecmp(cmdSplit[l], "highpas\0") == 0)
                {
                    oneCorner = true;
                }
                if (strcasecmp(cmdSplit[l], "bandpas\0") == 0 ||
                    strcasecmp(cmdSplit[l], "bandrej\0") == 0)
                {
                    oneCorner = false;
                }
            }
            for (l=1; l<nsplit; l++)
            {
                if (strcasecmp(cmdSplit[l], "corner\0") == 0 ||
                    strcasecmp(cmdSplit[l], "corners\0") == 0)
                {
                    if (oneCorner)
                    {
                        memset(c64, 0, 64*sizeof(char));
                        sprintf(c64, "dt %f corner %s", dt0, cmdSplit[l+1]);
                        strcat(cwork, c64);
                        l = l + 1;
                    }
                    else
                    {
                        memset(c64, 0, 64*sizeof(char));
                        sprintf(c64, "dt %f corners %s %s",
                                dt0, cmdSplit[l+1], cmdSplit[l+2]);
                        strcat(cwork, c64);
                        l = l + 2;
                    }
                }
                else
                {
                    strcat(cwork, cmdSplit[l]);
                }
                if (l < nsplit - 1){strcat(cwork, " ");}
            }
            for (l=0; l<nsplit; l++){free(cmdSplit[l]);}
            free(cmdSplit);
        }
        else if (strcasecmp(cmds[i], "cut\0") == 0)
        {
            *ierr = sacio_getEpochalStartTime(data.header, &epoch);
            if (*ierr != 0)
            {
                fprintf(stderr, "%s: Failed to get start time\n", __func__);
                goto ERROR;
            }
            *ierr = sacio_getFloatHeader(SAC_FLOAT_A, data.header,
                                         &ptime);
            if (*ierr != 0)
            {
                fprintf(stderr, "%s: Failed to get pick time\n", __func__);
                goto ERROR;
            }
            t0 = epoch + ptime + cut0; // superfluous; epoch will be removed
            t1 = epoch + ptime + cut1; // superfluous; epoch will be removed
            *ierr = cut_cutEpochalTimesToString(dt0, epoch, t0, t1, cwork);
            if (*ierr != 0)
            {
                fprintf(stderr, "%s: Failed to modify cut command\n", __func__);
                goto ERROR;
            }
        }
        else if (strcasecmp(cmds[i], "decimate\0") == 0)
        {
            *ierr = decimate_createDesignCommandsWithDT(dt0, targetDt,
                                                        laaFilter, lfixPhase,
                                                        cwork);
            if (*ierr != 0)
            {
                fprintf(stderr, "%s: Couldn't modify the decimate command\n",
                        __func__);
                goto ERROR;
            }
            dt0 = targetDt;
        }
        else if (strcasecmp(cmds[i], "downsample\0") == 0)
        {
            *ierr = downsample_downsampleTargetSamplingPeriodToString(
                        dt0, targetDt, cwork);
            if (*ierr != 0)
            {
                fprintf(stderr, "%s: Couldn't modify downsample command\n",
                        __func__);
                goto ERROR;
            }
            dt0 = targetDt;
        }
        else
        {
            strcpy(cwork, cmds[i]);
        }
        // Update the command
        strcpy(newCmds[i], cwork);
    } // Loop on commands
ERROR:;
    return newCmds;
}
