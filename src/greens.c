#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iniparser.h>
#include "parmt_utils.h"
#include "tdsearch_greens.h"
#include "tdsearch_commands.h"
#ifdef TDSEARCH_USE_INTEL
#include <mkl_blas.h>
#else
#include <cblas.h>
#endif
#include "tdsearch_struct.h"
#include "tdsearch_hudson.h"
#include "ispl/process.h"
#include "iscl/array/array.h"
#include "iscl/fft/fft.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"

static int getPrimaryArrival(const struct sacHeader_struct hdr,
                             double *time, char phaseName[8]);

/*!
 * @brief Reads the generic Green's functions pre-processing files from
 *        the ini file.
 *
 * @param[in] iniFile    Name of ini file.
 * @param[in] nobs       Number of observations.
 *
 * @param[in,out] grns   On input contains the number of observations.
 *                       On exit contains the generic Green's functions
 *                       processing commands for the observations.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 * @TODO Add an option to read a processing list file.
 *
 */
int tdsearch_greens_setPreprocessingCommandsFromIniFile(
    const char *iniFile,
    const int nobs,  
    struct tdSearchGreens_struct *grns)
{
    dictionary *ini;
    char **cmds;
    const char *s;
    char varname[128];
    size_t lenos;
    int ierr, k, ncmds, ncmdsWork;
    ierr = 0;
    grns->nobs = nobs;
    if (grns->nobs < 1){return 0;}
    if (!os_path_isfile(iniFile))
    {
        fprintf(stderr, "%s: Error ini file %s doesn't exist\n",
                __func__, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    ncmds = iniparser_getint(ini, "tdSearch:greens:nCommands\0", 0);
    grns->cmds = (struct tdSearchDataProcessingCommands_struct *)
                 calloc((size_t) nobs, //ncmds,
                        sizeof(struct tdSearchDataProcessingCommands_struct));
    //if (!luseProcessingList)
    {
        if (ncmds > 0)
        {
            ncmdsWork = ncmds;
            ncmds = 0;
            cmds = (char **) calloc((size_t) ncmdsWork, sizeof(char *));
            for (k=0; k<ncmdsWork; k++)
            {
                memset(varname, 0, 128*sizeof(char));
                sprintf(varname, "tdSearch:greens:command_%d", k+1);
                s = iniparser_getstring(ini, varname, NULL);
                if (s == NULL){continue;}
                lenos = strlen(s);
                cmds[ncmds] = (char *) calloc(lenos+1, sizeof(char));
                strcpy(cmds[ncmds], s);
                //printf("%s\n", cmds[ncmds]);
                ncmds = ncmds + 1;
            }
            // Attach these to the data processing structure
            for (k=0; k<grns->nobs; k++)
            {
                ierr = tdsearch_greens_attachCommandsToGreens(
                           k, ncmds, (const char **) cmds, grns);
            }
            // Free space
            for (k=0; k<ncmdsWork; k++)
            {
                if (cmds[k] != NULL){free(cmds[k]);}
            }
            free(cmds);
        }
    }
    iniparser_freedict(ini);
    return ierr;
}
//============================================================================//
/*!
 * @brief Attaches the Green's functions processing commands to the Green's
 *        structure.
 * 
 * @param[in] iobs      Observation number to attach commands to.
 * @param[in] ncmds     Number of commands.
 * @param[in] cmds      Pre-processing commands for the Green's functions
 *                      corresponding to the iobs'th observation [ncmds].
 *
 * @param[in,out] grns  On input has the Green's functions (depths and t*'s)
 *                      corresponding to the observations.
 *                      On output contains the pre-processing commands
 *                      for Green's functions corresponding to the
 *                      observations.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_attachCommandsToGreens(const int iobs, const int ncmds,
                                           const char **cmds,
                                           struct tdSearchGreens_struct *grns)
{
    int i;
    size_t lenos;
    // Make sure iobs is in bounds 
    if (iobs < 0 || iobs >= grns->nobs)
    {    
        fprintf(stderr, "%s: Error iobs=%d is out of bounds [0,%d]\n",
                __func__, iobs, grns->nobs);
        return -1;      
    }        
    // Try to handle space allocation if not already done
    if (grns->cmds == NULL && grns->nobs > 0)
    {    
        grns->cmds = (struct tdSearchDataProcessingCommands_struct *)
                     calloc((size_t) grns->nobs,
                        sizeof(struct tdSearchDataProcessingCommands_struct));
    }    
    if (grns->cmds == NULL)
    {    
        fprintf(stderr, "%s: Error grns->cmds is NULL\n", __func__);
        return -1;
    }
    if (ncmds == 0){return 0;}
    grns->cmds[iobs].ncmds = ncmds;
    grns->cmds[iobs].cmds = (char **) calloc((size_t) ncmds, sizeof(char *));
    if (cmds == NULL){printf("problem 1\n");}
    for (i=0; i<ncmds; i++)
    {
        if (cmds[i] == NULL){printf("problem 2\n");}
        lenos = strlen(cmds[i]);
        grns->cmds[iobs].cmds[i] = (char *) calloc(lenos+1, sizeof(char));
        strcpy(grns->cmds[iobs].cmds[i], cmds[i]);
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Convenience function which returns the index of the Green's
 *        function on the Green's function structure.
 *
 * @param[in] GMT_TERM    Name of the desired Green's function:
 *                        (G11_TERM, G22_TERM, ..., G23_TERM).
 * @param[in] iobs        Desired observation number (C numbering).
 * @param[in] itstar      Desired t* (C numbering).
 * @param[in] idepth      Desired depth (C numbering).
 * @param[in] grns        Contains the number of observations, depths, and t*'s
 *                        on the Green's functions structure.
 *
 * @result Negative indicates failure.  Otherwise, this is the index in 
 *         grns.grns corresponding to the desired 
 *         (iobs, idepth, itstar, G??_GRNS) coordinate.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_getGreensFunctionIndex(
    const enum prepmtGreens_enum GMT_TERM,
    const int iobs, const int itstar, const int idepth,
    const struct tdSearchGreens_struct grns)
{
    int igx, indx;
    indx =-1;
    igx = (int) GMT_TERM - 1;
    if (igx < 0 || igx > 5)
    {
        fprintf(stderr, "%s: Can't classify Green's functions index\n",
                __func__);
        return indx;
    }
    indx = iobs*(6*grns.ntstar*grns.ndepth)
         + idepth*(6*grns.ntstar)
         + itstar*6
         + igx;
    if (indx >= grns.ngrns)
    {
        fprintf(stdout, "%s: indx out of bounds - segfault is coming\n",
                __func__);
        return -1;
    }
    return indx;
}
//============================================================================//
/*!
 * @brief Convenience function for extracting the: 
 *        \$ \{ G_{xx}, G_{yy}, G_{zz}, G_{xy}, G_{xz}, G_{yz} \} \$
 *        Green's functions indices for the observation, t*, and depth.
 *
 * @param[in] iobs      Observation number.
 * @param[in] itstar    t* index.  This is C numbered.
 * @param[in] idepth    Depth index.  This is C numbered.
 * @param[in] grns      Contains the Green's functions.
 * @param[out] indices  Contains the Green's functions indices defining
 *                      the indices that return the:
 *                       \$ \{ G_{xx}, G_{yy}, G_{zz}, 
 *                             G_{xy}, G_{xz}, G_{yz} \} \$
 *                      for this observation, t*, and depth. 
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_getGreensFunctionsIndices(
    const int iobs, const int itstar, const int idepth,
    const struct tdSearchGreens_struct grns, int indices[6])
{
    int i, ierr;
    const enum prepmtGreens_enum mtTerm[6] = 
       {G11_GRNS, G22_GRNS, G33_GRNS, G12_GRNS, G13_GRNS, G23_GRNS};
    ierr = 0;
    for (i=0; i<6; i++)
    {
        indices[i] = tdsearch_greens_getGreensFunctionIndex(mtTerm[i],
                                                          iobs, itstar, idepth,
                                                          grns); 
        if (indices[i] < 0){ierr = ierr + 1;}
    } 
    return ierr;
}
//============================================================================//
/*!
 * @brief Releases memory on the Greens functions structure.
 *
 * @param[out] grns   On exit all memory has been freed and variables set to
 *                    0 or NULL.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_free(struct tdSearchGreens_struct *grns)
{
    int i, k;
    if (grns->grns != NULL && grns->ngrns > 0)
    {
        for (k=0; k<grns->ngrns; k++)
        {
            sacio_free(&grns->grns[k]);
        }
        free(grns->grns);
    }
    if (grns->cmds != NULL)
    {
        for (k=0; k<grns->nobs; k++)
        {
            if (grns->cmds[k].cmds != NULL)
            {
                for (i=0; i<grns->cmds[k].ncmds; i++)
                {
                    if (grns->cmds[k].cmds[i] != NULL)
                    {
                        free(grns->cmds[k].cmds[i]);
                        grns->cmds[k].cmds[i] = NULL;
                    }
                }
                free(grns->cmds[k].cmds);
                grns->cmds[k].cmds = NULL;
            }
        }
        free(grns->cmds);
        grns->cmds = NULL;
    }
    if (grns->cmdsGrns != NULL)
    {
        for (k=0; k<grns->ngrns; k++)
        {
            if (grns->cmdsGrns[k].cmds != NULL)
            {
                for (i=0; i<grns->cmdsGrns[k].ncmds; i++)
                {
                    if (grns->cmdsGrns[k].cmds != NULL)
                    {
                        free(grns->cmdsGrns[k].cmds[i]);
                        grns->cmdsGrns[k].cmds[i] = NULL;
                    }
                }
                free(grns->cmdsGrns[k].cmds);
                grns->cmdsGrns[k].cmds = NULL;
             }
        } 
        free(grns->cmdsGrns);
        grns->cmdsGrns = NULL;
    }
    memset(grns, 0, sizeof(struct tdSearchGreens_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Converts the fundamental faults Green's functions to Green's
 *        functions that can be used by tdsearch.
 *
 * @param[in] data    tdSearch data structure.
 * @param[in] ffGrns  fundamental fault Green's functions for every t* and
 *                    depth in the grid search for each observation. 
 *
 * @param[out] grns   Contains the Green's functions that can be applied to
 *                    a moment tensor to produce a synthetic for every t*
 *                    and depth in the grid search for each observation.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_ffGreensToGreens(const struct tdSearchData_struct data,
                                     const struct tdSearchHudson_struct ffGrns,
                                     struct tdSearchGreens_struct *grns)
{
    char knetwk[8], kstnm[8], kcmpnm[8], khole[8], phaseName[8],
         phaseNameGrns[8];
    double az, baz, cmpaz, cmpinc, cmpincSEED, dt0, epoch, epochNew,
           evla, evlo, o, pick, pickTime, pickTimeGrns, stel, stla, stlo;
    int i, icomp, id, ierr, idx, indx, iobs, it, kndx, l, npts;
    const char *kcmpnms[6] = {"GXX\0", "GYY\0", "GZZ\0",
                              "GXY\0", "GXZ\0", "GYZ\0"};
    const double xmom = 1.0;     // no confusing `relative' magnitudes 
    const double xcps = 1.e-20;  // convert dyne-cm mt to output cm
    const double cm2m = 1.e-2;   // cm to meters
    const double dcm2nm = 1.e+7; // magnitudes intended to be specified in
                                 // Dyne-cm but I work in N-m
    // Given a M0 in Newton-meters get a seismogram in meters
    const double xscal = xmom*xcps*cm2m*dcm2nm;
    const int nTimeVars = 11; 
    const enum sacHeader_enum pickVars[11]
       = {SAC_FLOAT_A,
          SAC_FLOAT_T0, SAC_FLOAT_T1, SAC_FLOAT_T2, SAC_FLOAT_T3,
          SAC_FLOAT_T4, SAC_FLOAT_T5, SAC_FLOAT_T6, SAC_FLOAT_T7,
          SAC_FLOAT_T8, SAC_FLOAT_T9};
    //memset(grns, 0, sizeof(struct tdSearchGreens_struct));
    grns->ntstar = ffGrns.ntstar;
    grns->ndepth = ffGrns.ndepth;
    grns->nobs = data.nobs;
    grns->ngrns = 6*grns->ntstar*grns->ndepth*grns->nobs;
    if (grns->ngrns < 1)
    {
        fprintf(stderr, "%s: Error grns is empty\n", __func__);
        return -1;
    }
    grns->grns = (struct sacData_struct *)
                 calloc((size_t) grns->ngrns, sizeof(struct sacData_struct));
    for (iobs=0; iobs<data.nobs; iobs++)
    {
        ierr = 0;
        ierr += sacio_getFloatHeader(SAC_FLOAT_AZ,
                                     data.obs[iobs].header, &az);
        ierr += sacio_getFloatHeader(SAC_FLOAT_BAZ,
                                     data.obs[iobs].header, &baz);
        ierr += sacio_getFloatHeader(SAC_FLOAT_CMPINC,
                                     data.obs[iobs].header, &cmpinc);
        ierr += sacio_getFloatHeader(SAC_FLOAT_CMPAZ,
                                     data.obs[iobs].header, &cmpaz);
        ierr += sacio_getFloatHeader(SAC_FLOAT_EVLA,
                                     data.obs[iobs].header, &evla);
        ierr += sacio_getFloatHeader(SAC_FLOAT_EVLO,
                                     data.obs[iobs].header, &evlo);
        ierr += sacio_getFloatHeader(SAC_FLOAT_STLA,
                                     data.obs[iobs].header, &stla);
        ierr += sacio_getFloatHeader(SAC_FLOAT_STLO,
                                     data.obs[iobs].header, &stlo);
        ierr += sacio_getFloatHeader(SAC_FLOAT_DELTA,
                                     data.obs[iobs].header, &dt0);
        ierr += sacio_getCharacterHeader(SAC_CHAR_KNETWK,
                                         data.obs[iobs].header, knetwk);
        ierr += sacio_getCharacterHeader(SAC_CHAR_KSTNM,
                                         data.obs[iobs].header, kstnm);
        ierr += sacio_getCharacterHeader(SAC_CHAR_KCMPNM,
                                         data.obs[iobs].header, kcmpnm);
        // this one isn't critical
        sacio_getFloatHeader(SAC_FLOAT_STEL,
                             data.obs[iobs].header, &stel);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error reading header variables\n", __func__);
            break;
        }
        // Station location code is not terribly important
        sacio_getCharacterHeader(SAC_CHAR_KHOLE, data.obs[iobs].header, khole);
        cmpincSEED = cmpinc - 90.0; // SAC to SEED convention
        // Get the primary arrival
        ierr = sacio_getEpochalStartTime(data.obs[iobs].header, &epoch);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error getting start time\n", __func__);
            break;
        }
        ierr += getPrimaryArrival(data.obs[iobs].header, &pickTime, phaseName);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error getting primary pick\n", __func__);
            break;
        }
        // Need to figure out the component
        icomp = 1;
        if (kcmpnm[2] == 'Z' || kcmpnm[2] == 'z' || kcmpnm[2] == '1')
        {
            icomp = 1;
        }
        else if (kcmpnm[2] == 'N' || kcmpnm[2] == 'n' || kcmpnm[2] == '2')
        {
            icomp = 2;
        }
        else if (kcmpnm[2] == 'E' || kcmpnm[2] == 'e' || kcmpnm[2] == '3')
        {
            icomp = 3;
        }
        else
        {
            fprintf(stderr, "%s: Can't classify component: %s\n",
                    __func__, kcmpnm);
        }
        // Process all Green's functions in this block
        for (id=0; id<ffGrns.ndepth; id++)
        {
            for (it=0; it<ffGrns.ntstar; it++)
            {
                idx = tdsearch_hudson_observationDepthTstarToIndex(iobs, id, it,
                                                                   ffGrns);
                kndx = 10*idx;
                sacio_getFloatHeader(SAC_FLOAT_O, ffGrns.grns[kndx].header, &o);
                getPrimaryArrival(ffGrns.grns[kndx].header,
                                  &pickTimeGrns, phaseNameGrns);
                if (strcasecmp(phaseNameGrns, phaseName) != 0)
                {
                    fprintf(stdout, "%s: Phase name mismatch %s %s\n",
                            __func__, phaseName, phaseNameGrns);
                }
                npts = ffGrns.grns[kndx].npts;
                indx = tdsearch_greens_getGreensFunctionIndex(G11_GRNS,
                                                              iobs, it, id,
                                                              *grns);
                for (i=0; i<6; i++)
                {
                    sacio_copy(ffGrns.grns[kndx], &grns->grns[indx+i]);
                    if (data.obs[iobs].pz.lhavePZ)
                    {
                        sacio_copyPolesAndZeros(data.obs[iobs].pz,
                                                &grns->grns[indx+i].pz);
                    }
                    sacio_setFloatHeader(SAC_FLOAT_AZ, az,
                                         &grns->grns[indx+i].header);
                    sacio_setFloatHeader(SAC_FLOAT_BAZ, baz,
                                         &grns->grns[indx+i].header);
                    sacio_setFloatHeader(SAC_FLOAT_CMPAZ, cmpaz,
                                         &grns->grns[indx+i].header);
                    sacio_setFloatHeader(SAC_FLOAT_CMPINC, cmpinc,
                                         &grns->grns[indx+i].header);
                    sacio_setFloatHeader(SAC_FLOAT_EVLA, evla,
                                         &grns->grns[indx+i].header);
                    sacio_setFloatHeader(SAC_FLOAT_EVLO, evlo,
                                         &grns->grns[indx+i].header);
                    sacio_setFloatHeader(SAC_FLOAT_STLA, stla,
                                         &grns->grns[indx+i].header);
                    sacio_setFloatHeader(SAC_FLOAT_STLO, stlo,
                                         &grns->grns[indx+i].header);
                    sacio_setFloatHeader(SAC_FLOAT_STEL, stel,
                                         &grns->grns[indx+i].header); 
                    sacio_setCharacterHeader(SAC_CHAR_KNETWK, knetwk,
                                             &grns->grns[indx+i].header);
                    sacio_setCharacterHeader(SAC_CHAR_KSTNM, kstnm,
                                             &grns->grns[indx+i].header);
                    sacio_setCharacterHeader(SAC_CHAR_KHOLE, khole,
                                             &grns->grns[indx+i].header);
                    sacio_setCharacterHeader(SAC_CHAR_KCMPNM, kcmpnms[i],
                                             &grns->grns[indx+i].header);
                    sacio_setCharacterHeader(SAC_CHAR_KEVNM, "SYNTHETIC\0",
                                             &grns->grns[indx+i].header);
                    // Set the start time by aligning on the arrival
                    epochNew = epoch + (pickTime - o) - pickTimeGrns;
                    sacio_setEpochalStartTime(epochNew,
                                              &grns->grns[indx+i].header);
                    // Update the pick times
                    for (l=0; l<11; l++)
                    {
                        ierr = sacio_getFloatHeader(pickVars[l],
                                                    grns->grns[indx+i].header,
                                                    &pick);
                        if (ierr == 0)
                        {
                            pick = pick + o;
                            sacio_setFloatHeader(pickVars[l],
                                                 pick,
                                                 &grns->grns[indx+i].header);
                        }
                    } 
                }
                ierr = parmt_utils_ff2mtGreens64f(npts, icomp,
                                                  az, baz,
                                                  cmpaz, cmpincSEED,
                                                  ffGrns.grns[kndx+0].data,
                                                  ffGrns.grns[kndx+1].data,
                                                  ffGrns.grns[kndx+2].data,
                                                  ffGrns.grns[kndx+3].data,
                                                  ffGrns.grns[kndx+4].data,
                                                  ffGrns.grns[kndx+5].data,
                                                  ffGrns.grns[kndx+6].data,
                                                  ffGrns.grns[kndx+7].data,
                                                  ffGrns.grns[kndx+8].data,
                                                  ffGrns.grns[kndx+9].data,
                                                  grns->grns[indx+0].data,
                                                  grns->grns[indx+1].data,
                                                  grns->grns[indx+2].data,
                                                  grns->grns[indx+3].data,
                                                  grns->grns[indx+4].data,
                                                  grns->grns[indx+5].data);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Failed to rotate Greens functions\n",
                            __func__);
                }
                // Fix the characteristic magnitude scaling in CPS 
                for (i=0; i<6; i++)
                {
                    cblas_dscal(npts, xscal, grns->grns[indx+i].data, 1); 
                }
            }
        }
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Repicks the Green's functions onset time with an STA/LTA picker.
 *
 * @param[in] sta        Short term average window length (seconds).
 * @param[in] lta        Long term average window length (seconds).
 * @param[in] threshPct  Percentage of max STA/LTA after which an arrival
 *                       is declared.
 * @param[in] iobs       C numbered observation index.
 * @param[in] itstar     C numbered t* index.
 * @param[in] idepth     C numbered depth index.
 *
 * @param[in,out] grns   On input contains the Green's functions.
 *                       On output the first arrival time has been modified
 *                       with an STA/LTA picker.
 *
 * @brief 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_repickGreensWithSTALTA(
    const double sta, const double lta, const double threshPct,
    const int iobs, const int itstar, const int idepth,
    struct tdSearchGreens_struct *grns)
{
    struct stalta_struct stalta;
    double *charFn, *g, *Gxx, *Gyy, *Gzz, *Gxy, *Gxz, *Gyz,
           *gxxPad, *gyyPad, *gzzPad, *gxyPad, *gxzPad, *gyzPad,
           charMax, dt, tpick;
    int indices[6], ierr, k, npts, nlta, npad, nsta, nwork, prePad;
    // Check STA/LTA 
    ierr = 0;
    memset(&stalta, 0, sizeof(struct stalta_struct));
    if (lta < sta || sta < 0.0)
    {
        if (lta < sta){fprintf(stderr,"%s: Error lta < sta\n", __func__);}
        if (sta < 0.0){fprintf(stderr,"%s: Error sta is negative\n", __func__);}
        return -1;
    }
    ierr = tdsearch_greens_getGreensFunctionsIndices(iobs, itstar, idepth,
                                                     *grns, indices);
    if (ierr != 0)
    {
        fprintf(stderr, "%s: Failed to get Greens functions indicies\n",
                 __func__);
        return -1;
    }
    ierr = sacio_getIntegerHeader(SAC_INT_NPTS,
                                  grns->grns[indices[0]].header, &npts);
    if (ierr != 0 || npts < 1)
    {
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error getting number of points from header\n",
                    __func__);
        }
        else
        {
            fprintf(stderr, "%s: ERror no data points\n", __func__);
        }
        return -1;
    }
    ierr = sacio_getFloatHeader(SAC_FLOAT_DELTA,
                                grns->grns[indices[0]].header, &dt);
    if (ierr != 0 || dt <= 0.0)
    {
        if (ierr != 0){fprintf(stderr, "%s: failed to get dt\n", __func__);}
        if (dt <= 0.0){fprintf(stderr, "%s: invalid sampling period\n", __func__);}
        return -1;
    }
    // Define the windows
    nsta = (int) (sta/dt + 0.5);
    nlta = (int) (lta/dt + 0.5);
    prePad = MAX(64, fft_nextpow2(nlta, &ierr));
    npad = prePad + npts;
    // Set space
    gxxPad = memory_calloc64f(npad);
    gyyPad = memory_calloc64f(npad);
    gzzPad = memory_calloc64f(npad);
    gxyPad = memory_calloc64f(npad);
    gxzPad = memory_calloc64f(npad);
    gyzPad = memory_calloc64f(npad);
    charFn = memory_calloc64f(npad);
    // Reference pointers
    Gxx = grns->grns[indices[0]].data;
    Gyy = grns->grns[indices[1]].data;
    Gzz = grns->grns[indices[2]].data;
    Gxy = grns->grns[indices[3]].data;
    Gxz = grns->grns[indices[4]].data;
    Gyz = grns->grns[indices[5]].data;
    // Pre-pad signals
    array_set64f_work(prePad, Gxx[0], gxxPad);
    array_set64f_work(prePad, Gyy[0], gyyPad);
    array_set64f_work(prePad, Gzz[0], gzzPad); 
    array_set64f_work(prePad, Gxy[0], gxyPad);
    array_set64f_work(prePad, Gxz[0], gxzPad);
    array_set64f_work(prePad, Gyz[0], gyzPad);
    // Copy rest of array
    array_copy64f_work(npts, Gxx, &gxxPad[prePad]);
    array_copy64f_work(npts, Gyy, &gyyPad[prePad]);
    array_copy64f_work(npts, Gzz, &gzzPad[prePad]);
    array_copy64f_work(npts, Gxy, &gxyPad[prePad]);
    array_copy64f_work(npts, Gxz, &gxzPad[prePad]);
    array_copy64f_work(npts, Gyz, &gyzPad[prePad]);
    // apply the sta/lta
    for (k=0; k<6; k++)
    {
        g = NULL;
        if (k == 0)
        {
            g = gxxPad;
        }
        else if (k == 1)
        {
            g = gyyPad;
        }
        else if (k == 2)
        {
            g = gzzPad;
        }
        else if (k == 3)
        {
            g = gxyPad;
        }
        else if (k == 4)
        {
            g = gxzPad;
        }
        else if (k == 5)
        {
            g = gyzPad;
        }
        ierr = stalta_setShortAndLongTermAverage(nsta, nlta, &stalta);
        if (ierr != 0)
        {
            printf("%s: Error setting STA/LTA\n", __func__);
            break;
        }
        ierr = stalta_setData64f(npad, g, &stalta);
        if (ierr != 0)
        {
            printf("%s: Error setting data\n", __func__);
            break;
        }
        ierr = stalta_applySTALTA(&stalta);
        if (ierr != 0)
        {
            printf("%s: Error applying STA/LTA\n", __func__);
            break;
        }
        ierr = stalta_getData64f(stalta, npad, &nwork, g);
        if (ierr != 0)
        {
            printf("%s: Error getting result\n", __func__);
            break;
        }
        cblas_daxpy(npad, 1.0, g, 1, charFn, 1);
        stalta_resetInitialConditions(&stalta);
        stalta_resetFinalConditions(&stalta);
        g = NULL;
    }
    // Compute the pick time
    charMax = array_max64f(npts, &charFn[prePad], &ierr);
    tpick =-1.0;
    for (k=prePad; k<npad; k++)
    {
        if (charFn[k] > 0.01*threshPct*charMax)
        {
            tpick = (double) (k - prePad)*dt;
            break;
        }
    }
    if (tpick ==-1.0)
    {
        tpick = (double) (array_argmax64f(npad, charFn, &ierr) - prePad)*dt;
    }
    // Overwrite the pick; by this point should be on SAC_FLOAT_A
    //double apick;
    //sacio_getFloatHeader(SAC_FLOAT_A, grns->grns[indices[0]].header, &apick);
    for (k=0; k<6; k++)
    {
        sacio_setFloatHeader(SAC_FLOAT_A, tpick,
                             &grns->grns[indices[k]].header);
    }
    // Dereference pointers and free space
    Gxx = NULL;
    Gyy = NULL;
    Gzz = NULL;
    Gxy = NULL;
    Gxz = NULL;
    Gyz = NULL;
    memory_free64f(&gxxPad);
    memory_free64f(&gyyPad);
    memory_free64f(&gzzPad);
    memory_free64f(&gxyPad);
    memory_free64f(&gxzPad);
    memory_free64f(&gyzPad);
    memory_free64f(&charFn);
    stalta_free(&stalta);
    return ierr;
}
//============================================================================//
/*!
 * @brief Modifies the Green's functions processing commands.
 */
int tdsearch_greens_modifyProcessingCommands(
    const int iodva,
    const double cut0, const double cut1, const double targetDt,
    struct tdSearchGreens_struct *grns)
{
    struct tdSearchModifyCommands_struct options;
    const char **cmds;
    char **newCmds;
    size_t lenos;
    int i, ierr, iobs, k, kndx1, kndx2, ncmds;
    ierr = 0;
    grns->cmdsGrns = (struct tdSearchDataProcessingCommands_struct *)
                     calloc((size_t) grns->ngrns,
                        sizeof(struct tdSearchDataProcessingCommands_struct));
    for (iobs=0; iobs<grns->nobs; iobs++)
    {
        newCmds = NULL;
        ncmds =  grns->cmds[iobs].ncmds;
        if (ncmds < 1){continue;} // Nothing to do
        cmds = (const char **) grns->cmds[iobs].cmds;
        options.cut0 = cut0;
        options.cut1 = cut1;
        options.targetDt = targetDt;
        options.ldeconvolution = false;
        options.iodva = iodva;
        kndx1 = iobs*(6*grns->ntstar*grns->ndepth);
        kndx2 = (iobs+1)*(6*grns->ntstar*grns->ndepth);
        newCmds = tdsearch_commands_modifyCommands(ncmds, (const char **) cmds,
                                                   options,
                                                   grns->grns[kndx1], &ierr);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to set processing commands\n", __func__);
            goto ERROR;
        }
        // Expand the processing commands 
        for (k=kndx1; k<kndx2; k++)
        {
            grns->cmdsGrns[k].ncmds = ncmds;
            grns->cmdsGrns[k].cmds = (char **)
                                     calloc((size_t) ncmds, sizeof(char *));
            for (i=0; i<ncmds; i++)
            {
                lenos = strlen(newCmds[i]);
                grns->cmdsGrns[k].cmds[i] = (char *)
                                            calloc(lenos+1, sizeof(char));
                strcpy(grns->cmdsGrns[k].cmds[i], newCmds[i]);
            }
        } 
        // Release the memory
        if (newCmds != NULL)
        {
            for (i=0; i<ncmds; i++)
            {
                if (newCmds[i] != NULL){free(newCmds[i]);}
            }
            free(newCmds);
        }
    }
ERROR:;
    return ierr;
}
//============================================================================//
/*!
 * @brief Processes the Green's functions.
 *
 * @param[in,out] grns    On input contains the Green's functions for all
 *                        observations, depths, and t*'s as well as the
 *                        processing chains.
 *                        On output contains the filtered Green's functions
 *                        for all observations, depths, and t*'s.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_process(struct tdSearchGreens_struct *grns)
{
    double *data;
    struct serialCommands_struct commands;
    struct parallelCommands_struct parallelCommands;
    double dt, dt0, epoch, epoch0, time;
    int dataPtr[7], i, i0, ierr, iobs, idep, it, kndx, l, npts, npts0, nq,
        nwork, ny, ns;
    bool lnewDt, lnewStartTime;
    const int nsignals = 6;
    const int nTimeVars = 12;
    const enum sacHeader_enum timeVars[12]
       = {SAC_FLOAT_A, SAC_FLOAT_O, 
          SAC_FLOAT_T0, SAC_FLOAT_T1, SAC_FLOAT_T2, SAC_FLOAT_T3,
          SAC_FLOAT_T4, SAC_FLOAT_T5, SAC_FLOAT_T6, SAC_FLOAT_T7,
          SAC_FLOAT_T8, SAC_FLOAT_T9};
    // Loop on the observations
    for (iobs=0; iobs<grns->nobs; iobs++)
    {
        for (idep=0; idep<grns->ndepth; idep++)
        {
            for (it=0; it<grns->ntstar; it++)
            {
                memset(&parallelCommands, 0,
                       sizeof(struct parallelCommands_struct)); 
                memset(&commands, 0, sizeof(struct serialCommands_struct));
                kndx = tdsearch_greens_getGreensFunctionIndex(G11_GRNS,
                                                              iobs, it, idep,
                                                              *grns);
                ierr = process_stringsToSerialCommandsOptions(
                                      grns->cmdsGrns[kndx].ncmds,
                                      (const char **) grns->cmdsGrns[kndx].cmds,
                                      &commands);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Error setting serial command string\n", __func__);
                    goto ERROR;
                }
                // Determine some characteristics of the processing
                sacio_getEpochalStartTime(grns->grns[kndx].header, &epoch0);
                sacio_getFloatHeader(SAC_FLOAT_DELTA,
                                     grns->grns[kndx].header, &dt0);
                lnewDt = false;
                lnewStartTime = false;
                epoch = epoch0;
                dt = dt0;
                for (i=0; i<commands.ncmds; i++)
                {
                     if (commands.commands[i].type == CUT_COMMAND)
                     {
                         i0 = commands.commands[i].cut.i0;
                         epoch = epoch + (double) i0*dt;
                         lnewStartTime = true;
                     }
                     if (commands.commands[i].type == DOWNSAMPLE_COMMAND)
                     {
                         nq = commands.commands[i].downsample.nq;
                         dt = dt*(double) nq;
                         lnewDt = true;
                     }
                     if (commands.commands[i].type == DECIMATE_COMMAND)
                     {
                         nq = commands.commands[i].decimate.nqAll;
                         dt = dt*(double) nq;
                         lnewDt = true;
                     }
                }
                dataPtr[0] = 0;
                for (i=0; i<nsignals; i++)
                {
                    dataPtr[i+1] = dataPtr[i] + grns->grns[kndx+i].npts;
                }
                process_setCommandOnAllParallelCommands(nsignals, commands,
                                                        &parallelCommands);
                data = memory_calloc64f(dataPtr[nsignals]);
                for (i=0; i<nsignals; i++)
                {
                    ierr = array_copy64f_work(grns->grns[kndx+i].npts,
                                              grns->grns[kndx+i].data,
                                              &data[dataPtr[i]]);
                    if (ierr != 0)
                    {
                         fprintf(stderr, "%s: Failed to copy data\n", __func__);
                         goto ERROR;
                    }
                }
                ierr = process_setParallelCommandsData64f(nsignals, dataPtr,
                                                          data,
                                                          &parallelCommands);
                if (ierr != 0)
                {
                    fprintf(stderr, "%s: Failed to set data\n", __func__);
                    goto ERROR;
                }
                ierr = process_applyParallelCommands(&parallelCommands);
                if (ierr != 0)
                {
                    fprintf(stderr ,"%s: Failed to process data\n", __func__);
                    goto ERROR;
                }
                // Get the data
                nwork = dataPtr[nsignals];
                ierr = process_getParallelCommandsData64f(parallelCommands,
                                                          -1, nsignals,
                                                          &ny, &ns,
                                                          dataPtr, data);
                if (ny > nwork)
                {
                    memory_free64f(&data);
                    data = memory_calloc64f(ny);
                }
                nwork = ny;
                ierr = process_getParallelCommandsData64f(parallelCommands,
                                                          nwork, nsignals,
                                                          &ny, &ns,
                                                          dataPtr, data);
                for (i=0; i<nsignals; i++)
                {
                    sacio_getIntegerHeader(SAC_INT_NPTS,  
                                           grns->grns[kndx+i].header, &npts0);
                    npts = dataPtr[i+1] - dataPtr[i];
                    // Resize event
                    if (npts != npts0)
                    {
                        sacio_freeData(&grns->grns[kndx+i]);
                        grns->grns[kndx+i].data = sacio_malloc64f(npts);
                        grns->grns[kndx+i].npts = npts;
                        sacio_setIntegerHeader(SAC_INT_NPTS, npts,
                                               &grns->grns[kndx+i].header);
                        ierr = array_copy64f_work(npts,
                                                  &data[dataPtr[i]],
                                                  grns->grns[kndx+i].data);
                    }
                    else
                    {
                        ierr = array_copy64f_work(npts,
                                                  &data[dataPtr[i]],
                                                  grns->grns[kndx+i].data);
                    }
                }
                // Update the times
                if (lnewStartTime)
                {
                    for (i=0; i<nsignals; i++)
                    {
                        // Update the picks
                        for (l=0; l<nTimeVars; l++)
                        {
                            ierr = sacio_getFloatHeader(timeVars[l],
                                              grns->grns[kndx+i].header, &time);
                            if (ierr == 0)
                            {
                                time = time + epoch0; // Turn to real time
                                time = time - epoch;  // Relative to new time
                                sacio_setFloatHeader(timeVars[l], time,
                                                    &grns->grns[kndx+i].header);
                            }
                        } // Loop on picks
                        sacio_setEpochalStartTime(epoch,
                                                  &grns->grns[kndx+i].header);
                    } // Loop on signals
                }
                // Update the sampling period
                if (lnewDt)
                {
                    for (i=0; i<nsignals; i++)
                    {
                        sacio_setFloatHeader(SAC_FLOAT_DELTA, dt,
                                             &grns->grns[kndx+i].header);
                    }
                }
                process_freeSerialCommands(&commands);
                process_freeParallelCommands(&parallelCommands);
                memory_free64f(&data);
//TODO fix this
tdsearch_greens_repickGreensWithSTALTA(2.0, 10.0, 80.0, iobs, it, idep, grns);
            }
        }
    }
ERROR:;
    return 0;
}
//============================================================================//
int tdsearch_greens_writeSelectGreensFunctions(
    const char *dirnm,
    const int iobs, const int itstar, const int idepth,
    const struct tdSearchGreens_struct grns) 
{
    char fileName[PATH_MAX], rootName[PATH_MAX];
    size_t lenos;
    int i, ierr, indx;
    memset(rootName, 0, PATH_MAX*sizeof(char));
    if (dirnm == NULL)
    {
        strcpy(rootName, "./\0");
    }
    else
    {
        lenos = strlen(dirnm);
        if (lenos > 0)
        {
            strcpy(rootName, dirnm);
            if (rootName[lenos-1] != '/'){strcat(rootName, "/\0");}
        }
        else
        {
            strcpy(rootName, "./\0");
        }
    }
    if (!os_path_isdir(rootName))
    {
        ierr = os_makedirs(rootName);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to make output directory %s\n",
                    __func__, rootName);
            return -1;
        }
    }
    // Get the indices to write
    indx =  tdsearch_greens_getGreensFunctionIndex(G11_GRNS,
                                                   iobs, itstar, idepth,
                                                   grns);
    if (indx < 0)
    {
        fprintf(stderr, "%s: Invalid index\n", __func__);
        return -1;
    }
    for (i=0; i<6; i++)
    {
        memset(fileName, 0, PATH_MAX*sizeof(char));
        sprintf(fileName, "%s%s.%s.%s.%s.DEPTH_%d.TSTAR_%d.SAC",
                rootName, grns.grns[indx+i].header.knetwk,
                grns.grns[indx+i].header.kstnm,
                grns.grns[indx+i].header.kcmpnm,
                grns.grns[indx+i].header.khole, idepth, itstar);
        sacio_writeTimeSeriesFile(fileName, grns.grns[indx+i]);
    }
    return 0;
}
//============================================================================//
static int getPrimaryArrival(const struct sacHeader_struct hdr,
                             double *time, char phaseName[8])
{
    const enum sacHeader_enum timeVars[11]
       = {SAC_FLOAT_A,
          SAC_FLOAT_T0, SAC_FLOAT_T1, SAC_FLOAT_T2, SAC_FLOAT_T3,
          SAC_FLOAT_T4, SAC_FLOAT_T5, SAC_FLOAT_T6, SAC_FLOAT_T7,
          SAC_FLOAT_T8, SAC_FLOAT_T9};
    const enum sacHeader_enum timeVarNames[11]
       = {SAC_CHAR_KA,
          SAC_CHAR_KT0, SAC_CHAR_KT1, SAC_CHAR_KT2, SAC_CHAR_KT3,
          SAC_CHAR_KT4, SAC_CHAR_KT5, SAC_CHAR_KT6, SAC_CHAR_KT7,
          SAC_CHAR_KT8, SAC_CHAR_KT9};
    int i, ifound1, ifound2; 
    memset(phaseName, 0, 8*sizeof(char));
    for (i=0; i<11; i++)
    {
        ifound1 = sacio_getFloatHeader(timeVars[i], hdr, time);
        ifound2 = sacio_getCharacterHeader(timeVarNames[i], hdr, phaseName); 
        if (ifound1 == 0 && ifound2 == 0){return 0;}
    }
    printf("%s: Failed to get primary pick\n", __func__);
    *time =-12345.0;
    memset(phaseName, 0, 8*sizeof(char));
    strcpy(phaseName, "-12345"); 
    return -1;
}
