#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iniparser.h>
#include "tdsearch_hudson.h"
#include "cps.h"
#include "iscl/array/array.h"
#include "iscl/fft/fft.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"
#include "iscl/sorting/sorting.h"

static int grd2ijk(const int igrd,
                   const int n1, const int n2, const int n3, 
                   int *i, int *j, int *k);
/*!
 * @brief Reads the ini file for Computer Programs in Seismology hudson96
 *        forward modeling variables.
 *
 * @param[in] iniFile    Name of ini file.
 *
 * @param[out] parms     The hudson96 parameters.
 *
 * @result 0 indicates success.
 *
 * ingroup tdsearch_hudson
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_hudson_readHudson96Parameters(const char *iniFile,
                                           struct hudson96_parms_struct *parms)
{
    const char *s;
    dictionary *ini;
    //------------------------------------------------------------------------//
    cps_setHudson96Defaults(parms);
    if (!os_path_isfile(iniFile))
    {
        fprintf(stderr, "%s: ini file: %s does not exist\n", __func__, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    if (ini == NULL)
    {
        fprintf(stderr, "%s: Cannot parse ini file\n", __func__);
        return -1;
    }
    // Teleseismic model
    s = iniparser_getstring(ini, "hudson96:modeltel", NULL);
    if (s != NULL)
    {
        strcpy(parms->modeltel, s);
    }
    else
    {
        strcpy(parms->modeltel, "tak135sph.mod\0");
    }
    // Receiver model
    s = iniparser_getstring(ini, "hudson96:modelrec", NULL);
    if (s != NULL)
    {
        strcpy(parms->modelrec, s);
    }
    else
    {
        strcpy(parms->modelrec, parms->modelrec);
    }
    // Source model
    s = iniparser_getstring(ini, "hudson96:modelsrc", NULL);
    if (s != NULL)
    {
        strcpy(parms->modelsrc, s);
    }
    else
    {
        strcpy(parms->modelsrc, parms->modelsrc);
    }
    parms->hs = iniparser_getdouble(ini, "hudson96:hs", 0.0);
    parms->dt = iniparser_getdouble(ini, "hudson96:dt", 1.0);
    parms->npts = iniparser_getint(ini, "hudson96:npts", 1024);
    parms->gcarc = iniparser_getdouble(ini, "hudson96:gcarc", 50.0);
    parms->offset = iniparser_getdouble(ini, "hudson96:offset", 10.0);
    parms->dosrc = iniparser_getboolean(ini, "hudson96:dosrc", 1);
    parms->dorec = iniparser_getboolean(ini, "hudson96:dorec", 1);
    parms->dotel = iniparser_getboolean(ini, "hudson96:dotel", 1);
    parms->dop = iniparser_getboolean(ini, "hudson96:dop", 1);
    parms->dokjar = iniparser_getboolean(ini, "hudson96:dokjar", 1);
    parms->loffsetdefault = iniparser_getboolean(ini,
                                                 "hudson96:loffsetdefault", 1);
    parms->verbose = iniparser_getboolean(ini, "hudson96:verbose", 0);
    parms->utstar = iniparser_getdouble(ini, "hudson96:utstar", -12345.0);
    parms->dottonly = iniparser_getboolean(ini, "hudson96:dottonly", 0);
    parms->zsrc = iniparser_getdouble(ini, "hudson96:zsrc", 100.0);
    parms->zrec = iniparser_getdouble(ini, "hudson96:zrec", 60.0);
    // Free ini dictionary
    iniparser_freedict(ini);
    return 0;
}
//============================================================================//
/*!
 * @brief Reads the ini file for the Computer Programs in Seismology hpulse96
 *        modeling variables.
 *
 * @param[in] iniFile     Name of ini file.
 *
 * @param[out] parms      The hpulse96 parameters.
 *
 * @result 0 indicates success.
 *
 * @ingroup tdsearch_hudson
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_hudson_readHpulse96Parameters(const char *iniFile,
                                           struct hpulse96_parms_struct *parms)
{
    const char *s;
    dictionary *ini;
    //------------------------------------------------------------------------//
    cps_setHpulse96Defaults(parms);
    if (!os_path_isfile(iniFile))
    {
        fprintf(stderr, "%s: ini file: %s does not exist\n", __func__, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    if (ini == NULL)
    {
        fprintf(stderr, "%s: Cannot parse ini file\n", __func__);
        return -1;
    }
    parms->ntau = iniparser_getint(ini, "hpulse96:ntau", -1);
    parms->alp = iniparser_getdouble(ini, "hpulse96:alp", -1.0);
    parms->ipt = iniparser_getint(ini, "hpulse96:ipt", -1);
    if (parms->ipt == 4)
    {
        s = iniparser_getstring(ini, "hpulse96:rfile", NULL);
        if (s != NULL)
        {
            strcpy(parms->rfile, s);
            if (!os_path_isfile(parms->rfile))
            {
                fprintf(stderr, "%s: Response file does not exist!\n", __func__);
                return -1;
            }
        }
        else
        {
            fprintf(stderr, "%s: Response file not specified!\n", __func__);
            return -1;
        }
    }
    parms->idva = iniparser_getint(ini, "hpulse96:idva", 1);
    parms->iodva = iniparser_getint(ini, "hpulse96:iodva", -1);
    parms->xmult = iniparser_getdouble(ini, "hpulse96:xmult", 1.0);
    parms->dozero = iniparser_getboolean(ini, "hpulse96:dozero", 0);

    if (parms->ipt >= 0 && parms->ipt <= 1 && parms->ntau <= 0)
    {
        parms->ntau = 1;
    }
    if (parms->ipt == 2 && parms->alp <= 0.0){parms->alp = 1.0;}
    if (parms->ipt == 2 && parms->alp < 0.0)
    {
        fprintf(stderr, "%s: No alpha for Ohnaka pulse\n", __func__);
        return -1;
    }
    if (parms->ipt < 0)
    {
        fprintf(stderr, "%s: No pulse shape defined\n", __func__);
        return -1;
    }
    if (parms->iodva < 0){parms->iodva = parms->idva;}
    iniparser_freedict(ini);
    return 0;
}
//============================================================================//
/*!
 * @brief Reads the forward modeling parameters for using hpulse and hudson.
 *
 * @param[in] iniFile    Name of ini file to read.
 *
 * @param[out] grns      On successful exit contains the forward modeling
 *                       parameters for hudson and hpulse.
 *
 * @result 0 indicates success.
 *
 * @ingroup tdsearch_hudson
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_hudson_initializeParametersFromIniFile(
    const char *iniFile,
    struct tdSearchHudson_struct *grns)
{
    const char *s;
    struct hudson96_parms_struct hudsonParms;
    struct hpulse96_parms_struct hpulseParms;
    dictionary *ini;
    int ierr, ierr1;
    ierr = 0;
    memset(grns, 0, sizeof(struct tdSearchHudson_struct));
    ierr1 = tdsearch_hudson_readHudson96Parameters(iniFile, &hudsonParms);
    if (ierr1 != 0)
    {
        fprintf(stderr, "%s: Error reading hudson96 parameters\n", __func__);
        ierr = ierr + 1;
    }
    ierr1 = tdsearch_hudson_readHpulse96Parameters(iniFile, &hpulseParms);
    if (ierr1 != 0)
    {   
        fprintf(stderr, "%s: Error reading hpulse96 parameters\n", __func__);
        ierr = ierr + 1;
    }
    // Set the modeling parameters on the modeling structure anyway
    tdsearch_hudson_setHudson96Parms(hudsonParms, grns);
    tdsearch_hudson_setHpulse96Parms(hpulseParms, grns);
    // Set a few more modeling parameters
    if (!os_path_isfile(iniFile))
    {
        fprintf(stderr, "%s: Error ini file doesn't exist\n", __func__);
        return -1;
    }
    ini = iniparser_load(iniFile);
    if (ini == NULL)
    {
        fprintf(stderr, "%s: Cannot parse ini file\n", __func__);
        return -1; 
    }
    grns->luseCrust1 = iniparser_getboolean(ini, "tdSearch:greens:useCrust\0",
                                            false);
    if (grns->luseCrust1)
    {
        s = iniparser_getstring(ini, "tdSearch:greens:crustDir\0", 
                                CPS_DEFAULT_CRUST1_DIRECTORY);
        if (!os_path_isdir(s))
        {
            fprintf(stderr, "%s: Error crust1.0 directory %s doesn't exist\n",
                     __func__, s);
            grns->luseCrust1 = false;
        }
        else
        {
            strcpy(grns->crust1Dir, s);
        }
    }
    grns->luseSrcModel
         = iniparser_getboolean(ini, "tdSearch:greens:useSourceModel\0", false);
    if (grns->luseSrcModel)
    {
        s = iniparser_getstring(ini, "tdSearch:greens:sourceModel\0", NULL);
        if (!os_path_isfile(s))
        {
            fprintf(stderr, "%s: Error source model %s doesn't exist\n", __func__, s);
            grns->luseSrcModel = false;
        }
        else
        {
            strcpy(grns->srcModelName, s);
        }
    }
    iniparser_freedict(ini);
    return ierr;
}
//============================================================================//
/*!
 * @brief Convenience utility for setting the hudson forward modeling
 *        parameters on the greens structure.
 *
 * @param[in] hudsonParms  The hudson96 forward modeling parameters.
 *
 * @param[out] grns        Contains the hudson forward modeling parameters.
 *
 * @result 0 indicates success.
 *
 * @ingroup tdsearch_hudson
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_hudson_setHudson96Parms(
    const struct hudson96_parms_struct hudsonParms,
    struct tdSearchHudson_struct *grns)
{
    cps_utils_copyHudson96ParmsStruct(hudsonParms,
                                      &grns->modelingParms.hudson96Parms);
    return 0;
}
//============================================================================//
/*!
 * @brief Convenience utility for setting the hpulse source time function
 *        parameters on the greens structure.
 *
 * @param[in] hpulseParms  The hpulse96 source time function parameters.
 *
 * @param[out] grns        Contains the hpulse STF parameters.
 *
 * @result 0 indicates success.
 *
 * @ingroup tdsearch_hudson
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_hudson_setHpulse96Parms(
    const struct hpulse96_parms_struct hpulseParms,
    struct tdSearchHudson_struct *grns) 
{
    cps_utils_copyHpulse96ParmsStruct(hpulseParms,
                                      &grns->modelingParms.hpulse96Parms);
    return 0;
}
//============================================================================//
/*
int tdsearch_hudson_setDistancesToModel(const double distMin,
                                        const double distMax,
                                        const struct tdSearchData_struct data,
                                        struct tdSearchHudson_struct *grns)
{
    double *gcarc;
    int *grnsToDataMap, *perm, ierr, k, nuse;
    const double tol = 0.001/111.195; // Stations are colocated at 1 m.
    if (data.nobs < 1)
    {
        fprintf(stderr, "%s: Error no data to model\n", __func__);
        return -1;
    }
    gcarc = memory_calloc64f(data.nobs+1);
    for (k=0; k<data.nobs; k++)
    {
        ierr = sacio_getFloatHeader(SAC_FLOAT_GCARC, data.obs[k].header,
                                    &gcarc[k]);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error - gcarc not set on waveform %d\n", __func__, k+1);
            gcarc[k] =-12345.0;
            continue;
        }
        if (gcarc[k] < distMin || gcarc[k] > distMax)
        {
            fprintf(stderr, "%s: Error gcarc %f not in modeling bounds [%f,%f]\n",
                       __func__, gcarc[k], distMin, distMax);
            gcarc[k] =-12345.0;
            continue;
        }
    }
    // Create a permutation
    perm = sorting_argsort64f(data.nobs, gcarc, SORT_ASCENDING, &ierr);
    grnsToDataMap = array_set32i(data.nobs, -1, &ierr);
    // Create a data to greens function map that is valid
    nuse = 0;
    for (k=0; k<data.nobs; k++)
    {
        if (gcarc[perm[k]] >= distMin && gcarc[perm[k]] <= distMax)
        {
            if (fabs(gcarc[perm[k]] - gcarc[perm[k+1]]) > tol)
            {
                grnsToDataMap[k] = nuse;
                nuse = nuse + 1;
            }
        }
    }
    memory_free32i(&perm);
    memory_free64f(&gcarc);
    return 0;
}
*/
//============================================================================//
//int tdsearch_hudson_setStructuralModels(struct tdSearchHudson_struct *grns)
//{
//    // Loop on the data
//}
//============================================================================//
/*!
 * @brief Releases memory on the tdSearchHudson structure.
 *
 * @param[out] grns   On exit all memory has been freed and variables
 *                    set to NULL or 0.
 *
 * @result 0 indicates success.
 *
 * @ingroup tdsearch_hudson
 *
 */
int tdsearch_hudson_free(struct tdSearchHudson_struct *grns)
{
    int i, nloop;
    nloop = 10*grns->nobs*grns->ntstar*grns->ndepth;
    memory_free64f(&grns->tstars);
    memory_free64f(&grns->depths);
    if (nloop > 0 && grns->grns != NULL)
    {
        for (i=0; i<nloop; i++)
        {
            sacio_free(&grns->grns[i]);
        }
        free(grns->grns);
    }
    return 0;
}
//============================================================================//
/*!
 * @brief Convenience routine to convert the observation, depth,
 *        and t* index to a global index in an array.
 *
 * @param[in] iobs   Observation index in the range [0,grns.nobs-1].
 * @param[in] idep   Depth index in the range [0,grns.ndepth-1]. 
 * @param[in] it     t* index in the range [0,grns.ntstar-1]. 
 * @param[in] grns   Green's functions structure which has the number
 *                   of observations, depths, and t*'s. 
 *
 * @result The index in the grns functions array corresponding to the
 *         given observation, t*, and depth. 
 *
 * @ingroup tdsearch_hudson
 *
 */
int tdsearch_hudson_observationDepthTstarToIndex(
    const int iobs, const int idep, const int it,
    const struct tdSearchHudson_struct grns) 
{
    int indx;
    indx = iobs*grns.ntstar*grns.ndepth + idep*grns.ntstar + it;
    return indx;
}
//============================================================================//
/*
int tdsearch_hudson_ffGreenToGreens(const struct tdSearchData_struct data,
                                    struct tdSearchHudson_struct *grns )
{
    double az, baz, cmpaz, cmpinc;
    int ierr, iobs;
    const double xmom = 1.0;     // no confusing `relative' magnitudes 
    const double xcps = 1.e-20;  // convert dyne-cm mt to output cm
    const double cm2m = 1.e-2;   // cm to meters
    const double dcm2nm = 1.e+7; // magnitudes intended to be specified in
                                 // Dyne-cm but I work in N-m
    // Given a M0 in Newton-meters get a seismogram in meters
    const double xscal = xmom*xcps*cm2m*dcm2nm;
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
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error reading header variables\n", __func__);
            break;
        }
        // Rescale the Green's functions

    }
    return 0;
}
*/
//============================================================================//
/*!
 * @brief Computes the fundamental fault Green's functions for the teleseismic
 *        body waves with hudson96.
 *
 * @param[in] data   Contains the data and SAC header information to which the
 *                   Green's functions will be contextualized.
 * @param[out] grns  Contains the Green's functions corresponding to the data.
 * 
 * @result 0 indicates success.
 *
 * @ingroup tdsearch_hudson
 *
 */
int tdsearch_hudson_computeGreensFF(const struct tdSearchData_struct data,
                                    struct tdSearchHudson_struct *grns)
{
    struct hudson96_parms_struct hudson96ParmsWork;
    struct hpulse96_parms_struct hpulse96ParmsWork;
    struct hpulse96_data_struct zresp;
    struct hwave_greens_struct *ffGrns;
    struct vmodel_struct *recmod, srcmod, telmod;
    char phaseName[8];
    double *lats, *lons, *dptr, cmpaz, cmpinc, dt, evla, evlo, gcarc, offset0,
           pickTime;
    int i, idep, ierr, ierrAll, ierr1, ierr2, indx, iobs, iobs0, ip, it,
        kndx, nloop, npts;
    bool lfound, lsh;
    const enum sacHeader_enum pickVars[11]
       = {SAC_FLOAT_A,
          SAC_FLOAT_T0, SAC_FLOAT_T1, SAC_FLOAT_T2, SAC_FLOAT_T3,
          SAC_FLOAT_T4, SAC_FLOAT_T5, SAC_FLOAT_T6, SAC_FLOAT_T7,
          SAC_FLOAT_T8, SAC_FLOAT_T9};
    const enum sacHeader_enum pickTypes[11]
       = {SAC_CHAR_KA,
          SAC_CHAR_KT0, SAC_CHAR_KT1, SAC_CHAR_KT2, SAC_CHAR_KT3,
          SAC_CHAR_KT4, SAC_CHAR_KT5, SAC_CHAR_KT6, SAC_CHAR_KT7,
          SAC_CHAR_KT8, SAC_CHAR_KT9};
    const char *cfaults[10] = {"ZDS\0", "ZSS\0", "ZDD", "ZEX\0",
                               "RDS\0", "RSS\0", "RDD", "REX\0",
                               "TDS\0", "TSS\0"};
    const int npTypes = 11;
    const int nrDist = 1;
    const int nrDepths = 1;
    const int nsDepths = 1;
    //------------------------------------------------------------------------//
    //
    // Set the modeling structures
    memset(&hudson96ParmsWork, 0, sizeof(struct hudson96_parms_struct));
    memset(&hpulse96ParmsWork, 0, sizeof(struct hpulse96_parms_struct));
    cps_setHudson96Defaults(&hudson96ParmsWork);
    cps_setHpulse96Defaults(&hpulse96ParmsWork);
    cps_utils_copyHudson96ParmsStruct(grns->modelingParms.hudson96Parms,
                                      &hudson96ParmsWork);
    cps_utils_copyHpulse96ParmsStruct(grns->modelingParms.hpulse96Parms,
                                      &hpulse96ParmsWork);
    offset0 = grns->modelingParms.hudson96Parms.offset;
    // Set the teleseismic model
    memset(&srcmod, 0, sizeof(struct vmodel_struct));
    memset(&telmod, 0, sizeof(struct vmodel_struct));
    recmod = (struct vmodel_struct *)
             calloc((size_t) data.nobs, sizeof(struct vmodel_struct));
    cps_globalModel_ak135f(&telmod);
    if (grns->luseCrust1)
    {
        printf("%s: Reading crust1.0...\n", __func__);
        lats = memory_calloc64f(data.nobs);
        lons = memory_calloc64f(data.nobs);
        ierr = sacio_getFloatHeader(SAC_FLOAT_EVLA, data.obs[0].header, &evla);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error - evla not set\n", __func__);
            return -1;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_EVLO, data.obs[0].header, &evlo);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error - evlo not set\n", __func__);
            return -1;
        }
        for (iobs=0; iobs<data.nobs; iobs++) 
        {
            sacio_getFloatHeader(SAC_FLOAT_STLA, data.obs[iobs].header,
                                 &lats[iobs]);
            sacio_getFloatHeader(SAC_FLOAT_STLO, data.obs[iobs].header,
                                 &lons[iobs]);
        }
        ierr = cps_crust1_getCrust1ForHerrmann(NULL, false,
                                               evla, evlo,
                                               lats, lons,
                                               data.nobs,
                                               &srcmod, recmod);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to load crust1.0 model\n", __func__);
            return -1;
        }
        memory_free64f(&lats);
        memory_free64f(&lons);
    }
    // Use teleseismic model for receiver model
    else
    {
        fprintf(stdout, "%s: Setting receiver models to teleseismic model\n", __func__);
        for (iobs=0; iobs<data.nobs; iobs++)
        {
            cps_utils_copyVmodelStruct(telmod, &recmod[iobs]);
        }
    }
    // Use source model
    if (grns->luseSrcModel)
    {
        cps_utils_freeVmodelStruct(&srcmod);
        ierr = cps_getmod(grns->srcModelName, &srcmod);
        if (ierr != 0){return EXIT_FAILURE;}
    }
    fprintf(stdout, "%s: Computing Green's functions...\n", __func__);
    // Otherwise source to teleseismic model
    if (!grns->luseSrcModel && !grns->luseCrust1)
    {
        fprintf(stdout, "%s: Setting source model to teleseismic model\n", __func__);
        cps_utils_copyVmodelStruct(telmod, &srcmod);
    }
    // Set space for output
    grns->nobs = data.nobs;
    nloop = data.nobs*grns->ndepth*grns->ntstar;
    grns->grns = (struct sacData_struct *)
                 calloc((size_t) (10*nloop), sizeof(struct sacData_struct));
    // Loop on the distances, depths, t*'s and compute greens functions
    iobs0 =-1;
    ierrAll = 0;
    for (indx=0; indx<nloop; indx++)
    {
        ffGrns = NULL;
        // Convert 3D grid into loop indices (i->n1 is fast; k->n3 is slow)
        ierr = grd2ijk(indx,
                       grns->ntstar, grns->ndepth, data.nobs,
                       &it, &idep, &iobs);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Failed to convert to grid\n", __func__);
            break;
        }
        if (data.lskip[iobs]){continue;}
        // New observation (station) -> may need to update velocity model
        if (iobs != iobs0)
        {
            iobs0 = iobs;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_DELTA,
                                    data.obs[iobs].header, &dt);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Could not get sampling period\n", __func__);
            break;
        }
        ierr = sacio_getFloatHeader(SAC_FLOAT_GCARC,
                                    data.obs[iobs].header, &gcarc);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Could not get gcarc\n", __func__);
            break;
        }
        ierr = sacio_getIntegerHeader(SAC_INT_NPTS,
                                      data.obs[iobs].header, &npts);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Could not get npts\n", __func__);
            break;
        }
        // Tie waveform modeling to first pick type 
        lfound = false;
        hudson96ParmsWork.dop = true;
        for (ip=0; ip<npTypes; ip++)
        {
            ierr1 = sacio_getCharacterHeader(pickTypes[ip],
                                             data.obs[iobs].header, phaseName);
            ierr2 = sacio_getFloatHeader(pickVars[ip], 
                                         data.obs[iobs].header, &pickTime);
            if (ierr1 == 0 && ierr2 == 0)
            {
                lfound = true;
                if (strncasecmp(phaseName, "P", 1) == 0)
                {
                    hudson96ParmsWork.dop = true;
                }
                else if (strncasecmp(phaseName, "S", 1) == 0)
                {
                    hudson96ParmsWork.dop = false;
                }
                else
                {
                    fprintf(stderr, "%s: Can't classify phase: %s\n", __func__, phaseName);
                    lfound = false;
                }
                break;
            }
        }
        if (!lfound)
        {
            fprintf(stdout, "%s: No pick available - won't be able to align\n", __func__);
            continue;
        }
        memset(&zresp, 0, sizeof(struct hpulse96_data_struct));
        hudson96ParmsWork.hs = grns->depths[idep];
        hudson96ParmsWork.dt = dt;
        hudson96ParmsWork.gcarc = gcarc;
        hudson96ParmsWork.utstar = grns->tstars[it];
        hudson96ParmsWork.npts = MIN(2048, 2*fft_nextpow2(MAX(1, npts), &ierr));
        //hudson96ParmsWork.offset = fmax(pickTime, offset0); //TODO: fmin or fmax?
        hudson96ParmsWork.offset = fmax(offset0,
                                  (double) (hudson96ParmsWork.npts - 1)*dt*0.2);
        if (pickTime < hudson96ParmsWork.offset)
        {
            hudson96ParmsWork.offset = (double) (int) (pickTime/dt + 0.5)*dt;
        }
        zresp = hudson96_interface(hudson96ParmsWork,
                                   telmod, recmod[iobs], srcmod, &ierr);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error calling hudson96\n", __func__);
            ierrAll = ierrAll + 1;
            goto NEXT_OBS;
        }
        ffGrns = hpulse96_interface(nrDist, nrDepths, nsDepths,
                                    hpulse96ParmsWork,
                                    &zresp, &ierr);
        if (ierr != 0)
        {
            fprintf(stderr, "%s: Error calling hpulse96 interface\n", __func__);
            ierrAll = ierrAll + 1;
            goto NEXT_OBS;
        }
        // Set the fundamental fault Green's functions
        for (i=0; i<10; i++)
        {
            dptr = NULL;
            lsh = false;
            if (i == 0)
            {
                dptr = ffGrns->zds;
            }
            else if (i == 1)
            {
                dptr = ffGrns->zss;
            }
            else if (i == 2)
            {
                dptr = ffGrns->zdd;
            }
            else if (i == 3)
            {
                dptr = ffGrns->zex;
            }
            else if (i == 4)
            {
                dptr = ffGrns->rds;
                cmpinc = 90.0;
            }
            else if (i == 5)
            {
                dptr = ffGrns->rss;
                cmpinc = 90.0;
            }
            else if (i == 6)
            {
                dptr = ffGrns->rdd;
                cmpinc = 90.0;
            }
            else if (i == 7)
            {
                dptr = ffGrns->rex;
                cmpinc = 90.0;
            }
            else if (i == 8)
            {
                dptr = ffGrns->tds;
                lsh = true;
                cmpaz = 90.0;
                cmpinc = 90.0;
            }
            else if (i == 9)
            {
                dptr = ffGrns->tss;
                lsh = true;
                cmpaz = 90.0;
                cmpinc = 90.0;
            }
            kndx = indx*10 + i;
            sacio_setDefaultHeader(&grns->grns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NPTS, ffGrns->npts,
                                   &grns->grns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZYEAR, 1970,
                                   &grns->grns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZJDAY, 1,
                                   &grns->grns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZHOUR, 0,
                                   &grns->grns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZMIN, 0,
                                   &grns->grns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZSEC, 0,
                                   &grns->grns[kndx].header);
            sacio_setIntegerHeader(SAC_INT_NZMSEC, 0,
                                   &grns->grns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_DELTA, ffGrns->dt,
                                 &grns->grns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_GCARC, ffGrns->dist/111.195,
                                 &grns->grns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_DIST, ffGrns->dist,
                                 &grns->grns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_EVDP, grns->depths[idep],
                                 &grns->grns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_STEL, ffGrns->stelel,
                                 &grns->grns[kndx].header); 
            sacio_setFloatHeader(SAC_FLOAT_AZ,  0.0,
                                 &grns->grns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_BAZ, 180.0,
                                 &grns->grns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_O, -ffGrns->t0, 
                                 &grns->grns[kndx].header); 
            sacio_setFloatHeader(SAC_FLOAT_B, ffGrns->t0,
                                 &grns->grns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_CMPAZ, cmpaz,
                                 &grns->grns[kndx].header);
            sacio_setFloatHeader(SAC_FLOAT_CMPINC, cmpinc,
                                 &grns->grns[kndx].header);
            if (hudson96ParmsWork.dop)
            {
                sacio_setFloatHeader(SAC_FLOAT_A, ffGrns->timep,
                                     &grns->grns[kndx].header);
                sacio_setCharacterHeader(SAC_CHAR_KA, "P",
                                         &grns->grns[kndx].header);
            }
            else
            {
                if (!lsh)
                {
                    sacio_setFloatHeader(SAC_FLOAT_A, ffGrns->timesv,
                                         &grns->grns[kndx].header);
                }
                else
                {
                    sacio_setFloatHeader(SAC_FLOAT_A, ffGrns->timesh,
                                         &grns->grns[kndx].header);
                }
                sacio_setCharacterHeader(SAC_CHAR_KA, "S",
                                         &grns->grns[kndx].header);
            }
            sacio_setCharacterHeader(SAC_CHAR_KO, "O\0",
                                     &grns->grns[kndx].header);
            sacio_setCharacterHeader(SAC_CHAR_KEVNM, "SYNTHETIC\0",
                                     &grns->grns[kndx].header);
            sacio_setCharacterHeader(SAC_CHAR_KCMPNM, cfaults[i],
                                     &grns->grns[kndx].header);
            grns->grns[kndx].npts = ffGrns->npts;
            grns->grns[kndx].data = sacio_malloc64f(ffGrns->npts);
            if (dptr != NULL)
            {
                array_copy64f_work(ffGrns->npts, dptr, grns->grns[kndx].data);
            }
            else
            {
                array_set64f_work(ffGrns->npts, 0.0, grns->grns[kndx].data);
            }
            dptr = NULL;
        }
NEXT_OBS:;
        if (ffGrns != NULL)
        {
            cps_utils_freeHwaveGreensStruct(ffGrns);
            free(ffGrns);
        }
        cps_utils_freeHpulse96DataStruct(&zresp);
    }
    // Release memory
    for (iobs=0; iobs<data.nobs; iobs++)
    {
        cps_utils_freeVmodelStruct(&recmod[iobs]);
    }
    free(recmod);
    cps_utils_freeVmodelStruct(&srcmod);
    cps_utils_freeVmodelStruct(&telmod);
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the t*, depth modeling grid on the greens structure.
 *
 * @param[in] ntstar   Number of t*'s (> 0).
 * @param[in] tstars   t* attenutation factors for use in hudson96. 
 * @param[in] ndepth   Number of depths (> 0).
 * @param[in] depths   Depths (km) at which to compute synthetics.
 *
 * @param[out] grns    Contains the t* and depths at which to model waveforms.
 *
 * @result 0 indicates success.
 *
 * @ingroup tdsearch_hudson
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_hudson_setGrid(const int ntstar, const double *__restrict__ tstars,
                            const int ndepth, const double *__restrict__ depths,
                            struct tdSearchHudson_struct *grns)
{
    int ierr;
    if (ntstar < 1 || ndepth < 1 || tstars == NULL || depths == NULL)
    {
        if (ntstar < 1)
        {
            fprintf(stderr, "%s: Error no t*'s to model\n", __func__);
        }
        if (ndepth < 1)
        {
            fprintf(stderr, "%s: Error no depths to model\n", __func__);
        }
        if (tstars == NULL)
        {
            fprintf(stderr, "%s: Error tstars cannot be NULL\n", __func__);
        }
        if (depths == NULL)
        {
            fprintf(stderr, "%s: Error depths cannot be NULL\n", __func__);
        }
        return -1;
    }
    if (array_min64f(ndepth, depths, &ierr) < 0.0)
    {
        fprintf(stderr, "%s: Error all depths must be >= 0\n", __func__);
        return -1;
    }
    grns->ntstar = ntstar;
    grns->ndepth = ndepth;
    grns->tstars = array_copy64f(ntstar, tstars, &ierr);
    grns->depths = array_copy64f(ndepth, depths, &ierr);
    return 0;
}
//============================================================================//
static int grd2ijk(const int igrd,
                   const int n1, const int n2, const int n3, 
                   int *i, int *j, int *k) 
{
    int ierr, n12;
    ierr = 0;
    n12 = n1*n2;
    *k = (igrd)/n12;
    *j = (igrd - *k*n12)/n1;
    *i =  igrd - *k*n12 - *j*n1;
    if (*i < 0 || *i > n1 - 1){ierr = ierr + 1;} 
    if (*j < 0 || *j > n2 - 1){ierr = ierr + 1;} 
    if (*k < 0 || *k > n3 - 1){ierr = ierr + 1;} 
    return ierr;
}
