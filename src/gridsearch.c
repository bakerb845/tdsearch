#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <math.h>
#include <float.h>
#include <iniparser.h>
#include "tdsearch_gridsearch.h"
#ifdef TDSEARCH_USE_INTEL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif 
#include "sacio.h"
#include "parmt_utils.h"
#include "iscl/array/array.h"
#include "iscl/log/log.h"
#include "iscl/memory/memory.h"
#include "iscl/os/os.h"
#include "iscl/signal/convolve.h"

#define NTSTAR 12   /*!< Default number of t*'s in the grid search. */
#define TSTAR0 0.1  /*!< Default first t* in grid search. */
#define TSTAR1 1.2  /*!< Default last t* in grid search. */
#define NDEPTH 12   /*!< Default number of depths in grid search. */
#define DEPTH0 0.0  /*!< First depth in grid search (km). */
#define DEPTH1 3.0  /*!< Last depth in grid search (km). */
#define DOUBLE_PAGE_SIZE sysconf(_SC_PAGE_SIZE)/sizeof(double)
#define PAGE_SIZE sysconf(_SC_PAGE_SIZE)
#define LDG 8

static int computePageSizePadding64f(const int n);

/*!
 * @brief Sets the default grid search parameters.
 *
 * @param[out] tds   On exit contains the default grid search settings such
 *                   as the t* range, depth range, and toggles off the max
 *                   time shift.
 *
 * @result 0 indicates success.
 *
 */
int tdsearch_gridSearch_setDefaults(struct tdSearch_struct *tds)
{
    int ierr;
    memset(tds, 0, sizeof(struct tdSearch_struct));
    tds->maxShiftTime = DBL_MAX;
    tds->lhaveMaxShiftTime = false;
    tds->lrecompute = true;
    ierr = tdsearch_gridSearch_defineTstarGrid(NTSTAR, TSTAR0, TSTAR1, tds);
    ierr = tdsearch_gridSearch_defineDepthGrid(NDEPTH, DEPTH0, DEPTH1, tds);
    return 0;
}
//============================================================================//
/*!
 * @brief Frees memory on the grid-search structure.
 *
 * @param[out] tds    On exit all memory has been freed and parameters set to
 *                    0 or NULL.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_gridSearch_free(struct tdSearch_struct *tds)
{
    memory_free64f(&tds->delayTime);
    memory_free64f(&tds->data);
    memory_free64f(&tds->Gmat);
    memory_free64f(&tds->tstar);
    memory_free64f(&tds->depths);
    memory_free64f(&tds->synthetic);
    memory_free64f(&tds->xcorr);
    memset(tds, 0, sizeof(struct tdSearch_struct));
    return 0;
}
//============================================================================//
/*!
 * @brief Defines the t* grid in grid search. 
 *
 * @param[in] ntstars    Number of t*'s in grid search (> 0).
 * @param[in] tstar0     First t* in grid search (>= 0).
 * @param[in] tstar1     Final t* in grid search (>= tstar0).
 *
 * @param[out] tds       Grid search structure with t* grid.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_gridSearch_defineTstarGrid(const int ntstars, const double tstar0,
                                        const double tstar1,
                                        struct tdSearch_struct *tds)
{
    const char *fcnm = "tdsearch_gridSearch_defineTstarGrid\0";
    int ierr;
    tds->ntstar = 0;
    tds->lrecompute = true;
    if (ntstars < 1 || tstar0 < 0.0 || tstar0 > tstar1)
    {
        if (ntstars < 1)
        {
            log_errorF("%s: Error no t*'s %d in search\n", fcnm, ntstars);
        }
        if (tstar0 < 0.0)
        {
            log_errorF("%s: Error tstar0 %f must be non-negative\n",
                       fcnm, tstar0);
        }
        if (tstar0 > tstar1)
        {
            log_errorF("%s: Error tstar0 %f > tstar1 %f\n",
                       fcnm, tstar0, tstar1);
        }
        return -1;
    }
    tds->ntstar = ntstars;
    memory_free64f(&tds->tstar);
    tds->tstar = array_linspace64f(tstar0, tstar1, ntstars, &ierr);
    return ierr;
}
//============================================================================//
/*!
 * @brief Defines the depth grid in grid search. 
 *
 * @param[in] ndepths    Number of depths in grid search (> 0).
 * @param[in] depth0     First depth (km) in grid search (>= 0).
 * @param[in] depth1     Final depth (km) in grid search (>= depth0).
 *
 * @param[out] tds       Grid search structure with depth grid.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_gridSearch_defineDepthGrid(const int ndepths, const double depth0,
                                        const double depth1,
                                        struct tdSearch_struct *tds)
{
    const char *fcnm = "tdsearch_gridSearch_defineDepthGrid\0";
    int ierr;
    tds->ndepth = 0;
    tds->lrecompute = true;
    if (ndepths < 1 || depth0 < 0.0 || depth0 > depth1)
    {   
        if (ndepths < 1)
        {
            log_errorF("%s: Error no depthss %d in search\n", fcnm, ndepths);
        }
        if (depth0 < 0.0)
        {
            log_errorF("%s: Error depth0 %f must be non-negative\n",
                       fcnm, depth0);
        }
        if (depth0 > depth1)
        {
            log_errorF("%s: Error depth0 %f > depth1 %f\n",
                       fcnm, depth0, depth1);
        }
        return -1; 
    }
    tds->ndepth = ndepths;
    memory_free64f(&tds->depths);
    tds->depths = array_linspace64f(depth0, depth1, ndepths, &ierr);
    return ierr;
}
//============================================================================//
/*!
 * @brief Reads the grid search parameters from the ini file.
 *
 * @param[in] iniFile    Name of ini file with parameters.
 *
 * @param[out] tds       On successful exit contains the grid search parameters.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_gridSearch_initializeFromFile(const char *iniFile,
                                           struct tdSearch_struct *tds)
{
    const char *fcnm = "tdsearch_gridSearch_initializeFromFile\0";
    double depth0, depth1, tstar0, tstar1;
    dictionary *ini;
    int ierr, ndepth, ntstar;
    // Set the defaults
    tdsearch_gridSearch_setDefaults(tds);
    // Verify the file exists
    if (!os_path_isfile(iniFile))
    {   
        log_errorF("%s: Error ini file %s doesn't exist\n", fcnm, iniFile);
        return -1; 
    }   
    ini = iniparser_load(iniFile);
    // Unpack the parameter file
    ierr = 1;
    tds->maxShiftTime
      = iniparser_getdouble(ini, "tdSearch:gridSearch:maxShiftTime\0", DBL_MAX);
    if (tds->maxShiftTime < DBL_MAX)
    {
        if (tds->maxShiftTime <= 0.0)
        {
            log_errorF("%s: Max shift time must be positive\n", fcnm);
            ierr = 1;
            goto ERROR;
        }
        tds->lhaveMaxShiftTime = true;
    }
    // Number of t*'s
    ntstar = iniparser_getint(ini, "tdSearch:gridSearch:ntstar\0", NTSTAR);
    if (ntstar < 1)
    {
        log_errorF("%s: number of t*'s must be positive\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    // Number of depths
    ndepth = iniparser_getint(ini, "tdSearch:gridSearch:ndepths\0", NDEPTH);
    if (ndepth < 1)
    {
        log_errorF("%s: number of depths must be positive\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    // Depth grid search range
    depth0 = iniparser_getdouble(ini, "tdSearch:gridSearch:depth0\0", DEPTH0);
    depth1 = iniparser_getdouble(ini, "tdSearch:gridSearch:depth1\0", DEPTH1);
    if (depth0 >= depth1)
    {
        log_errorF("%s: depth0 >= depth1 is invalid\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    // t* grid search range
    tstar0 = iniparser_getdouble(ini, "tdSearch:gridSearch:tstar0\0", TSTAR0);
    tstar1 = iniparser_getdouble(ini, "tdSearch:gridSearch:tstar1\0", TSTAR1);
    if (tstar0 >= tstar1)
    {
        log_errorF("%s: tstar0 >= tstar1 is invalid\n", fcnm);
        ierr = 1;
        goto ERROR;
    }
    // Create the grid
    ierr = tdsearch_gridSearch_defineTstarGrid(ntstar,
                                               tstar0, tstar1,
                                               tds);
    ierr = tdsearch_gridSearch_defineDepthGrid(ndepth,
                                               depth0, depth1,
                                               tds);
ERROR:;
    iniparser_freedict(ini);
    return 0;
}
//============================================================================//
/*!
 * @brief Returns the number of depths in the grid search. 
 *
 * @param[in] tds     Holds the number of t*'s.
 *
 * @result The number of t*'s in the grid search.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_gridSearch_getNumberOfTstars(const struct tdSearch_struct tds)
{
    return tds.ntstar;
}
//============================================================================//
/*!
 * @brief Returns the number of depths in the grid search.
 *
 * @param[in] tds     Holds the number of depths.
 *
 * @result The number of depths in the grid search.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_gridSearch_getNumberOfDepths(const struct tdSearch_struct tds)
{
    return tds.ndepth;
}
//============================================================================//
/*!
 * @brief Convenience utility for computing the index corresponding to the 
 *        2D row major depth x t* grid.
 *
 * @param[in] it      The t* index [0,ntstar-1].
 * @param[in] id      The depth indx [0,ndepth-1].
 * @param[in] ntstar  The number of t*'s in the grid search.
 *
 * @result The grid search index (= id*ntstar + it).
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_gridSearch_gridToIndex2(const int it, const int id,
                                     const int ntstar)
{
    int idt;
    idt = id*ntstar + it;
    return idt;
}
//============================================================================//
/*!
 * @brief Convenience utility for computing the index corresponding to the 
 *        2D row major depth x t* grid.
 *
 * @param[in] it      The t* index [0,ntstar-1].
 * @param[in] id      The depth indx [0,ndepth-1].
 * @param[in] tds     Grid search structure holding the number of ntstar's
 *                    in the grid search.
 *
 * @result The grid search index (= id*tds.ntstar + it).
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_gridSearch_gridToIndex(const int it, const int id,
                                    const struct tdSearch_struct tds)
{
    int idt, ntstar;
    ntstar = tdsearch_gridSearch_getNumberOfTstars(tds);
    idt = tdsearch_gridSearch_gridToIndex2(it, id, ntstar);
    return idt;
}
//============================================================================//
int tdsearch_gridSearch_setGmatFromSAC(const int ngrns,
                                       const struct sacData_struct data,
                                       const struct sacData_struct *ZDS,
                                       const struct sacData_struct *ZSS,
                                       const struct sacData_struct *ZDD,
                                       const struct sacData_struct *ZEX,
                                       const struct sacData_struct *RDS,
                                       const struct sacData_struct *RSS,
                                       const struct sacData_struct *RDD,
                                       const struct sacData_struct *REX,
                                       const struct sacData_struct *TDS,
                                       const struct sacData_struct *TSS,
                                       struct tdSearch_struct *tds)
{
    const char *fcnm = "tdsearch_gridSearch_setGmatFromSAC\0";
    double *Gxx, *Gyy, *Gzz, *Gxy, *Gxz, *Gyz, az, baz, cmpaz, cmpinc;
    char kcmpnm[8];
    int icomp, id, idt, ierr, indx, it, ndepth, ndt, npts, npts0, ntstar;
    ntstar = tdsearch_gridSearch_getNumberOfTstars(*tds);
    ndepth = tdsearch_gridSearch_getNumberOfDepths(*tds);
    ndt = ntstar*ndepth;
    if (ndt < 1)
    {
        log_errorF("%s: Error no points in grid search\n", fcnm);
        return -1;
    }
    if (ndt != ngrns)
    {
        log_errorF("%s: Error grid-search inconsistent with ngrns %d %d\n",
                   fcnm, ndt, ngrns);
        return -1;
    }
    // Space inquiry
    npts0 = 0;
    for (idt=0; idt<ngrns; idt++)
    {
        if (idt == 0)
        {
            ierr = sacio_getIntegerHeader(SAC_INT_NPTS, ZDS[idt].header, &npts0);
        }
        ierr = sacio_getIntegerHeader(SAC_INT_NPTS, ZDS[idt].header, &npts);
        if (ierr != 0)
        {
            log_errorF("%s: Error npts not set on grns %d\n", fcnm, idt); 
            return -1;
        }
        if (npts != npts0)
        {
            log_errorF("%s: Error inconsistent sizes on grns\n", fcnm);
            return -1; 
        }
    }
    // Extract the requisite header information from the data
    ierr = 0;
    ierr += sacio_getFloatHeader(SAC_FLOAT_BAZ, data.header, &baz);
    ierr += sacio_getFloatHeader(SAC_FLOAT_AZ,  data.header, &az);
    ierr += sacio_getFloatHeader(SAC_FLOAT_CMPAZ,  data.header, &cmpinc);
    ierr += sacio_getFloatHeader(SAC_FLOAT_CMPINC, data.header, &cmpaz);
    ierr += sacio_getCharacterHeader(SAC_CHAR_KCMPNM, data.header, kcmpnm);
    if (ierr != 0)
    {
        log_errorF("%s: Failed to get header information from data\n", fcnm);
        return -1;
    }
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
        log_errorF("%s: Can't classify component\n", fcnm);
        return -1;
    }
    cmpinc = cmpinc - 90.0; // SAC to SEED convention
    // Set space for Gmat
    memory_free64f(&tds->Gmat);
    tds->ldg = LDG;
    tds->ngrns = ngrns;
    tds->nptsGrns = npts;
    tds->mptsGrns = npts + computePageSizePadding64f(npts);
    tds->Gmat = memory_calloc64f(tds->ngrns*tds->mptsGrns*tds->ldg);
    // Set workspace
    Gxx = memory_calloc64f(npts);
    Gyy = memory_calloc64f(npts);
    Gzz = memory_calloc64f(npts);
    Gxy = memory_calloc64f(npts);
    Gxz = memory_calloc64f(npts);
    Gyz = memory_calloc64f(npts);
    for (id=0; id<ndepth; id++)
    {
        for (it=0; it<ntstar; it++)
        {
            idt = tdsearch_gridSearch_gridToIndex(it, id, *tds);
            sacio_getIntegerHeader(SAC_INT_NPTS, ZDS[idt].header, &npts);
            // Rotate into correct Green's functions
            ierr = parmt_utils_ff2mtGreens64f(npts, icomp, az, baz,
                                              cmpaz, cmpinc,
                                              ZDS[idt].data, ZSS[idt].data,
                                              ZDD[idt].data, ZEX[idt].data,
                                              RDS[idt].data, RSS[idt].data,
                                              RDD[idt].data, REX[idt].data,
                                              TDS[idt].data, TSS[idt].data,
                                              Gxx, Gyy, Gzz,
                                              Gxy, Gxz, Gyz);
            if (ierr != 0)
            {
                log_errorF("%s: Error mapping fund. faults to Greens fns\n",
                           fcnm);
                return -1;
            }
            indx = idt*tds->mptsGrns*tds->ldg;
            cblas_dcopy(npts, Gxx, 1, &tds->Gmat[indx+0], tds->ldg); 
            cblas_dcopy(npts, Gyy, 1, &tds->Gmat[indx+1], tds->ldg);
            cblas_dcopy(npts, Gzz, 1, &tds->Gmat[indx+2], tds->ldg);
            cblas_dcopy(npts, Gxy, 1, &tds->Gmat[indx+3], tds->ldg);
            cblas_dcopy(npts, Gxz, 1, &tds->Gmat[indx+4], tds->ldg);
            cblas_dcopy(npts, Gyz, 1, &tds->Gmat[indx+5], tds->ldg);
        }
    }
    memory_free64f(&Gxx);
    memory_free64f(&Gyy);
    memory_free64f(&Gzz);
    memory_free64f(&Gxy);
    memory_free64f(&Gxz);
    memory_free64f(&Gyz); 
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the synthetic seismogram on tdSearch structure corresponding
 *        to the previously set Green's functions and moment tensor for the 
 *        it'th t* and id'th depth. 
 *
 * @param[in] it          t* index [0,ntstar-1].
 * @param[in] id          Depth index [0,ndepth-1].
 *
 * @param[in,out] tds     On input the Green's functions matrix and moment
 *                        tensor has been set.
 *                        On output the synthetic has been computed.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */ 
int tdSearch_computeSynthetic(const int it, const int id,
                              struct tdSearch_struct *tds)
{
    const char *fcnm = "tdSearch_computeSynthetic\0";
    int idt, indx, jndx;
    idt = tdsearch_gridSearch_gridToIndex(it, id, *tds);
    if (!tds->lhaveMT || !tds->lhaveGrns ||
        it < 0 || it > tds->ntstar - 1 ||
        id < 0 || id > tds->ndepth - 1)
    {
        if (!tds->lhaveMT){log_errorF("%s: moment tensor not set\n", fcnm);}
        if (!tds->lhaveGrns){log_errorF("%s: grns fns not set\n", fcnm);}
        if (it < 0 || it > tds->ntstar - 1)
        {
            log_errorF("%s: Error invalid t* index %d\n", fcnm, it);
        } 
        if (id < 0 || id > tds->ndepth - 1)
        {
            log_errorF("%s: Error invalid depth index %d\n", fcnm, id);
        }
        return -1;
    }
    if (idt < 0 || idt > tds->ngrns - 1)
    {
        log_errorF("%s: Well this is strange; idt is wrong\n", fcnm);
        return -1;
    }
    indx = idt*tds->mptsGrns*tds->ldg;
    jndx = idt*tds->mptsGrns;
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                tds->nrowsGmat, tds->ncolsGmat, 1.0, &tds->Gmat[indx], tds->ldg,
                tds->mt, 1, 0.0, &tds->synthetic[jndx], 1);
    return 0;
}

int tdSearch_gridSearch_computeAllSynthetics(struct tdSearch_struct *tds)
{
    const char *fcnm = "tdSearch_gridSearch_computeAllSynthetics\0";
    int ndepth, ngrid, ntstar, nrowsAll;
    if (!tds->lhaveMT || !tds->lhaveGrns)
    {
        if (!tds->lhaveMT){log_errorF("%s: moment tensor not set\n", fcnm);}
        if (!tds->lhaveGrns){log_errorF("%s: grns fns not set\n", fcnm);}
        return -1; 
    }   
    ndepth = tdsearch_gridSearch_getNumberOfDepths(*tds);
    ntstar = tdsearch_gridSearch_getNumberOfTstars(*tds);
    ngrid = ndepth*ntstar;
    nrowsAll = ngrid*tds->mptsGrns;
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                nrowsAll, tds->ncolsGmat, 1.0, tds->Gmat, tds->ldg,
                tds->mt, 1, 0.0, tds->synthetic, 1); 
    return 0;
}

/*!
 * @brief Sets the moment tensor on the tdSearch structure
 *
 * @param[in] m6    moment tensor in NED system packed: 
 *                   \f$ \{ m_{xx}, m_{yy}, m_{zz}, 
 *                         m_{xy}, m_{xz}, m_{yz} \} \f$ 
 *
 * @param[out] tds  on successful exit contains the moment tensor
 *
 * @result 0 indicates success
 *
 * @author Ben Baker, ISTI
 *
 */
int tdSearch_gridSearch_setMomentTensor(const double *__restrict__ m6,
                                        struct tdSearch_struct *tds)
{
    const char *fcnm = "tdSearch_setMomentTensor\0";
    tds->lhaveMT = false;
    array_zeros64f_work(8, tds->mt);
    if (m6 == NULL)
    {
        log_errorF("%s: Error moment tensor cannot be NULL\n", fcnm);
        return -1;
    }
    array_copy64f_work(6, m6, tds->mt);
    tds->lhaveMT = true;
    return 0;
}
//============================================================================//
int tdSearch_gridSearch_getSynthetic(const int it, const int id,
                                     const struct tdSearch_struct tds,
                                     const int nwork, int *npts,
                                     double *__restrict__ synth)
{
    const char *fcnm = "tdSearch_gridSearch_getSynthetic\0";
    int idt, jndx;
    *npts = tds.nptsGrns;
    if (nwork < 0){return 0;}
    if (it < 0 || it > tds.ntstar - 1 || id < 0 || id > tds.ndepth - 1)
    {
        log_errorF("%s: Invalid indicies %d %d\n", fcnm, it, id);
        return -1;
    }
    if (*npts > nwork || synth == NULL)
    {
        if (*npts > nwork)
        {
            log_errorF("%s: nwork %d is too small - resize to %d\n",
                       fcnm, nwork, *npts);
            return -1;
        }
        if (synth == NULL)
        {
            log_errorF("%s: Error synth cannot be NULL\n", fcnm);
            return -1;
        }
    }
    idt = tdsearch_gridSearch_gridToIndex(it, id, tds);
    jndx = idt*tds.mptsGrns;
    array_copy64f_work(*npts, &tds.synthetic[jndx], synth);
    return 0;
}
//============================================================================//
/*!
 * @brief Computes the cross-correlation grid search over t* and depth
 *
 * @param[in,out] tds     On input contains Green's functions, data, and
 *                        grid search parameters.
 *                        On output contains the cross-correlation evaluated
 *                        at all points.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdSearch_gridSearch_computeGridSearch(struct tdSearch_struct *tds)
{
    const char *fcnm = "tdSearch_gridSearch_computeGridSearch\0";
    double *synth, *xcsig, ccmax, maxShiftTime, sqrtd2, sqrtg2;
    int id, idelay, idt, ierr, ierr1, imax, it, jndx, lc;
    //------------------------------------------------------------------------//
    //
    // make sure we have everything we need
    ierr = 0;
    if (!tds->lhaveData || !tds->lhaveGrns ||
        !tds->lhaveMT   || !tds->lhaveTD)
    {
        if (!tds->lhaveData){log_errorF("%s: data not set\n", fcnm);}
        if (!tds->lhaveGrns){log_errorF("%s: greens fns not set\n", fcnm);}
        if (!tds->lhaveMT){log_errorF("%s: moment tensor not set\n", fcnm);}
        if (!tds->lhaveTD){log_errorF("%s: gridsearch size not set\n", fcnm);}
        return -1;
    }
    maxShiftTime = DBL_MAX;
    if (tds->lhaveMaxShiftTime){maxShiftTime = tds->maxShiftTime;}
    // compute all synthetics
    tdSearch_gridSearch_computeAllSynthetics(tds);
    // get the energy in the data for the normalization 
    sqrtd2 = cblas_dnrm2(tds->npts, tds->data, 1);
    lc = tds->npts + tds->nptsGrns - 1;
    // Loop on the t* and depths
/*
    #pragma omp parallel \
     firstprivate (sqrtd2, lc) \
     private (ccmax, id, idelay, idt, ierr1, imax, jndx, it, nptsGrns, sqrtg2, xcsig) \
     shared (fcnm, maxShiftTime, tds) \
     default (none) reduction(+:ierr)
    {
*/
    xcsig = memory_calloc64f(lc);
    synth = memory_calloc64f(tds->nptsGrns);
/*
    #pragma for collapse(2)
*/
    for (id=0; id<tds->ndepth; id++)
    {
        for (it=0; it<tds->ntstar; it++)
        {
            idt = id*tds->ntstar + it;
            tds->xcorr[idt] =-1.0; // Set to a failure
 /*
            ierr = tdSearch_gridSearch_getSynthetic(it, id, *tds, tds->nptsGrns, &nptsGrns, synth);
 */
            // Compute the synthetic
            ierr1 = tdSearch_computeSynthetic(it, id, tds);
            if (ierr1 != 0)
            {
                log_errorF("%s: Error computing synthetic for (t*,d)=(%d,%d)\n",
                           fcnm, it, id);
                ierr = ierr + 1;
                continue;
            }
            // Compute the cross-correlation
            jndx = idt*tds->mptsGrns;
            ierr1 = signal_convolve_correlate64f_work(tds->npts,
                                                      tds->data,
                                                      tds->nptsGrns,
                                                      &tds->synthetic[jndx],
                                                      CONVCOR_FULL,
                                                      lc, xcsig);
            if (ierr1 != 0)
            {
                log_errorF("%s: Error computing xcorr for (t*,d)=(%d,%d)\n",
                           fcnm, it, id);
                ierr = ierr + 1;
                continue;
            }
            // Take the maximum positive cross-correlation
            imax = array_argmax64f(lc, xcsig);
            idelay = tds->nptsGrns - 1 - imax;
            sqrtg2 = cblas_dnrm2(tds->nptsGrns, &tds->synthetic[jndx], 1);
            ccmax = xcsig[imax]/(sqrtg2*sqrtd2);
            tds->delayTime[idt] = (double) idelay*tds->dt; 
            if (tds->delayTime[idt] > maxShiftTime) 
            {
                log_warnF("%s: Lag %f exceed max allowable time %f\n",
                          fcnm, tds->delayTime[idt], tds->maxShiftTime);
                ccmax =-1.0;
            }
            tds->xcorr[idt] = ccmax;
        } // Loop on t*
    } // Loop on depths
    memory_free64f(&synth);
    memory_free64f(&xcsig);
/*
    }
*/
    if (ierr != 0)
    {
        log_errorF("%s: Errors encountered during gridsearch\n", fcnm);
        return -1;
    }
    return 0;
}

static int computePageSizePadding64f(const int n)
{
    size_t mod, pad; 
    int ipad;
    // Set space and make G matrix
    pad = 0; 
    mod = ((size_t) n*sizeof(double))%PAGE_SIZE;
    if (mod != 0)
    {    
        pad = (PAGE_SIZE - mod)/sizeof(double);
    }    
    ipad = (int) pad; 
    return ipad;
}
