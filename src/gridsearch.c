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
#include "tdsearch_greens.h"
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
#define LDG 6

/*!
 * @brief Sets the default grid search parameters.
 *
 * @param[out] tds   On exit contains the default grid search settings such
 *                   as the t* range, depth range, and toggles off the max
 *                   time shift.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_gridSearch_setDefaults(struct tdSearch_struct *tds)
{
    const char *fcnm = "tdsearch_gridSearch_setDefaults\0";
    int ierr;
    memset(tds, 0, sizeof(struct tdSearch_struct));
    tds->maxShiftTime = DBL_MAX;
    ierr = tdsearch_gridSearch_defineTstarGrid(NTSTAR, TSTAR0, TSTAR1, tds);
    if (ierr != 0)
    {
        log_errorF("%s: Failed to set t* grid\n", fcnm);
        return -1;
    } 
    ierr = tdsearch_gridSearch_defineDepthGrid(NDEPTH, DEPTH0, DEPTH1, tds);
    if (ierr != 0)
    {
        log_errorF("%s: Failed to set depth grid\n", fcnm);
        return -1;
    }
    return ierr;
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
    memory_free64f(&tds->G);
    memory_free64f(&tds->Gxc);
    memory_free64f(&tds->tstar);
    memory_free64f(&tds->depths);
    memory_free64f(&tds->synthetic);
    memory_free64f(&tds->xcorr);
    memory_free32i(&tds->lags);
    memory_free32i(&tds->grnsMatrixPtr);
    memory_free32i(&tds->grnsXCMatrixPtr);
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
    ierr = 0;
    tdsearch_gridSearch_setDefaults(tds);
    // Verify the file exists
    if (!os_path_isfile(iniFile))
    {   
        log_errorF("%s: Error ini file %s doesn't exist\n", fcnm, iniFile);
        return -1; 
    }   
    ini = iniparser_load(iniFile);
    // Unpack the parameter file
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
    if (ierr != 0)
    {
        log_errorF("%s: Error defining t* grid\n", fcnm);
        goto ERROR;
    }
    ierr = tdsearch_gridSearch_defineDepthGrid(ndepth,
                                               depth0, depth1,
                                               tds);
    if (ierr != 0)
    {
        log_errorF("%s: Error defining depth grid\n", fcnm);
        goto ERROR;      
    }
ERROR:;
    iniparser_freedict(ini);
    return ierr;
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
int tdsearch_gridSearch_setForwardModelingMatrices(
    const int iobs, 
    const struct tdSearchData_struct data,
    const struct tdSearchGreens_struct grns,
    struct tdSearch_struct *tds)
{
    const char *fcnm = "tdsearch_gridSearch_setForwardModelingMatrices\0";
    double *G, *Gxc, *work;
    int *grnsMatrixPtr, *grnsXCMatrixPtr, indices[6],
        i, idep, idt, ierr, indx, it, jndx, kndx,
        lxc, lxcMax, lxc0, ntstar, ndepth, npgrns, npts;
    // Reset pointers
    memory_free64f(&tds->G);
    memory_free64f(&tds->Gxc); 
    memory_free32i(&tds->grnsMatrixPtr);
    memory_free32i(&tds->grnsXCMatrixPtr);
    // Check sizes
    ndepth = grns.ndepth;
    ntstar = grns.ntstar;
    if (ndepth != tds->ndepth || ntstar != tds->ntstar)
    {
        if (ndepth != tds->ndepth)
        {
            log_errorF("%s: Error ndepth %d != tds->ndepth %d\n",
                       fcnm, ndepth, tds->ndepth);
        }
        if (ntstar != tds->ntstar)
        {
            log_errorF("%s: Error ntstar %d != tds->ntstar %d\n",
                       fcnm, ntstar, tds->ntstar);
        }
        return -1;
    }
    npts = data.obs[iobs].npts;
    if (npts < 1)
    {
        log_errorF("%s: Error data has no points\n", fcnm);
        return -1;
    }
    tds->dnorm = cblas_dnrm2(npts, data.obs[iobs].data, 1);
    if (fabs(tds->dnorm) < DBL_EPSILON)
    {
        log_errorF("%s: Error - data is identically zero\n", fcnm);
        return -1;
    } 
    tds->npts = npts;
    sacio_getFloatHeader(SAC_FLOAT_DELTA, data.obs[iobs].header, &tds->dt);
    // Set the Green's functions forward modeling matrix size
    lxcMax = 0;
    grnsMatrixPtr = memory_calloc32i(ndepth*ntstar+1);
    grnsXCMatrixPtr = memory_calloc32i(ndepth*ntstar+1);
    grnsMatrixPtr[0] = 0;
    grnsXCMatrixPtr[0] = 0;
    for (idep=0; idep<ndepth; idep++)
    {
        for (it=0; it<ntstar; it++)
        {
            // Get the Green's functions indices
            ierr = tdsearch_greens_getGreensFunctionsIndices(iobs, it,
                                                             idep, grns,
                                                             indices);
            if (ierr != 0)
            {
                log_errorF("%s: Failed to get indices\n", fcnm);
                return -1;
            }
            // Compute the cross-correlation length at (dep, t*)
            idt = tdsearch_gridSearch_gridToIndex2(it, idep, ntstar);
            lxc0 = 0;
            for (i=0; i<6; i++)
            {
                indx = indices[i];
                npgrns = grns.grns[indx].npts;
                if (npgrns < 1)
                {
                    log_errorF("%s: Error - no points in grns fn\n", fcnm);
                    return -1;
                }
                lxc = npgrns + npts - 1;
                lxcMax = MAX(lxcMax, lxc);
                if (i == 0)
                {
                    lxc0 = lxc;
                    grnsMatrixPtr[idt+1] = grnsMatrixPtr[idt] + npgrns;
                    grnsXCMatrixPtr[idt+1] = grnsXCMatrixPtr[idt] + lxc;
                }
                if (lxc != lxc0)
                {
                    log_errorF("%s: Inconsistent column sizes\n", fcnm);
                    return -1;
                }
            }
        }
    }
    if (lxcMax < 1)
    {
        log_errorF("%s: Error - cross-correlations are empty\n", fcnm);
        return -1;
    }
    // Set the forward modeling matrices
    tds->nptsGrns = grnsMatrixPtr[ndepth*ntstar];
    G = memory_calloc64f(LDG*grnsMatrixPtr[ndepth*ntstar]);
    Gxc = memory_calloc64f(LDG*grnsXCMatrixPtr[ndepth*ntstar]);
    work = memory_calloc64f(lxcMax);
    for (idep=0; idep<ndepth; idep++)
    {
        for (it=0; it<ntstar; it++)
        {
            idt = tdsearch_gridSearch_gridToIndex2(it, idep, ntstar);
            // Get the Green's functions indices
            ierr = tdsearch_greens_getGreensFunctionsIndices(iobs, it, 
                                                             idep, grns,
                                                             indices);
            // Compute the cross-correlation length at (dep, t*)
            idt = tdsearch_gridSearch_gridToIndex2(it, idep, ntstar);
            jndx = LDG*grnsMatrixPtr[idt];
            kndx = LDG*grnsXCMatrixPtr[idt];
            for (i=0; i<6; i++)
            {
                indx = indices[i];
                npgrns = grns.grns[indx].npts;
                signal_convolve_correlate64f_work(npts, data.obs[iobs].data,
                                                  npgrns, grns.grns[indx].data,
                                                  CONVCOR_FULL,
                                                  lxc, work);
                //printf("%e %e\n",cblas_dnrm2(npgrns, grns.grns[indx].data, 1),
                //               cblas_dnrm2(lxc, work, 1));
                cblas_dcopy(npgrns, grns.grns[indx].data, 1, &G[jndx+i],   LDG);
                cblas_dcopy(lxc,    work,                 1, &Gxc[kndx+i], LDG);
            }
            //printf("%d\n", grnsMatrixPtr[idt+1]);
        }
    }
    tds->G = G;
    tds->Gxc = Gxc;
    tds->grnsMatrixPtr = grnsMatrixPtr;
    tds->grnsXCMatrixPtr = grnsXCMatrixPtr;
    memory_free64f(&work);
    return 0;
}
//============================================================================//
/*!
 * @brief Sets the moment tensor on the tdSearch structure.
 *
 * @param[in] m6    moment tensor in NED system packed: 
 *                   \f$ \{ m_{xx}, m_{yy}, m_{zz}, 
 *                         m_{xy}, m_{xz}, m_{yz} \} \f$.
 *
 * @param[out] tds  On successful exit contains the moment tensor.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdSearch_gridSearch_setMomentTensor(const double *__restrict__ m6,
                                        struct tdSearch_struct *tds)
{
    const char *fcnm = "tdSearch_setMomentTensor\0";
    tds->lhaveMT = false;
    array_zeros64f_work(8, tds->mt); // NOTE - this is padded
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
int tdSearch_gridSearch_setMomentTensorFromElements(
    const double m11, const double m22, const double m33,
    const double m12, const double m13, const double m23,
    const enum compearthCoordSystem_enum basis,
    struct tdSearch_struct *tds)
{
    const char *fcnm = "tdSearch_gridSearch_setMomentTensorFromElements\0";
    double mt[6], mtNED[6];
    int ierr;
    tds->lhaveMT = false;
    mt[0] = m11;
    mt[1] = m22;
    mt[2] = m33;
    mt[3] = m12;
    mt[4] = m13;
    mt[5] = m23;
    ierr = compearth_convertMT(basis, NED, mt, mtNED);
    if (ierr != 0)
    {
        printf("%s: Failed to set moment tensor\n", fcnm);
        return -1;
    }
    ierr = tdSearch_gridSearch_setMomentTensor(mtNED, tds);
    if (ierr != 0)
    {
        printf("%s: Failed to set moment tensor\n", fcnm);
        return -1;
    }
    return 0;
}
//============================================================================//
int tdSearch_gridSearch_setMomentTensorFromEvent(
    const struct tdSearchEventParms_struct event,
    struct tdSearch_struct *tds)
{
    const char *fcnm = "tdSearch_gridSearch_setMomentTensorFromEvent\0";
    int ierr;
    ierr = tdSearch_gridSearch_setMomentTensorFromElements(
                    event.m11, event.m22, event.m33,
                    event.m12, event.m13, event.m23,
                    event.basis, tds);
    if (ierr != 0)
    {
        log_errorF("%s: Error setting moment tensor\n", fcnm);
        return -1;
    }
    return 0;
}
//============================================================================//
int tdSearch_gridSearch_performGridSearch(struct tdSearch_struct *tds)
{
    const char *fcnm = "tdSearch_gridSearch_performGridSearch\0";
    double *u, *ud, dnorm, unorm;
    int id, idt, indx, iopt, it, jndx, l1, l2, lxc, maxlags,
        ndepth, ndt, nlag, npgrns, nrows, ntstar;
    const int ncols = 6;
    enum isclError_enum ierr;
    tds->itopt =-1;
    tds->idopt =-1;
    ndepth = tds->ndepth;
    ntstar = tds->ntstar;
    ndt = ntstar*ndepth;
    dnorm = tds->dnorm;
    memory_free64f(&tds->xcorr);
    memory_free64f(&tds->synthetic);
    memory_free32i(&tds->lags);
    // Some easy checks
    if (!tds->lhaveMT)
    {
        log_errorF("%s: Error moment tensor not set\n", fcnm);
        return -1;
    }
    if (tds->G == NULL || tds->Gxc == NULL ||
        tds->grnsMatrixPtr == NULL || tds->grnsXCMatrixPtr == NULL)
    {
        log_errorF("%s: You need to set the forward modeling matrices\n", fcnm);
        return -1;
    }
    maxlags =-1;
    if (tds->maxShiftTime < DBL_MAX)
    {
        maxlags = (int) (tds->maxShiftTime/tds->dt + 0.5);
    }
    // Set output space
    tds->xcorr = memory_calloc64f(ndt);
    tds->lags = memory_calloc32i(ndt);
    // Compute the cross-correlation
    nrows = tds->grnsXCMatrixPtr[ndt];
    ud = memory_calloc64f(nrows);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                nrows, ncols, 1.0, tds->Gxc, LDG,
                tds->mt, 1, 0.0, ud, 1);
    // Compute the synthetics
    nrows = tds->grnsMatrixPtr[ndt]; 
    u = memory_calloc64f(nrows);
    cblas_dgemv(CblasRowMajor, CblasNoTrans,
                nrows, ncols, 1.0, tds->G, LDG,
                tds->mt, 1, 0.0, u, 1);
    // Now do the grid search
    for (id=0; id<ndepth; id++)
    {
        for (it=0; it<ntstar; it++)
        {
            idt = tdsearch_gridSearch_gridToIndex(it, id, *tds);
            indx = tds->grnsMatrixPtr[idt];
            jndx = tds->grnsXCMatrixPtr[idt];
            lxc = tds->grnsXCMatrixPtr[idt+1]
                - tds->grnsXCMatrixPtr[idt];
            npgrns = tds->grnsMatrixPtr[idt+1]
                   - tds->grnsMatrixPtr[idt];
            unorm = cblas_dnrm2(npgrns, &u[indx], 1);
            if (maxlags < 0)
            {
                iopt = array_argmax64f(lxc, &ud[jndx], &ierr);
            }
            else
            {
                l1 = MAX(0,  npgrns - maxlags);
                l2 = MIN(lxc-1, npgrns + maxlags);
                nlag = l2 - l1 + 1;
                iopt = l1 + array_argmax64f(nlag, &ud[jndx+l1], &ierr);
 /*
for (int k=0; k<lxc; k++)
{
printf("%d %e\n",k, ud[jndx+k]);
}
getchar();
*/
            }
            tds->xcorr[idt] = ud[jndx+iopt]/(dnorm*unorm);
            // Compute the lag where the times range from [-npts+1:npgrns-1]
            tds->lags[idt] =-npgrns + 1 + iopt;
            //printf("%f %f %f %d\n", tds->depths[id], tds->tstar[it], xc, idelay);
        }
        //printf("\n");
    }
    tds->synthetic = u;
    memory_free64f(&ud);
    // Compute the optimum indices
    iopt = array_argmax64f(ndt, tds->xcorr, &ierr);
    tds->idopt = iopt/tds->ntstar;
    tds->itopt = iopt - tds->idopt*tds->ntstar;
    if (tdsearch_gridSearch_gridToIndex(tds->itopt, tds->idopt, *tds) != iopt)
    {
        log_warnF("%s: Switching to slow optimum computation\n", fcnm);
        for (id=0; id<ndepth; id++)
        {
            for (it=0; it<ntstar; it++)
            {
                idt = tdsearch_gridSearch_gridToIndex(it, id, *tds);
                if (iopt == idt)
                {
                    tds->idopt = id;
                    tds->itopt = it;
                } 
            }
        }
    }
    return 0;
}
int tdsearch_gridSearch_writeHeatMap(const char *dirnm, const char *froot,
                                     const struct tdSearch_struct tds)
{
    const char *fcnm = "tdsearch_gridSearch_writeHeatMap\0";
    FILE *ofl;
    char fname[PATH_MAX], fgnu[PATH_MAX], cmd[256];
    int id, idt, it;
    memset(fname, 0, PATH_MAX*sizeof(char));
    if (!tds.xcorr || !tds.lags)
    {
        log_errorF("%s: Error gridsearch likely not performed\n", fcnm);
        return -1;
    }
    if (!tds.tstar || !tds.depths)
    {
        log_errorF("%s: Error tds likely not initialized\n", fcnm);
        return -1;
    }
    if (dirnm != NULL)
    {
        sprintf(fname, "%s/%s.surf_txt", dirnm, froot);
    }
    else
    {
        sprintf(fname, "./%s.surf_txt", froot);
    }
    ofl = fopen(fname, "w");
    for (id=0; id<tds.ndepth; id++)
    {
        for (it=0; it<tds.ntstar; it++)
        {
            idt = tdsearch_gridSearch_gridToIndex(it, id, tds);
            fprintf(ofl, "%f %f %e %d\n",
                    tds.tstar[it], tds.depths[id],
                    tds.xcorr[idt], tds.lags[idt]);
        }
       fprintf(ofl, "\n");
    } 
    fclose(ofl);
    // and the gnuplot file
    memset(fgnu, 0, PATH_MAX*sizeof(char));
    if (dirnm != NULL)
    {
        sprintf(fgnu, "%s/%s.gnu", dirnm, froot);
    }
    else
    {
        sprintf(fgnu, "./%s.gnu", froot);
    }
    ofl = fopen(fgnu, "w");
    fprintf(ofl, "#!/usr/bin/gnuplot -persist\n");
    fprintf(ofl, "set pm3d map\n");
    fprintf(ofl, "set xlabel 't*'\n");
    fprintf(ofl, "set ylabel 'Depth (km)'\n");
    fprintf(ofl, "set title 'Cross-Correlation'\n");
    fprintf(ofl, "set grid\n");
    fprintf(ofl, "set xrange [%f:%f]\n", tds.tstar[0],
                                         tds.tstar[MAX(0, tds.ntstar-1)]);
    fprintf(ofl, "set yrange [%f:%f]\n", tds.depths[MAX(0, tds.ndepth-1)],
                                         tds.depths[0]);
    fprintf(ofl, "splot '%s.surf_txt' u 1:2:3\n", froot);
    fclose(ofl);
    // change permissions
    memset(cmd, 0, 256*sizeof(char));
    sprintf(cmd, "chmod 0755 %s", fgnu);
    system(cmd);
    return 0;
}
//============================================================================//
int tdSearch_gridSearch_makeSACSynthetic(
    const int iobs, const int it, const int id,
    const struct tdSearchData_struct data,
    const struct tdSearchGreens_struct grns,
    const struct tdSearch_struct tds,
    struct sacData_struct *synth)
{
    const char *fcnm = "tdSearch_gridSearch_getSynthetic\0";
    double epoch;
    int idt, ierr, igrns, npGrns;
    memset(synth, 0, sizeof(struct sacData_struct));
    // Pick off the header
    if (it < 0 || it >= tds.ntstar || id < 0 || id >= tds.ndepth)
    {
        log_errorF("%s: Invalid (t*,depth) index (%d,%d)\n", fcnm, it, id);
        return - 1;
    }
    if (tds.synthetic == NULL || tds.lags == NULL)
    {
        if (tds.synthetic == NULL)
        {
            log_errorF("%s: Error tds.synethetic is NULL\n", fcnm);
        }
        if (tds.lags == NULL)
        {
            log_errorF("%s: Error tds.lags is NULL\n", fcnm);
        }
        return -1;
    }
    idt = tdsearch_gridSearch_gridToIndex(it, id, tds);
    igrns = tdsearch_greens_getGreensFunctionIndex(G11_GRNS, iobs,
                                                   it, id, grns);
    if (idt < 0 || igrns < 0)
    {
        log_errorF("%s: Failed getting an index\n", fcnm);
        return -1;
    }
    // Copy header
    sacio_copyHeader(grns.grns[igrns].header, &synth->header); 
    // Copy corresponding synthetic
    synth->npts = grns.grns[igrns].npts;
    synth->data = sacio_malloc64f(synth->npts);
    npGrns = tds.grnsMatrixPtr[idt+1] - tds.grnsMatrixPtr[idt];
    if (npGrns != synth->npts)
    {
        log_errorF("%s: Size mismatch\n", fcnm);
        sacio_free(synth);
        return -1; 
    }
    array_copy64f_work(npGrns, &tds.synthetic[tds.grnsMatrixPtr[idt]],
                       synth->data);
    // Get the component inline with the estimate
    sacio_setCharacterHeader(SAC_CHAR_KCMPNM, data.obs[iobs].header.kcmpnm,
                             &synth->header);
    // Fix the timing
    ierr = sacio_getEpochalStartTime(synth->header, &epoch);
    if (ierr != 0)
    {
        log_errorF("%s: Failed to get start time\n", fcnm);
        return -1;
    }
    epoch = epoch + (double) tds.lags[idt]*tds.dt;
    sacio_setEpochalStartTime(epoch, &synth->header);
    return 0;
}
//============================================================================//
