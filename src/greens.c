#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parmt_utils.h"
#ifdef TDSEARCH_USE_INTEL
#include <mkl_blas.h>
#else
#include <cblas.h>
#endif
#include "tdsearch_struct.h"
#include "tdsearch_hudson.h"
#include "iscl/log/log.h"

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
    const enum tdSearchGreens_enum GMT_TERM,
    const int iobs, const int itstar, const int idepth,
    const struct tdSearchGreens_struct grns)
{
    const char *fcnm = "tdsearch_greens_getGreensFunctionIndex\0";
    int igx, indx;
    indx =-1;
    igx = (int) GMT_TERM - 1;
    if (igx < 0 || igx > 5)
    {
        log_errorF("%s: Can't classify Green's functions index\n", fcnm);
        return indx;
    }
    indx = iobs*(6*grns.ntstar*grns.ndepth)
         + idepth*(6*grns.ntstar)
         + itstar*6
         + igx;
    return indx;
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
    int k;
    if (grns->grns != NULL && grns->ngrns > 0)
    {
        for (k=0; k<grns->ngrns; k++)
        {
            sacio_free(&grns->grns[k]);
        }
        free(grns->grns);
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
    const char *fcnm = "tdsearch_greens_ffGreenToGreens\0";
    char kcmpnm[8];
    double *G, az, baz, cmpaz, cmpinc, cmpincSEED;
    int i, icomp, id, ierr, idx, indx, iobs, it, kndx, npts;
    size_t nalloc;
    const char *kcmpnms[6] = {"GXX\0", "GYY\0", "GZZ\0",
                              "GXY\0", "GXZ\0", "GYZ\0"};
    const double xmom = 1.0;     // no confusing `relative' magnitudes 
    const double xcps = 1.e-20;  // convert dyne-cm mt to output cm
    const double cm2m = 1.e-2;   // cm to meters
    const double dcm2nm = 1.e+7; // magnitudes intended to be specified in
                                 // Dyne-cm but I work in N-m
    // Given a M0 in Newton-meters get a seismogram in meters
    const double xscal = xmom*xcps*cm2m*dcm2nm;
    memset(grns, 0, sizeof(struct tdSearchGreens_struct));
    grns->ntstar = ffGrns.ntstar;
    grns->ndepth = ffGrns.ndepth;
    grns->nobs = data.nobs;
    grns->ngrns = 6*grns->ntstar*grns->ndepth*grns->nobs;
    if (grns->ngrns < 1)
    {
        log_errorF("%s: Error grns is empty\n", fcnm);
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
        ierr += sacio_getCharacterHeader(SAC_CHAR_KCMPNM,
                                         data.obs[iobs].header, kcmpnm);
        if (ierr != 0)
        {
            log_errorF("%s: Error reading header variables\n", fcnm);
            break;
        }
        cmpincSEED = cmpinc - 90.0; // SAC to SEED convention
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
            log_errorF("%s: Can't classify component: %s\n", fcnm, kcmpnm);
        }
        // Process all Green's functions in this block
        for (id=0; id<ffGrns.ndepth; id++)
        {
            for (it=0; it<ffGrns.ntstar; it++)
            {
                idx = tdsearch_hudson_observationDepthTstarToIndex(iobs, id, it,
                                                                   ffGrns);
                kndx = 10*idx;
                npts = ffGrns.grns[kndx].npts;
                indx = tdsearch_greens_getGreensFunctionIndex(G11_GRNS,
                                                              iobs, it, id,
                                                              *grns);
                for (i=0; i<6; i++)
                {
                    sacio_copy(ffGrns.grns[kndx], &grns->grns[indx+i]);
                    sacio_setFloatHeader(SAC_FLOAT_AZ, az,
                                         &grns->grns[indx].header);
                    sacio_setFloatHeader(SAC_FLOAT_BAZ, baz,
                                         &grns->grns[indx].header);
                    sacio_setFloatHeader(SAC_FLOAT_CMPAZ, cmpaz,
                                         &grns->grns[indx].header);
                    sacio_setFloatHeader(SAC_FLOAT_CMPINC, cmpinc,
                                         &grns->grns[indx].header);
                    sacio_setCharacterHeader(SAC_CHAR_KCMPNM, kcmpnms[i],
                                             &grns->grns[indx].header);
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
                    log_errorF("%s: Failed to rotate Greens functions\n", fcnm);
                }
                // Fix the scaling issue from Programs in Seismology
                for (i=0; i<6; i++)
                {
                    cblas_dscal(npts, xscal, grns->grns[indx+i].data, 1); 
                }
            }
        } 
    }
    return 0;
}
