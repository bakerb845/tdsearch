#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parmt_utils.h"
#include "tdsearch_struct.h"
#include "tdsearch_hudson.h"
#include "iscl/log/log.h"

int tdsearch_greens_ffGreensToGreens(const struct tdSearchData_struct data,
                                     const struct tdSearchHudson_struct ffGrns )
{
    const char *fcnm = "tdsearch_greens_ffGreenToGreens\0";
    char kcmpnm[8];
    double *G, az, baz, cmpaz, cmpinc, cmpincSEED;
    int i, icomp, id, ierr, idx, iobs, it, kndx, npts;
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
double Gxx[512], Gyy[512], Gzz[512], Gxy[512], Gxz[512], Gyz[512];
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
                                                  Gxx, //grns->Gxx[idx].data,
                                                  Gyy, //grns->Gyy[idx].data,
                                                  Gzz, //grns->Gzz[idx].data,
                                                  Gxy, //grns->Gxy[idx].data,
                                                  Gxz, //grns->Gxz[idx].data,
                                                  Gyz); //grns->Gyz[idx].data);
                if (ierr != 0)
                {
                    log_errorF("%s: Failed to rotate Greens functions\n", fcnm);
                }
            }
        } 
    }
    return 0;
}
