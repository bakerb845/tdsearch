#ifndef TDSEARCH_HUDSON_H__
#define TDSEARCH_HUDSON_H__ 1
#include "tdsearch_struct.h"

#ifdef __cplusplus
extern "C"
{
#endif

int tdsearch_hudson_computeGreensFF(const struct tdSearchData_struct data,
                                    struct tdSearchHudson_struct *grns);
int tdsearch_hudson_setDistancesToModel(const double distMin,
                                        const double distMax,
                                        const struct tdSearchData_struct data,
                                        struct tdSearchHudson_struct *grns);
int tdsearch_hudson_setGrid(const int ntstar, const double *__restrict__ tstars,
                            const int ndepth, const double *__restrict__ depths,
                            struct tdSearchHudson_struct *grns);
int tdsearch_hudson_initializeParametersFromIniFile(
    const char *iniFile,
    struct tdSearchHudson_struct *grns);
int tdsearch_hudson_setHudson96Parms(
    const struct hudson96_parms_struct hudsonParms,
    struct tdSearchHudson_struct *grns);
int tdsearch_hudson_setHpulse96Parms(
    const struct hpulse96_parms_struct hpulseParms,
    struct tdSearchHudson_struct *grns);
int tdsearch_hudson_free(struct tdSearchHudson_struct *grns);
int tdsearch_hudson_observationDepthTstarToIndex(
    const int iobs, const int idep, const int it,
    const struct tdSearchHudson_struct grns);

#ifdef __cplusplus
}
#endif
#endif
