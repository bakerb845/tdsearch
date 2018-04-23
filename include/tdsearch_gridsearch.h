#ifndef TDSEARCH_GRIDSEARCH_H__
#define TDSEARCH_GRIDSEARCH_H__ 1
#include "tdsearch_struct.h"

/*!
 * @brief Functions for configuring and applying the
 *        cross-correlation-based grid search.
 * @defgroup tdsearch_gridsearch Grid Search
 */

#ifdef __cplusplus
extern "C"
{
#endif

int tdsearch_gridSearch_setForwardModelingMatrices(
    const int iobs, 
    const struct tdSearchData_struct data,
    const struct tdSearchGreens_struct grns,
    struct tdSearch_struct *tds);
int tdSearch_gridSearch_setMomentTensorFromElements(
    const double m11, const double m22, const double m33,
    const double m12, const double m13, const double m23,
    const enum compearthCoordSystem_enum basis,
    struct tdSearch_struct *tds);
int tdSearch_gridSearch_setMomentTensorFromEvent(
    const struct tdSearchEventParms_struct event,
    struct tdSearch_struct *tds);
int tdSearch_gridSearch_makeSACSynthetic(
    const int iobs, const int it, const int id,
    const struct tdSearchData_struct data,
    const struct tdSearchGreens_struct grns,
    const struct tdSearch_struct tds,
    struct sacData_struct *synth);
int tdsearch_gridSearch_writeHeatMap(const char *dirnm, const char *fname,
                                     const struct tdSearch_struct tds);
int tdSearch_gridSearch_performGridSearch(struct tdSearch_struct *tds);
int tdsearch_gridSearch_defineTstarGrid(const int ntstars, const double tstar0,
                                        const double tstar1,
                                        struct tdSearch_struct *tds);
int tdsearch_gridSearch_defineDepthGrid(const int ndepths, const double depth0,
                                        const double depth1,
                                        struct tdSearch_struct *tds);
int tdsearch_gridSearch_initializeFromFile(const char *iniFile,
                                           struct tdSearch_struct *tds);
int tdsearch_gridSearch_free(struct tdSearch_struct *tds);

#ifdef __cplusplus
}
#endif
#endif
