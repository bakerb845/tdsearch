#ifndef TDSEARCH_GRIDSEARCH_H__
#define TDSEARCH_GRIDSEARCH_H__ 1
#include "tdsearch_struct.h"

#ifdef __cplusplus
extern "C"
{
#endif

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
