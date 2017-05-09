#ifndef TDSEARCH_DATA_H__
#define TDSEARCH_DATA_H__ 1
#include "tdsearch_struct.h"

#ifdef __cplusplus
extern "C"
{
#endif
int tdsearch_data_free(struct tdSearchData_struct *data);
int tdsearch_data_initializeFromFile(const char *iniFile,
                                     struct tdSearchData_struct *data);
int tdsearch_data_readFiles(const int nfiles,
                            const char **sacfls,
                            const char **pzfls,
                            struct tdSearchData_struct *data);
int tdsearch_data_setEventInformation(const double evla,
                                      const double evlo,
                                      const double evdp,
                                      const double evtime,
                                      struct tdSearchData_struct *data);
int tdsearch_data_computeTheoreticalPPickTimes(
    const char *dirnm, const char *model,
    const struct tdSearchData_struct data,
    const int nwork, double *__restrict__ ptimes);
int tdsearch_data_setPPickTimeFromTheoreticalTime(
    const char *dirnm, const char *model,
    const enum sacHeader_enum pickHeaderTime,
    const enum sacHeader_enum pickHeaderName,
    struct tdSearchData_struct *data);
int tdsearch_data_setPickTime(const double pickTime,
                              const char *phaseName,
                              const enum sacHeader_enum pickHeaderTime,
                              const enum sacHeader_enum pickHeaderName,
                              struct sacData_struct *data);

#ifdef __cplusplus
}
#endif

#endif
