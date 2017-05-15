#ifndef TDSEARCH_GREENS_H__
#define TDSEARCH_GREENS_H__ 1
#include "tdsearch_struct.h"

#ifdef __cplusplus
extern "C"
{
#endif



int tdsearch_greens_writeSelectGreensFunctions(
    const char *dirnm,
    const int iobs, const int itstar, const int idepth,
    const struct tdSearchGreens_struct grns);
int tdsearch_greens_free(struct tdSearchGreens_struct *grns);
int tdsearch_greens_ffGreensToGreens(const struct tdSearchData_struct data,
                                     const struct tdSearchHudson_struct ffGrns,
                                     struct tdSearchGreens_struct *grns);
int tdsearch_greens_getGreensFunctionIndex(
    const enum tdSearchGreens_enum GMT_TERM,
    const int iobs, const int itstar, const int idepth,
    const struct tdSearchGreens_struct grns);

#ifdef __cplusplus
}
#endif

#endif
