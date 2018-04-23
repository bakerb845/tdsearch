#ifndef TDSEARCH_GREENS_H__
#define TDSEARCH_GREENS_H__ 1
#include "tdsearch_struct.h"

/*!
 * @brief Green's function configuration and generation.
 * @defgroup tdsearch_greens Green's Functions
 * @ingroup tdsearch
 */

#ifdef __cplusplus
extern "C"
{
#endif


int tdsearch_greens_setPreprocessingCommandsFromIniFile(
    const char *iniFile,
    const int nobs,  
    struct tdSearchGreens_struct *grns);
int tdsearch_greens_modifyProcessingCommands(
    const int iodva,
    const double cut0, const double cut1, const double targetDt,
    struct tdSearchGreens_struct *grns);
int tdsearch_greens_repickGreensWithSTALTA(
    const double sta, const double lta, const double threshPct,
    const int iobs, const int itstar, const int idepth,
    struct tdSearchGreens_struct *grns);
int tdsearch_greens_getGreensFunctionsIndices(
    const int iobs, const int itstar, const int idepth,
    const struct tdSearchGreens_struct grns, int indices[6]);
int tdsearch_greens_attachCommandsToGreens(const int iobs, const int ncmds,
                                           const char **cmds,
                                           struct tdSearchGreens_struct *grns);
int tdsearch_greens_process(struct tdSearchGreens_struct *grns);
int tdsearch_greens_writeSelectGreensFunctions(
    const char *dirnm,
    const int iobs, const int itstar, const int idepth,
    const struct tdSearchGreens_struct grns);
int tdsearch_greens_attachCommandsToGreens(const int iobs, const int ncmds,
                                           const char **cmds,
                                           struct tdSearchGreens_struct *grns);
int tdsearch_greens_free(struct tdSearchGreens_struct *grns);
int tdsearch_greens_ffGreensToGreens(const struct tdSearchData_struct data,
                                     const struct tdSearchHudson_struct ffGrns,
                                     struct tdSearchGreens_struct *grns);
int tdsearch_greens_getGreensFunctionIndex(
    const enum prepmtGreens_enum GMT_TERM,
    const int iobs, const int itstar, const int idepth,
    const struct tdSearchGreens_struct grns);

#ifdef __cplusplus
}
#endif

#endif
