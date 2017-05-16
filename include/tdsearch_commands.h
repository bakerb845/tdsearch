#ifndef TDSEARCH_COMMANDS_H__
#define TDSEARCH_COMMANDS_H__ 1
#include "tdsearch_struct.h"

#ifdef __cplusplus
extern "C"
{
#endif

char **tdsearch_commands_modifyCommands(
    const int ncmds, const char **cmds,
    const struct tdSearchModifyCommands_struct options,
    const struct sacData_struct data,
    int *ierr);

#ifdef __cplusplus
}
#endif
#endif
