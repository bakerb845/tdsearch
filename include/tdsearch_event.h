#ifndef TDSEARCH_EVENT_H__
#define TDSEARCH_EVENT_H__ 1
#include "tdsearch_struct.h"

#ifdef __cplusplus
extern "C"
{
#endif

/*! Set the default event information */
int tdsearch_event_setDefaults(struct tdSearchEventParms_struct *event);
/*! Read the event properties from the ini file */
int tdsearch_event_initializeFromIniFile(
    const char *fname,
    struct tdSearchEventParms_struct *event);

#ifdef __cplusplus
}
#endif
#endif
