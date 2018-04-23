#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <iniparser.h>
#include "tdsearch_event.h"
#include "iscl/os/os.h"
#include "iscl/time/time.h"

#define M11 1.0        /*! Default to explosion - mxx or mtt */
#define M22 1.0        /*! Default to explosion - myy or mpp */
#define M33 1.0        /*! Default to explosion - mzz or mrr */
#define M12 0.0        /*! Default to explosion - */
#define M13 0.0        /*! Default to explosion - */
#define M23 0.0        /*! Default to explosion - */
#define EVLA 41.30800  /*! Default to North Korea test site latitude */
#define EVLO 129.0760  /*! Default to North Korea test site longitude */
#define EVDP 1.0       /*! Default to shallow depth */
#define EVTIME 0.0     /*! Default origin time is unknown */

/*!
 * @brief Sets the default event information for North Korea.
 *
 * @param[out] event   Default event is an explosion in North Korea.
 *
 * @result 0 indicates success.
 *
 * @ingroup tdsearch_event
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_event_setDefaults(struct tdSearchEventParms_struct *event)
{
    memset(event, 0, sizeof(struct tdSearchEventParms_struct));
    event->m11 = M11;
    event->m22 = M22;
    event->m33 = M33;
    event->m12 = M12;
    event->m13 = M13;
    event->m23 = M23;
    event->basis = CE_NED;
    event->latitude = EVLA;
    event->longitude = EVLO;
    event->depth = EVDP;
    event->time = EVTIME; 
    return 0;
}
//============================================================================//
/*!
 * @brief Reads the event information from the ini file.
 *
 * @param[in] fname    Name of ini file.
 *
 * @param[out] event   Event information (lat, lon, depth, origin time, and
 *                     moment tensor).
 *
 * @result 0 indicates success.
 *
 * @ingroup tdsearch_event
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_event_initializeFromIniFile(
    const char *fname,
    struct tdSearchEventParms_struct *event)
{
    const char *s;
    double second;
    int dom, month, nzhour, nzmin, nzsec, nzmusec, nzyear;
    bool lread;
    dictionary *ini;
    //------------------------------------------------------------------------//
    //
    // Set the defaults
    tdsearch_event_setDefaults(event);
    // Verify the file exists
    if (!os_path_isfile(fname))
    {
        fprintf(stderr, "%s: Error ini file %s doesn't exist\n",
                __func__, fname);
        return -1;
    }
    ini = iniparser_load(fname);
    // Read the moment tensor terms
    lread = false;
    s = iniparser_getstring(ini, "eventInfo:basis\0", NULL);
    if (s != NULL)
    {
        if (strcasecmp(s, "USE\0") == 0)
        {
            lread = true;
            event->basis = CE_USE;
        }
        if (strcasecmp(s, "NED\0") == 0)
        {
            lread = true;
            event->basis = CE_NED;
        }
    }
    if (lread && event->basis == CE_USE)
    {
        event->m11 = iniparser_getdouble(ini, "eventInfo:mrr\0", M11);
        event->m22 = iniparser_getdouble(ini, "eventInfo:mtt\0", M22);
        event->m33 = iniparser_getdouble(ini, "eventInfo:mpp\0", M33);
        event->m12 = iniparser_getdouble(ini, "eventInfo:mrt\0", M12);
        event->m13 = iniparser_getdouble(ini, "eventInfo:mrp\0", M13);
        event->m23 = iniparser_getdouble(ini, "eventInfo:mtp\0", M23);
    }
    if (lread && event->basis == CE_NED)
    {
        event->m11 = iniparser_getdouble(ini, "eventInfo:mxx\0", M11);
        event->m22 = iniparser_getdouble(ini, "eventInfo:myy\0", M22);
        event->m33 = iniparser_getdouble(ini, "eventInfo:mzz\0", M33);
        event->m12 = iniparser_getdouble(ini, "eventInfo:mxy\0", M12);
        event->m13 = iniparser_getdouble(ini, "eventInfo:mxz\0", M13);
        event->m23 = iniparser_getdouble(ini, "eventInfo:myz\0", M23);
    }
    event->latitude  = iniparser_getdouble(ini, "eventInfo:latitude\0",  EVLA);
    event->longitude = iniparser_getdouble(ini, "eventInfo:longitude\0", EVLO);
    event->depth     = iniparser_getdouble(ini, "eventInfo:depth\0",     EVDP);
    s = iniparser_getstring(ini, "eventInfo:time\0", NULL);
    if (s != NULL)
    {
        sscanf(s, "%d-%d-%d:%d:%d:%lf",
               &nzyear, &month, &dom, &nzhour, &nzmin, &second);
        nzsec = (int) second;
        nzmusec = (int) ((second - (double) nzsec)*1.e6);
        event->time = time_calendar2epoch2(nzyear, month, dom, nzhour,
                                           nzmin, nzsec, nzmusec);
    }
    iniparser_freedict(ini);
    return 0;
}
