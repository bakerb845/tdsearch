#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iniparser.h>
#include "tdsearch_hudson.h"
#include "cps.h"
#include "iscl/log/log.h"
#include "iscl/os/os.h"

/*!
 * @brief Reads the ini file for Computer Programs in Seismolgoy hudson96
 *        forward modeling variables.
 *
 * @param[in] ini_file    Name of ini file.
 *
 * @param[out] parms      The hudson96 parameters.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_readHudson96Parameters(const char *iniFile,
                                           struct hudson96_parms_struct *parms)
{
    const char *fcnm = "tdsearch_hudson_readHudson96Parameters\0";
    const char *s;
    dictionary *ini;
    //------------------------------------------------------------------------//
    cps_setHudson96Defaults(parms);
    if (!os_path_isfile(iniFile))
    {
        log_errorF("%s: ini file: %s does not exist\n", fcnm, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    if (ini == NULL)
    {
        log_errorF("%s: Cannot parse ini file\n", fcnm);
        return -1;
    }
    // Teleseismic model
    s = iniparser_getstring(ini, "hudson96:modeltel", NULL);
    if (s != NULL)
    {
        strcpy(parms->modeltel, s);
    }
    else
    {
        strcpy(parms->modeltel, "tak135sph.mod\0");
    }
    // Receiver model
    s = iniparser_getstring(ini, "hudson96:modelrec", NULL);
    if (s != NULL)
    {
        strcpy(parms->modelrec, s);
    }
    else
    {
        strcpy(parms->modelrec, parms->modelrec);
    }
    // Source model
    s = iniparser_getstring(ini, "hudson96:modelsrc", NULL);
    if (s != NULL)
    {
        strcpy(parms->modelsrc, s);
    }
    else
    {
        strcpy(parms->modelsrc, parms->modelsrc);
    }
    parms->hs = iniparser_getdouble(ini, "hudson96:hs", 0.0);
    parms->dt = iniparser_getdouble(ini, "hudson96:dt", 1.0);
    parms->npts = iniparser_getint(ini, "hudson96:npts", 1024);
    parms->gcarc = iniparser_getdouble(ini, "hudson96:gcarc", 50.0);
    parms->offset = iniparser_getdouble(ini, "hudson96:offset", 10.0);
    parms->dosrc = iniparser_getboolean(ini, "hudson96:dosrc", 1);
    parms->dorec = iniparser_getboolean(ini, "hudson96:dorec", 1);
    parms->dotel = iniparser_getboolean(ini, "hudson96:dotel", 1);
    parms->dop = iniparser_getboolean(ini, "hudson96:dop", 1);
    parms->dokjar = iniparser_getboolean(ini, "hudson96:dokjar", 1);
    parms->loffsetdefault = iniparser_getboolean(ini,
                                                 "hudson96:loffsetdefault", 1);
    parms->verbose = iniparser_getboolean(ini, "hudson96:verbose", 0);
    parms->utstar = iniparser_getdouble(ini, "hudson96:utstar", -12345.0);
    parms->dottonly = iniparser_getboolean(ini, "hudson96:dottonly", 0);
    parms->zsrc = iniparser_getdouble(ini, "hudson96:zsrc", 100.0);
    parms->zrec = iniparser_getdouble(ini, "hudson96:zrec", 60.0);
    // Free ini dictionary
    iniparser_freedict(ini);
    return 0;
}
//============================================================================//
/*!
 * @brief Reads the ini file for the Computer Programs in Seismology hpulse96
 *        modeling variables.
 *
 * @param[in] iniFile     Name of ini file.
 *
 * @param[out] parms      The hpulse96 parameters.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_readHpulse96Parameters(const char *iniFile,
                                           struct hpulse96_parms_struct *parms)
{
    const char *fcnm = "readini_hpulse\0";
    const char *s;
    dictionary *ini;
    //------------------------------------------------------------------------//
    cps_setHpulse96Defaults(parms);
    if (!os_path_isfile(iniFile))
    {
        log_errorF("%s: ini file: %s does not exist\n", fcnm, iniFile);
        return -1;
    }
    ini = iniparser_load(iniFile);
    if (ini == NULL)
    {
        log_errorF("%s: Cannot parse ini file\n", fcnm);
        return -1;
    }
    parms->ntau = iniparser_getint(ini, "hpulse96:ntau", -1);
    parms->alp = iniparser_getdouble(ini, "hpulse96:alp", -1.0);
    parms->ipt = iniparser_getint(ini, "hpulse96:ipt", -1);
    if (parms->ipt == 4)
    {
        s = iniparser_getstring(ini, "hpulse96:rfile", NULL);
        if (s != NULL)
        {
            strcpy(parms->rfile, s);
            if (!os_path_isfile(parms->rfile))
            {
                log_errorF("%s: Response file does not exist!\n", fcnm);
                return -1;
            }
        }
        else
        {
            log_errorF("%s: Response file not specified!\n", fcnm);
            return -1;
        }
    }
    parms->idva = iniparser_getint(ini, "hpulse96:idva", 1);
    parms->iodva = iniparser_getint(ini, "hpulse96:iodva", -1);
    parms->xmult = iniparser_getdouble(ini, "hpulse96:xmult", 1.0);
    parms->dozero = iniparser_getboolean(ini, "hpulse96:dozero", 0);

    if (parms->ipt >= 0 && parms->ipt <= 1 && parms->ntau <= 0)
    {
        parms->ntau = 1;
    }
    if (parms->ipt == 2 && parms->alp <= 0.0){parms->alp = 1.0;}
    if (parms->ipt == 2 && parms->alp < 0.0)
    {
        log_errorF("%s: No alpha for Ohnaka pulse\n", fcnm);
        return -1;
    }
    if (parms->ipt < 0)
    {
        log_errorF("%s: No pulse shape defined\n", fcnm);
        return -1;
    }
    if (parms->iodva < 0){parms->iodva = parms->idva;}
    iniparser_freedict(ini);
    return 0;
}
//============================================================================//
/*!
 * @brief Reads the hpulse and hudson forward modeling parameters from the
 *        ini file.
 *
 * @param[in] iniFile    Name of ini file to read.
 *
 * @param[out] grns      On successful exit contains the forward modeling
 *                       parameters for hudson and hpulse. 
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_initializeParametersFromIniFile(
    const char *iniFile,
    struct tdSearchGreens_struct *grns)
{
    const char *fcnm = "tdsearch_hudson_initializeParametersFromIniFile\0";
    struct hudson96_parms_struct hudsonParms;
    struct hpulse96_parms_struct hpulseParms;
    int ierr, ierr1;
    ierr = 0;
    memset(grns, 0, sizeof(struct tdSearchGreens_struct));
    ierr1 = tdsearch_greens_readHudson96Parameters(iniFile, &hudsonParms);
    if (ierr1 != 0)
    {
        log_errorF("%s: Error reading hudson96 parameters\n", fcnm);
        ierr = ierr + 1;
    }
    ierr1 = tdsearch_greens_readHpulse96Parameters(iniFile, &hpulseParms);
    if (ierr1 != 0)
    {   
        log_errorF("%s: Error reading hpulse96 parameters\n", fcnm);
        ierr = ierr + 1;
    }
    // Set the modeling parameters on the modeling structure anyway
    tdsearch_greens_setHudson96Parms(hudsonParms, grns);
    tdsearch_greens_setHpulse96Parms(hpulseParms, grns);
    return ierr;
}
//============================================================================//
/*!
 * @brief Convenience utility for setting the hudson forward modeling
 *        parameters on the greens structure.
 *
 * @param[in] hudsonParms  The hudson96 forward modeling parameters.
 *
 * @param[out] grns        Contains the hudson forward modeling parameters.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_setHudson96Parms(
    const struct hudson96_parms_struct hudsonParms,
    struct tdSearchGreens_struct *grns)
{
    cps_utils_copyHudson96ParmsStruct(hudsonParms,
                                      &grns->modelingParms.hudson96Parms);
    return 0;
}
//============================================================================//
/*!
 * @brief Convenience utility for setting the hpulse source time function
 *        parameters on the greens structure.
 *
 * @param[in] hpulseParms  The hpulse96 source time function parameters.
 *
 * @param[out] grns        Contains the hpulse STF parameters.
 *
 * @result 0 indicates success.
 *
 * @author Ben Baker, ISTI
 *
 */
int tdsearch_greens_setHpulse96Parms(
    const struct hpulse96_parms_struct hpulseParms,
    struct tdSearchGreens_struct *grns) 
{
    cps_utils_copyHpulse96ParmsStruct(hpulseParms,
                                      &grns->modelingParms.hpulse96Parms);
    return 0;
}
//============================================================================//
int tdsearch_greens_computeGreens( )
{
    int ierr;

    return 0;
}
