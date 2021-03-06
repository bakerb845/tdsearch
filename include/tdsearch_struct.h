#ifndef TDSEARCH_STRUCT_H__
#define TDSEARCH_STRUCT_H__ 1
#include <sacio.h>
#include <cps.h>
#include <compearth.h>
#include "prepmt/prepmt_struct.h"

//enum tdSearchGreens_enum
//{
//    G11_GRNS = 1,  /*!< Green's function that scales xx moment tensor term */
//    G22_GRNS = 2,  /*!< Green's function that scales yy moment tensor term */
//    G33_GRNS = 3,  /*!< Green's function that scales zz moment tensor term */
//    G12_GRNS = 4,  /*!< Green's function that scales xy moment tensor term */
//    G13_GRNS = 5,  /*!< Green's function that scales xz moment tensor term */
//    G23_GRNS = 6   /*!< Green's function that scales yz moment tensor term */
//};

struct tdSearchEventParms_struct
{
    double latitude;      /*!< Event latitude (degrees). */
    double longitude;     /*!< Event longitude (degrees). */
    double depth;         /*!< Event depth (km). */
    double time;          /*!< Origin epochal time (UTC seconds). */
    double m11;           /*!< Moment tensor component for xx or rr. */
    double m22;           /*!< Moment tensor component for yy or tt. */
    double m33;           /*!< Moment tensor component for zz or pp. */
    double m12;           /*!< Moment tensor component for xy or rt. */
    double m13;           /*!< Moment tensor component for xz or rp. */
    double m23;           /*!< Moment tensor component for yz or tp. */
    enum compearthCoordSystem_enum
         basis;           /*!< Moment tensor basis - NED or USE. */
    char pad[4];
};

struct tdSearchModifyCommands_struct
{
    double targetDt;
    double cut0;
    double cut1;
    int iodva;
    bool ldeconvolution;
};

struct tdSearchGreens_struct
{
    struct sacData_struct 
        *grns;            /*!< Green's functions.  These are packed
                               in row major format with the final row
                               corresponding to the
                               \f$ \{Gxx, Gyy, Gzz, Gxy, Gxz, Gyz\}  \f$
                               Green's functions which are
                               in the NED coordinate system and scaled so that
                               they can be applied to a moment tensor specified
                               in Newton-meters [nobs x ndepth x ntstar x 6]. */
    struct tdSearchDataProcessingCommands_struct
           *cmds;         /*!< Data processing commands for each observation
                               [nobs]. */
    struct tdSearchDataProcessingCommands_struct
           *cmdsGrns;     /*!< Green's functions commands that are to be
                               applied to each Green's functions. */
    int ntstar;           /*!< Number of t*'s (should correspodn to number of
                               t*'s in grid search). */
    int ndepth;           /*!< Number of depths (should correspond to number of
                               depths in grid search). */
    int nobs;             /*!< Number of observations (should correspond to
                                data). */
    int ngrns;            /*!< Number of Green's functions. */
};

struct tdsearchHudsonParms_struct
{
    struct hudson96_parms_struct
         hudson96Parms; /*!< Hudson96 modeling paramters. */
    struct hpulse96_parms_struct
        hpulse96Parms; /*!< Hpulse96 modeling parameters. */
};

struct tdSearchHudson_struct
{
    struct tdsearchHudsonParms_struct
         modelingParms; /*!< Waveform modeling parameters */
    struct sacData_struct *grns;
    char srcModelName[PATH_MAX]; /*!< Name of source model. */
    char crust1Dir[PATH_MAX];    /*!< Name of crust1.0 directory. */
    double *tstars;     /*!< t*'s at which to perform modeling [nstar]. */
    double *depths;     /*!< Source depths (km) at which to perform modeling
                             [ndepth]  */
    //struct vmodel_struct
    //       *teleModels; /*!< Teleseismic source models (ak135). [ndists] */
    //struct vmodel_struct
    //       *recvModels; /*!< Receiver source models [ndists]. */
    //struct vmodel_struct
    //       srcModels;   /*!< Source model. */
    int ntstar;         /*!< Number of t*'s. */
    int ndepth;         /*!< Number of depths. */
    int nobs;           /*!< Number of observations. */
    bool luseCrust1;    /*!< If true then use crust1.0 models at source
                             and receiver.  Otherwise use ak135.  */
    bool luseSrcModel;  /*!< If true then use the local source model
                             (srcModelName).  This supersedes the local
                             source model if luseCrust1 is true.  If
                             luseCrust1 and luseSrcModel are both false
                             then the program will default to ak135. */ 
    
};

//struct tdSearchGridSearchParms_struct
//{
//    double *tstar;           /*!< t* in grid search. */
//    double *depths;          /*!< Depths (km) in grid search. */
//    double maxShiftTime;     /*!< Max allowable shift time allowable (s). */
//    int ntstar;              /*!< Number of t*'s. */
//    int ndepth;              /*!< Number of depths. */
//    bool luseMaxShiftTime;   /*!< If true then the max-shift time is used in
//                                  the inversion. */
//    char pad[3];
//};

struct tdSearchDataProcessingCommands_struct
{
    char **cmds;   /*!< String commands for proessing data [ncmds]. */
    int ncmds;     /*!< Number of data processing commands. */
};

struct tdSearchData_struct
{
    struct sacData_struct *obs; /*!< Observed waveforms [maxobs]. */
    struct tdSearchDataProcessingCommands_struct *cmds;
    bool *lskip; /*!< If true then skip the [iobs]'th waveform. */
    int nobs;    /*!< Number of observations. */
    int maxobs;  /*!< Maximum space allotted to obs. */
};

struct tdSearch_struct
{
    double *G;           /*!< Green's functions forward modeling matrix
                              with columns ordered according to mt. */
    double *Gxc;         /*!< Cross-correlation forward modeling matrix
                              with columns ordered accordeding to mt. */
    double *synthetic;   /*!< Synthetic seismograms. */ 
    int *lags;           /*!< Synthetic waveform lags. */
    double *xcorr;       /*!< Cross-correlation map [ntstar x ndepth]. */
    double *tstar;       /*!< Holds the t* values in grid search [ntstar]. */
    double *depths;      /*!< Holds the depths in grid search [ndepth]. */
    double mt[8] __attribute__ ((aligned(64)));
                         /*!< Holds the moment tensor in NED format packed
                              \f$ \{ m_{xx}, m_{yy}, m_{zz}, 
                                     m_{xy}, m_{xz}, m_{yz} \} \f$.  */
    int *grnsMatrixPtr;   /*!< Points from (t*,depth) to start row index
                               of G [ntstar*ndepth+1]. */
    int *grnsXCMatrixPtr; /*!< Points from (t*,depth) to start row index
                               of Gxc [ntstar*ndepth+1]. */
    double maxShiftTime;  /*!< Max acceptable shift time (seconds) in
                               xcorr grid search. */
    double dt;            /*!< Sampling period (seconds) of data. */
    double dnorm;         /*!< Two norm of the data. */
    int npts;             /*!< Number of data points. */
    int nptsGrns;        /*!< Number of points in Green's functions. */
    int ngrns;           /*!< Total number of Green's functions - should
                              be ndepths*ntstar. */
    int ntstar;          /*!< Number of t* points in grid search. */
    int ndepth;          /*!< Number of depths in grid search. */
    int idopt;           /*!< Optimum depth index. */
    int itopt;           /*!< Optimum t* index. */
    bool lhaveMT;        /*!< If true then the moment tensor has been set. */
    char pad[3];
};

#endif
