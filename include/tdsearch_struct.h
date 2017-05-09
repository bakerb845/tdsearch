#ifndef TDSEARCH_STRUCT_H__
#define TDSEARCH_STRUCT_H__ 1
#include <sacio.h>
#include <compearth.h>

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

struct tdSearchGridSearchParms_struct
{
    double *tstar;           /*!< t* in grid search. */
    double *depths;          /*!< Depths (km) in grid search. */
    double maxShiftTime;     /*!< Max allowable shift time allowable (s). */
    int ntstar;              /*!< Number of t*'s. */
    int ndepth;              /*!< Number of depths. */
    bool luseMaxShiftTime;   /*!< If true then the max-shift time is used in
                                  the inversion. */
    char pad[3];
};

struct tdSearchData_struct
{
    struct sacData_struct *obs; /*!< Observed waveforms [maxobs]. */
    int nobs;   /*!< Number of observations. */
    int maxobs; /*!< Maximum space allotted to obs. */
};

struct tdSearch_struct
{
    double *data;        /*!< Data [npts]. */
    double *Gmat;        /*!< Green's functions matrix wth columns ordered
                              according to mt [nptsGrns*ldg]. */
    double *synthetic;   /*!< Estimate seismogram [nptsGrns]. */
    double *xcorr;       /*!< Cross-correlation map [ntstar x ndepth]. */
    double *delayTime;   /*!< Delay time (seconds) to shift each Green's
                              function to the observed [ntstar x ndepth]. */
    double *tstar;       /*!< Holds the t* values in grid search [ntstar]. */
    double *depths;      /*!< Holds the depths in grid search [ndepth]. */
    double mt[8] __attribute__ ((aligned(64)));
                         /*!< Holds the moment tensor in NED format packed
                              \f$ \{ m_{xx}, m_{yy}, m_{zz}, 
                                     m_{xy}, m_{xz}, m_{yz} \} \f$.  */
    double maxShiftTime; /*!< Max acceptable shift time (seconds) in
                              xcorr grid search. */
    double dt;           /*!< Sampling period (seconds) of data. */
    int npts;            /*!< Number of data points. */
    int nptsGrns;        /*!< Number of points in Green's functions. */
    int ldg;             /*!< Leading dimension of Gmat (should be 8). */
    int mptsGrns;        /*!< Leading dimension of (>= nptsGrns). */
    int nrowsGmat;       /*!< Number of rows in Gmat (=mptsGrns). */
    int ncolsGmat;       /*!< Number of columns in Gmat (should be 6). */
    int ngrns;           /*!< Total number of Green's functions - should
                              be ndepths*ntstar. */
    int ntstar;          /*!< Number of t* points in grid search. */
    int ndepth;          /*!< Number of depths in grid search. */
    int idopt;           /*!< Optimum depth index. */
    int itopt;           /*!< Optimum t* index. */
    bool lhaveData;      /*!< If true then the data has been set. */
    bool lhaveGrns;      /*!< If true then the Green's functions are set. */
    bool lhaveMT;        /*!< If true then the moment tensor has been set. */
    bool lhaveTD;        /*!< If true then the t* and depths defining the
                              grid search were. */
    bool lhaveMaxShiftTime; /*!< If true then the max delay time was set. */
    bool lrecompute;        /*!< If true then the grid search must be
                                 recomputed. */
    char pad[4];
};

#endif
