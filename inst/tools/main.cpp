/*####################################################################
 ###  TROLL
 ###  Individual-based forest dynamics simulator
 ###    Version 1: Jerome Chave
 ###    Version 2: Isabelle Marechaux & Jerome Chave
 ###
 ###  History:
 ###    version 0.1 --- JC - 22 Sep 97
 ###    version 0.2 --- JC - 06 Oct 97
 ###    version 0.3 --- JC - 11-14 Nov 97
 ###    version 1.0 --- JC - stable version Chave, Ecological Modelling (1999)
 ###    version 1.1 --- JC - 02-30 Sep 98
 ###    version 1.2 --- JC - 22 Jan 00
 ###	version 1.3 --- JC - 28 Sep 01 stable version Chave, American Naturalist (2001)
 ###
 ###	version 2.0 --- JC - 23 Mar 11 (physiology-based version, translation of comments into English)
 ###    version 2.01 --- IM - oct-nov 13
 ###    version 2.02 --- IM - apr-may 2015
 ###    version 2.03 --- IM - jul 2015
 ###    version 2.04 --- IM - jul 2015 (monthly timestep)
 ###    version 2.1 --- IM - dec 2015 (undef/defined alternative versions)
 ###    version 2.11 --- IM - jan 2016 (timestep better used and two input files - one for species, one for climate and environment parameters)
 ###    version 2.12 --- JC - jan 2016 porting to GitHub for social coding, check of the MPI routines and update, new header for code, trivia: reindentation (orphan lines removed)
 ###    version 2.2 --- IM - may 2016 core changes in: daily coupling with environment; respiration; treefall module
 ###    version dev --- SS - febr 2017 disturbance module implementation
 ###
 ####################################################################*/


/*
 Glossary: MPI = Message Passing Interface. Software for sharing information
 across processors in parallel computers. If global variable MPI is not defined,
 TROLL functions on one processor only.
 */

#undef MPI              /* If flag MPI defined, parallel routines (MPI software) are switched on */
#undef easyMPI         /* If flag easyMPI defined, parallel routine for fast search of parameter space are switched on */
#undef outputMPI         /* Reduced set of ouput files (not necessarily used with the MPI option) */

/* two parameters below describe two possible light treatment inside the canopy */
#define ONLYVERTICAL	/* if defined: overhead light vertical only, otherwise diffusion */
#undef INHOMOGENE       /* initial conditions */

#define BASICTREEFALL    /* if defined: treefall is a source of tree death (and if TREEFALL not defined, this is modeled through simple comparison between tree height and a threshold t_Ct, if not defined, treefall is not represented as a separated and independent source of death, but instead, all tree death are due to the deathrate value */
#undef TREEFALL         /* computation of the force field if TREEFALL is defined, neighboring trees contribute to fell each tree */

#define DAILYLIGHT      /* if defined: the rate of carbon assimilation integrates an average daily fluctuation of light (thanks to GPPDaily). Should be defined to ensure an appropriate use of Farquhar model */
#define CROWNGRADIENT   /* if defined: the gradient of light within a crown depth is taken into account (otherwise, the light is supposed constant through the crown depth, equal to the incident light at the top of the tree) */
#undef SEEDTRADEOFF     /* if defined: the number of seeds produced by each tree is determined by the tree NPP allocated to reproduction and the species seed mass, otherwise the number of seeds is fixed; besides, seedling recruitment in one site is not made by randomly and 'equiprobably' picking one species among the seeds present at that site but the probability of recruitment among the present seeds is proportional to the number of seeds (in s_Seed[site]) time the seed mass of each species */
#undef NDD              /*if defined, negative density dependant processes affect both the probablilty of seedling recruitment and the local tree death rate. The term of density-dependance is computed as the sum of conspecific tree basal area divided by their distance to the focal tree within a neighbourhood (circle of radius 15m) */

#define DISTURBANCE      /* if defined: implementation of a basic perturbance module at a given iteration step in waiting for a more sofisticated sylviculture module */

/* Libraries */
# include <cstdio>
# include <iostream>
# include <fstream>
# include <cstdlib>
# include <string>
# include <limits>
# include <ctime>
# include <cmath>
# ifdef MPI
# include "mpi.h"
# endif
# ifdef easyMPI
# include "mpi.h"
# endif

using namespace std;

/* Global constants */
# define Pis2 1.570796327
# define PI 3.141592654
# define deuPi 6.2831853071
# define iPi 0.3183099         // 1/pi
# define mask 255
# define PERSNB 500
int col0, row0, col, row, xx, yy, dm, haut;
char buffer[256], nomfi[256], nomf0[256], *bufi(0), *buf(0);


/* random number generators */
double genrand2(void);
void sgenrand2(unsigned long);
unsigned long genrand2i(void);
void sgenrand2i(unsigned long);


#define NBFOUT 158      // Number of output files (nb: this number is currently large because of a large number of diagnostic output files, generated to track errors or model behavior, and to generate a regeneration animation movie, but can be decreased.)


/* file output streams */

fstream out,out2;
fstream sor[NBFOUT+1], sch[NBFOUT+1];


/********************************/
/* Parameters of the simulation */
/********************************/

int sites,          /* number of sites */
cols,               /* nb of columns */
rows,               /* nb of lines */
numesp,             /* nb of species */
iterperyear,        /* nb of iter in a year (=12 if monthly timestep, =365 if daily timestep) */
nbiter,             /* total nb of timesteps */
iter,               /* current timestep */
nbout,              /* nb of outputs */
freqout;            /* frequency HDF outputs */

int HEIGHT,         /* max height (in m) */
dbhmaxincm,         /* max DBH times 100 (ie dbh in cm *100 = in meters) */
RMAX,               /* max crown radius */
SBORD;              /* RMAX*cols */

float NV,           /* nb cells per m (vertical) */
NH,                 /* nb cells per m (horizontal) */
LV,                 /* LV = 1.0/NV */
LH,                 /* LH = 1.0/NH */
timestep;           /* duration of one timestep (in years)=1/iterperyear */

float p,            /* ratio of non-vertical incident light */
Cseedrain,          /* constant used to scale total seed rain per hectare across species */
nbs0,               /* nb of seeds produced and dispersed by each mature tree when SEEDTRADEOFF is not defined */
Cair;               /* atmosphericCO2 concentration, if we aim at making CO2 vary (scenarios), CO2 will have to have the same status as other climatic variables  */

/* new version 2.2 */
float daily_light[24];    /* normalized (ie between 0 and 1) daily light variation (used if DAILYLIGHT defined) */
float daily_vpd[24];      /* normalized (ie between 0 and 1) daily vpd variation (used if DAILYLIGHT defined) */
float daily_T[24];        /* normalized (ie between 0 and 1) daily T variation (used if DAILYLIGHTdefined) */

/* new disturbance module */
int disturb_iter;       /* iteration step where the disturbation occure */
float disturb_intensity;      /* intensity of disturbance in percent of BA */


/*********************************************/
/* Environmental variables of the simulation */
/*********************************************/

/* Climate input data  */
/* these climate input data are given in the input file, its structure depends on the timestep and scenario used for the simulation */
/* new version 2.2 */
float *Temperature(0);                      /* in degree Celsius */
float *DailyMaxTemperature(0);              /* in degree Celsius */
float *NightTemperature(0);                 /* in degree Celsius */
float *Rainfall(0);                         /* in mm */
float *WindSpeed(0);                        /* in m/s */
float *MaxIrradiance(0);                    /* in W/m2 */
float *MeanIrradiance(0);                   /* in W/m2 */
float *SaturatedVapourPressure(0);          /* in kPa */
float *VapourPressure(0);                   /* in kPa */
float *VapourPressureDeficit(0);            /* in kPa */
float *DailyVapourPressureDeficit(0);       /* in kPa */
float *DailyMaxVapourPressureDeficit(0);    /* in kPa */


/***** Environmental variables, changed at each timestep *****/
float temp,         /* Temperature */
tmax,               /* Daily max temperature */
tnight,             /* Night mean temperature */
precip,             /* Rainfall  */
WS,                 /* WindSpeed */
Wmax,               /* Daily max irradiance (average for timestep) (in micromol PAR photon/m^2/s)*/
/* used in the photosynthesis part. see if it would not be better to have value in the right unit in the input file, however W/m2 is the common unit of meteo station */
/* below: new version 2.2 */
Wmean,              /* mean irradiance (in W/m2)*/
e_s,                /* SaturatedVapourPressure */
e_a,                /* VapourPressure*/
VPD,                /* VapourPressureDeficit */
VPDday,             /* DailyVapourPressureDeficit */
VPDmax;             /* DailyMaxVapourPressureDeficit */


/****************************************/
/*  Common variables for the species    */
/* (simplifies initial version 170199)  */
/****************************************/

//int Cm;           /* Treefall threshold */
float klight,       /* light absorption rate or extinction coefficient used in Beer-Lambert law to compute light within the canopy */
phi,                /* apparent quantum yield (in micromol C/micromol photon). phi is the quantum yield multiplied by leaf absorbance (in the literature, the quantum yield is often given per absorbed light flux, so one should multiply incident PPFD by leaf absorbance, see Poorter et al American Journal of Botany (assumed to be 0.91 for tropical tree species). Even though holding phi constant across species is widely assumed, in reality, phi varies across species and environmental conditions (see eg Domingues et al 2014 Plant Ecology & Diversity) */
alpha,              /* apparent quantum yield to electron transport in mol e-/mol photons, equal to the true quantum yield multiplied by the leaf absorbance -- v.2.2 */
vC,                 /* variance of treefall threshold */
H0,                 /* initial height (in m) */
DBH0,               /* initial DBH (in m) */
de0,                /* initial crown Crown_Depth (in m) */
de1,                /* Crown_Depth/height slope */
/* fallocwood and falloccanopy new -- v.2.2 */
fallocwood,         /* fraction of biomass allocated to above ground wood (branches+stem) */
falloccanopy,       /* fraction of biomass allocated to canopy (leaves + reproductive organs + twigs) */
dens,               /* initial leaf density (in m^2/m^2) *//*in m^2/m^3 no?*/
ra1,                /* crown radius - dbh slope */
ra0,                /* initial crown radius (in m) */
m,                  /* basal death rate */
m1;                 /* deathrate-wsg slope -- new v.2.2 */


float **LAI3D(0);   /* leaf density (per volume unit) */
unsigned short *Thurt[3];            /* Treefall field */

int    *ESP_N (0);
#ifdef SEEDTRADEOFF
float  *PROB_S (0);
#endif
float tempRday;      /* temporary variable used for the computation of leaf day respiration -- new v.2.2 */


/***************/
/* Diagnostics */
/***************/

int nbm_n1,             /* nb deaths other than treefall dbh > 1 cm, computed at each timestep */
nbm_n10,                /* nb deaths other than treefall dbh > 10 cm, computed at each timestep */
/* nb: BASICTREEFALL new in v.2.2 */
#ifdef BASICTREEFALL
nbm_c1,                 /* nb deaths caused by a treefall dbh > 1 cm, computed at each timestep*/
nbm_c10,                /* nb deaths caused by a treefall dbh > 10 cm, computed at each timestep */
#endif
nblivetrees;            /* nb live trees for each timestep  */
#ifdef BASICTREEFALL
int nbTreefall1,            /* nb treefalls for each timestep (dbh > 1cm)*/
nbTreefall10;           /* nb treefalls for each timestep (dbh > 10 cm)*/
#endif
//long int *persist;	/* persistence histogram */
int *nbdbh(0);          /* dbh size distribution */
float *layer(0);        /* vertical LAI histogram */


/**************/
/* Processors */
/**************/

int mpi_rank,mpi_size;
int easympi_rank,easympi_size;


/******************/
/* MPI procedures */
/******************/

#ifdef MPI
unsigned short **LAIc[2];
void
MPI_ShareSeed(unsigned char **,int),
MPI_ShareField(unsigned short **,unsigned short ***,int),
MPI_ShareTreefall(unsigned short **,int);
#endif


/**********************/
/* Simulator routines */
/**********************/

void
Initialise(void),
AllocMem(void),
BirthInit(void),
Evolution(void),
Disturbance(void), /* new disturbance module */ 
UpdateField(void),
#ifdef BASICTREEFALL /* new v.2.2 */
UpdateTreefall(void),
#endif
UpdateTree(void),
Average(void),
OutputField(void),
//  SorCh(unsigned short*, int, int, fstream &),
FreeMem(void);


/****************************/
/* Various inline functions */
/****************************/


inline float flor(float f) {
    if(f>0.) return f;
    else return 0.;
}
inline float florif(int i) {
    if(i>0) return float(i);
    else return 0.;
}
inline float maxf(float f1, float f2) {
    if(f1>f2) return f1;
    else return f2;
}
inline float minf(float f1, float f2) {
    if(f1<f2) return f1;
    else return f2;
}
inline int min(int i1, int i2) {
    if(i1<i2) return i1;
    else return i2;
}
inline int max(int i1, int i2) {
    if(i1>i2) return i1;
    else return i2;
}
inline int sgn(float f) {
    if(f>0.0) return 1;
    else return -1;
}


/*############################################
 ############################################
 ###########     Species  class   ###########
 ############################################
 ############################################*/

class Species {
    
public:
    int    s_nbind,			/* nb of individuals per species */
    s_dormDuration,         /* seed dormancy duration -- not used in v.2.2 */
    //s_freqgerm,           /* germination frequency -- deprecated v.2.2, we now assume seed production at each timestep */
    s_nbext;                /* total number of incoming seeds in the simulated plot at each timestep (seed rain) -- v.2.2 */
    char	s_name[256];	/* species name */
    float  s_LCP,			/* light compensation point  (in micromol photon/m^2/s) */
    s_Rdark,                /* dark respiration rate (at PPFD = 0) in micromol C/m^2/s) */
    s_ds,                   /* average dispersal distance */
    //  de1,                /* (crown depth) - height slope deprecated v.2.1 */
    s_dmax,                 /* maximal dbh (given in m) */
    s_hmax,                 /* maximal height (given in m) */
    s_Vcmax,                /* maximal rate of carboxylation, on an area basis (in micromolC/m^2/s) */
    s_Vcmaxm,               /* maximal rate of carboxylation, on a mass basis (in micromolC/g-1/s) */
    s_Jmax,                 /* maximal rate of electron transport, on an area basis (in micromol/m^2/s) */
    s_g1,                   /* parameter of the stomatal conductance model of Medlyn et al 2011 -- V.2.2 */
    s_fci,                  /* fraction of CO2 partial pressure in intercellular spaces divided by ambiant CO2 partial pressure (both in microbar, or ppm = micromol/mol) */
    s_Gamma,                /* compensation point for the carboxylation rate, here NORMALIZED by atm CO2 concentration (Cair) */
    s_Km,                   /* apparent kinetic constant for the rubiscco = Kc*(1+[O]/Ko), here normalized by atm CO2 concentration (Cair) */
    //s_d13C,                 /* isotopic carbon discrimination NOW normalized at zero height -- deprecated v.2.2 */
    s_LMA,                  /* leaf mass per area (in g/m^2) */
    s_Nmass,                /* leaf nitrogen concentration (in g/g) v.2.01 */
    s_Pmass,                /* leaf phosphorous concentration (in g/g) v.2.01 */
    s_wsg,                  /* wood specific gravity (in g/cm^3) */
    s_ah,                   /* parameter for allometric height-dbh equation */
    //s_seedmass,             /* seed mass, in g (from Baraloto & Forget 2007 dataset, in classes) deprecated v.2.2 */
    //s_iseedmass,            /* inverse of seed mass deprecated v.2.2 */
    //s_factord13Ca, /* deprecated v.2.2 -- factor used for a previous version of ci/ca ratio computation, from d13C value */
    //s_factord13Cb,        /* deprecated v.2.2 -- factor used for a previous version of ci/ca ratio computation, from d13C value */
    /* Below: new in v.2.2 */
    s_leaflifespan,         /* average leaf lifespan, in month */
    s_time_young,           /* leaf resident time in the young leaf class */
    s_time_mature,          /* leaf resident time in the mature leaf class */
    s_time_old,             /* leaf resident time in the old leaf class */
    s_chsor[24];            /* scalar output fields per species (<24) */
    
    
    unsigned char *s_Seed;	/* presence/absence of seeds; if def SEEDTRADEOFF, the number of seeds */
#ifdef MPI
    unsigned char *s_Gc[4]; /* MPI: seeds on neighboring procs */
#endif
    
    
    Species() {
        s_nbind=0;
        s_Seed=0;
    };                              /* constructor */
    
    virtual ~Species() {
        delete [] s_Seed;
    };                              /* destructor */
    
    void Init(int,fstream&);        /* init Species class */
    void FillSeed(int,int);         /* fills s_Seed field (and s_Gc (MPI)) */
    void UpdateSeed(int,int);       /* Updates s_Seed field */
#ifdef MPI
    void AddSeed(void);             /* MPI: adds fields s_Gc  to field s_Seed */
#endif
#ifdef NDD
    inline float DeathRate(float, float, float);
#else
    inline float DeathRate(float, float, int);  /* actual death rate -- new argument int v.2.2 */
#endif
    inline float GPPleaf(float, float, float);    /* Computation of the light-limited leaf-level NPP per m^2 (in micromol/m^2/s) -- two new arguments float v.2.2 */
    /* Farquhar von Caemmerer Berry model */
    inline float NPP(float, float, float); /* two new arguments float v.2.2 */
#ifdef DAILYLIGHT
    inline float dailyGPPleaf(float, float, float);    /* computation of the daily average assimilation rate, taking into account the daily variation in light, VPD and temperature two new arguments float v.2.2 */
#endif
};


/*############################
 ###  Initialize Species  ###
 ###    Species::Init     ###
 ############################*/

void Species::Init(int nesp,fstream& is) {
    
    int site;
    float regionalfreq;     // "regional" relative abundance of the species -- new name v.2.2
    float SLA;              // specific leaf area = 1/LMA
    
    /*** Read parameters ***/
    
    is  >> s_name >> s_Nmass  >> s_LMA >> s_wsg >> s_dmax >> s_hmax  >> s_ah
    >> s_dormDuration >> regionalfreq >> s_ds >>  s_Pmass >>  s_g1 ; /* Narea => Nmass v.2.0 */
    
    s_fci = 0.0;
    /* start changed in v.2.2 */
    /* Module of leaf lifespan as a function of functional traits */
    s_g1=3.77;
    //s_g1=-3.97*s_wsg+6.53;
    //s_iseedmass=1/s_seedmass;
    s_nbext = int(regionalfreq*Cseedrain*(sites*LH*LH/10000));
    //s_nbext = int(0.05*Cseedrain*(sites*LH*LH/10000));
    //s_freqgerm=1;           // germination frequency (in number of timesteps)
    
    //s_leaflifespan=0.5+pow(10,(-2.509+1.71*log10(s_LMA)));                  //this is the expression from Wright et al 2004 Nature (leaf economics spectrum).
    s_leaflifespan=1.5+pow(10,(7.18+3.03*log10(s_LMA*0.0001)));           //this is the expression from Reich et al 1991 Oecologia (San Carlos Rio Negro).
    
    s_time_young=1;
    s_time_mature=s_leaflifespan/3;
    s_time_old=s_leaflifespan-s_time_mature-s_time_young;
    /* End changed in v.2.2 */
    
    /*** Normalization of the parameters ***/
    /* vertical (NV) and horizontal (NH) scales */
    
    s_ah *= NV*LH;
    //if(!p_rank && (H0 < 1.0)) cerr << " probleme dans la routine TREEFALL\n";
    s_ds *= NH;
    s_hmax *= NV;
    s_dmax *= NH;
    
    s_nbind=0;
    
    SLA=10000/s_LMA;   // computation of specific leaf area in cm^2/g for use in Domingues et al 2010 model of photosynthetic capacities
    
    s_Vcmaxm=pow(10, minf((-1.56+0.43*log10(s_Nmass*1000)+0.37*log10(SLA)), (-0.80+0.45*log10(s_Pmass*1000)+0.25*log10(SLA))));    // this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7) which made better fits than equ1 (without LMA)
    s_Vcmax=s_Vcmaxm*s_LMA;
    //s_Vcmax=23.0*s_Narea-7.02;     // Relationship from Domingues et al 2005
    
    s_Jmax=pow(10, minf((-1.50+0.41*log10(s_Nmass*1000)+0.45*log10(SLA)), (-0.74+0.44*log10(s_Pmass*1000)+0.32*log10(SLA))))*s_LMA;         /* added as a Species member variable 14-04-2015; this is equ 2 in Domingues et al 2010 PCE (coefficients from fig7) */
    //s_Jmax=1.67*s_Vcmax;        /* this is Medlyn et al 2002 Jmax-Vcmax correlation */
    
    //s_Km = 565/Cair;          /* this value correponds to the formula Kc*(1+O/Ko), with the values: Kc=260microbar and Ko=179mbar (von caemmerer 2000 page 45) and O=210mbar (Farquhar et al 1980). deprecated v.2.2 */
    
    s_Rdark=s_LMA*(8.5341-130.6*s_Nmass-567.0*s_Pmass-0.0137*s_LMA+11.1*s_Vcmaxm+187600*s_Nmass*s_Pmass)*0.001;   // From Table 6 in Atkin et al 2015 New phytologist v.2.0 */
    //s_Rdark=(82.36*(s_LMA*1e-3)-0.1561)*(s_LMA*1e-3);  // (from Domingues et al 2007)
    //s_Rdark=0.01*s_Vcmax;          /* parameterization of Rdark commonly used in vegetation models */
    //s_Rdark=0.02*s_Vcmax-0.01;     /* parameterization of von Caemmerer 2000 Table 2.3 page 45 */
    
    s_Gamma = 38/Cair;      /* s_Gamma computed according to von Caemmerer 2000 formula: gamma=Kc*O*0.25/(2*Ko), with Kc=260 microbar, Ko=179mbar and O=210 mbar (the last value is from Fraquhar et al 1980, the first two one are from von Caemmerer 2000 table 2.3 page 45). gamma is set to 36.9 on Atkin et al 2015. Could be a global variable. v.2.0 */
    
    s_LCP = s_Rdark/phi;    /* Computation of the light compensation point from dark respiration and the quantum yield phi. By definition, Rdark is in micromolC/m^2/s and it is used in the Species::NPP() routine */
    
    cerr <<"s_Rdark for spp\t" <<s_name << "\t"<< s_Rdark<<"\ts_LCP\t"<< s_LCP << "\ts_nbext\t"<< s_nbext << "\ts_Vcmax\t" << s_Vcmax << "\ts_leaflifespan\t" << s_leaflifespan << "\ts_time_young\t" << s_time_young << "\ts_time_mature\t" << s_time_mature << "\ts_time_old\t" << s_time_old << endl;
    
    //s_factord13Ca = 1.0/(log(Wtot)-log(s_LCP)); /* deprecated v.2.0 */
    //s_factord13Cb = log(s_LCP); /* deprecated v.2.2 */
    
    
    /*** s_Seed field initialization ***/
    if (NULL==(s_Seed = new unsigned char[sites])) cerr<<"!!! Mem_Alloc\n";
    for(site=0;site<sites;site++) s_Seed[site]=0;
    
#ifdef MPI
    for(int i=0;i<4;i++) {
        if (NULL==(s_Gc[i] = new unsigned char[sites])) cerr<<"!!! Mem_Alloc\n";
        for(site=0;site<sites;site++) s_Gc[i][site]=0;
    }
#endif
}


/*############################
 ###     Species Seeds     ###
 ###   Species::FillSeed   ###
 ###  Species::UpdateSeed  ###
 ###   Species::AddSeed    ###
 #############################*/

/*###  Species::FillSeed  ###*/
/* creates one seed at point (col,row) */


void Species::FillSeed(int col, int row) {
    
    int site;
    
    if(col < cols) {
        
        if((row >=0) && (row < rows)) {
            site=col+cols*row;
#ifdef SEEDTRADEOFF
            s_Seed[site]++;                     /* ifdef SEETRADEOFF, s_Seed[site] is the number of seeds of this species at that site */
#else
            if(s_Seed[site]!=1) s_Seed[site]=1;  /* If s_Seed[site] = 0, site is not occupied, if s_Seed[site] > 1,
                                                  s_Seed[site] is the age of the youngest seed  */
#endif
        }
        
#ifdef MPI                                       /* on each processor a stripe of forest is simulated.
Nearest neighboring stripes are shared. Rque, this version is not valid ifdef SEEDTRADEOFF */
        else if((row+rows >=0) && (row < 0)) {
            site=col+cols*(row+rows);
            if(s_Gc[0][site]!=1) s_Gc[0][site]=1;
        }
        else if((row >=rows) && (row < 2*rows)) {
            site=col+cols*(row-rows);
            if(s_Gc[1][site]!=1) s_Gc[1][site]=1;
        }
#endif
    }
}


/*### Species::UpdateSeed ###*/
/* updates s_Seed field at a site*/

void Species::UpdateSeed(int age, int site) {
    
# ifdef MPI
    s_Gc[0][site]=s_Gc[1][site]=s_Gc[2][site]=s_Gc[3][site]=0;
#endif
#ifdef SEEDTRADEOFF
    s_Seed[site]=0;                              /* rk: with monthly timestep, should probably be modified, since as implemented now seeds are erased every month, same for freqgerm... a tree produces and disperses seeds every month : to be discussed */
#else
    if(age>0) s_Seed[site]=0;                   /* If site is occupied by a tree, the seed bank is emptied */
    else {                                      /* else, seed bank ages or disappears */
        if(s_Seed[site]==s_dormDuration) s_Seed[site]=0;
        else s_Seed[site]++;
    }
#endif
}


#ifdef MPI
/*########################################
 ###  Calculation of shared fields s_Gc ###
 ########################################*/
void Species::AddSeed() {
    /* Stripes shared by several processors are redefined */
    for(int site=0;site<sites;site++) {
        if(p_rank){
            if(!s_Seed[site]) s_Seed[site] = s_Gc[2][site];
            if(s_Seed[site]>1)
                if(s_Gc[2][site]) s_Seed[site] = min(s_Seed[site],s_Gc[2][site]);
        }
        if(p_rank<size-1){
            if(!s_Seed[site]) s_Seed[site] = s_Gc[3][site];
            if(s_Seed[site]>1)
                if(s_Gc[3][site]) s_Seed[site] = min(s_Seed[site],s_Gc[3][site]);
        }
    }
}
#endif


/*############################
 ###   Species::DeathRate  ###
 #############################*/

/* Here we assume a light-dependent version of the mortality.
 basal is the minimal species death rate, depending on the species wood density.
 When PPFD is smaller than the light compensation point, mortality risk is increased.
 When NDD is defined, this death rate is increased by density-dependence effect that impact survival of trees with a dbh<10cm . */


/* v.2.2 Simplify function Species::DeathRate -- JChave */
#ifdef NDD

/* Changed v.2.2*/
inline float Species::DeathRate(float PPFD, float dbh, float ndd) {
    
    float dr=0;
    //float basal=0.001+m*(1-s_wsg/0.85);
    //float densdepend;
    //float lightlimit = 1.0-PPFD/s_LCP;
    float add=2/(0.01-0.3*s_dmax);
    float bdd=-add*0.3*s_dmax;
    float dd=ndd*(add*dbh+bdd);        //the death rate is incremented by a negative density dependance term, which effect decreases whith the focal tree size (cf. Zhu et al 2015 Journal of Ecology) (size factor so that equals 0 at dbh=30cm and equals 2.5 at dbh=1 cm)
    //float senescent=0.01*dbh/s_dmax;
    
    /*if (nblivetrees > 0 && dbh<0.1*LH)
     densdepend=float(s_nbind)/float(nblivetrees);
     if (nblivetrees > 0)
     densdepend=float(s_nbind)/float(nblivetrees)*0.05*LH/dbh;                       // other try to make the density-dependance effect vary with tree size. 20-04-2015
     else densdepend=0;*/
    
    //dr=0.01*(1-s_wsg);                  // 5/08/2015: try with a basal rate function of wood density
    
    dr=basal;
    if (dbh<0.3*s_dmax) dr += dd;
    //if(lightlimit>0) dr += lightlimit;
    //if (dbh>s_dmax) dr += senescent;
    
    /* v.2.2  */
    if (iter == int(nbiter-1))  {
        sor[149]<< s_wsg << "\t" << basal << "\t"  << ndd << "\t" << dd << "\t" << dbh << "\t"  << dr   <<  "\n";
    }
    
    return dr*timestep;
}

#else

inline float Species::DeathRate(float PPFD, float dbh, int nppneg) {
    
    float dr=0;
    float basal=m-m1*s_wsg;
    //float lightlimit = 1.0-PPFD/s_LCP;
    
    //float juvenile=-7.5*dbh+0.75;
    //float juvenile=-1.67*dbh+0.167;
    
    dr=basal;
    //if(lightlimit>0) dr += lightlimit/timestep;
    // if (dbh<0.1) dr +=juvenile;
    
    if (nppneg > s_leaflifespan) {
        dr+=1.0/timestep;
    }
    
    if (iter == int(nbiter-1))  {
        sor[149]<< s_wsg << "\t" << basal << "\t"  << dbh << "\t"  << dr   <<  "\n";
    }
    
    return dr*timestep;
}

#endif


/*#############################################
 ###   Farquhar von Caemmerer Berry model  ###
 ###           Species:: NPP               ###
 #############################################*/

/* This function returns the leaf-level carbon assimilation rate  in micromoles C/m^2/s according to Farquhar-von Caemmerer-Berry model -- version 2.2 */

inline float Species::GPPleaf(float PPFD, float D, float T) {
    
    float A=0.0;
    
    /* ### Temperature dependencies ### */
    
    /* explicit real-time T dependence v.2.2 */
    //float theta=0.76+0.018*T-0.00037*T*T;         // this is the temperature dependency of theta provided by Bernacchi et al 2003 PCE
    float theta=0.7;                                // this is the fixed value of theta used by von Caemmerer 2000
    float denTB=1/(0.00831*(T+273.15));             // this is the denominator of the exponential term used for temperature dependencies provided by Bernacchi et al 2003 PCE
    //float KcT=exp(38.05-79.43*denTB);
    //float iKoT=exp(-20.30+36.38*denTB);
    float vCfactor=(25-T)/(298*8.314*(273+T));      // this is the factor of the exponential term used for temperature dependencies provided by von Caemmerer 2000
    float KcT=404*exp(vCfactor*59.36);
    float iKoT=1/248*exp(-vCfactor*35.94);
    float KmT=KcT*(1+210*iKoT)/Cair;
    //float GammaT=exp(19.02-37.83*denTB)/Cair;
    //float GammaT=KcT*iKoT*210*0.125/Cair;
    float GammaT=37*exp(vCfactor*23.4)/Cair;
    float VcmaxT=s_Vcmax*exp(26.35-65.33*denTB);
    //float VcmaxT=s_Vcmax*exp(vCfactor*58.52);
    float JmaxT=s_Jmax*exp(17.57-43.54*denTB);
    //float JmaxT=s_Jmax*exp(vCfactor*37);
    tempRday+=s_Rdark*exp((T-25)*0.1*log(3.09-0.0215*(25+T)));
    
    /* ### ci/ca ratio computation using stomatal conductance model from Medlyn et al 2011 ### */
    
    /* s_fci is the ratio of intercellular CO2 partial pressure divided by ambient partial pressure */
    //s_fci = minf(-0.04*s_d13C-0.3*(log(PPFD)-s_factord13Cb)*s_factord13Ca-0.57, 1.0);               //this was the initial version of ci/ca computation from d13C (see cernusak et al 2013) without explicit model of stomatal conductance; min added in order to prevent ci:ca bigger than 1 (even though Ehleringer et al 1986 reported some values above 1 (Fig3)
    s_fci=s_g1/(s_g1+sqrt(D));                                      /* this is the stomatal conductance model developped by Medlyn et al . (2011), see also Prentice et al 2014 Ecology Letters and Lin et al 2015 Nature Climate Change -- v.2.2 */
    
    if (iter == int(nbiter/2))  {
        sor[14]<< s_name << "\t" << PPFD << "\t" << s_LCP << "\t" << s_fci << "\n";
    }
    
    /* ### FvCB model ### */
    
    float I=alpha*PPFD;
    float J = (JmaxT+I-sqrt((JmaxT+I)*(JmaxT+I)-4*theta*JmaxT*I))*0.5/theta;
    A = minf(VcmaxT/(s_fci+KmT),0.25*J/(s_fci+2.0*GammaT))*(s_fci-GammaT);
    
    return A;
}


#ifdef DAILYLIGHT

/* dailyGPPleaf returns the assimilation rate (computed from Species::GPPleaf) averaged across the daily fluctuations in climatic conditions (light, VPD and T), in micromoles C/m^2/s */

/* new version v.2.2 */
inline float Species::dailyGPPleaf(float PPFD, float D, float T) {
    float dailyA=0;
    for (int i=0; i<24; i++) {
        dailyA+=GPPleaf(PPFD*daily_light[i], D*daily_vpd[i], T*daily_T[i]);    //daily_light is the averaged (across one year, meteo station Nouragues DZ) and normalized (from 0 to 1) daily fluctuation of light, with half-hour time step, during the day time (from 7am to 7pm, ie 12 hours in total), same for daily_vpd and daily_T. Taking into account these daily variation is necessary considering the non-linearity of FvCB model
    }
    dailyA*=0.0417;         // 0.0417=1/24 (24=12*2 = number of half hours in the 12 hours of daily light)
    tempRday*=0.0417;
    return dailyA;
}

#endif


/* Deprecated v.2.2. NPPleaf returns the net assimilation rate of carbon, ie the gross assimilation minus the respiration rate */
/*inline float Species::NPPleaf(float PPFD) {
 
 float NPP=0.0;
 
 #ifdef DAILYLIGHT
 NPP= dailyGPPleaf(PPFD);
 #else
 NPP=GPPleaf(PPFD);
 #endif
 
 NPP-=2*s_Rdark;   // Rdark is the leaf respiration.
 23/03/15: try to multiply Rdark by 2, in order to take stem and root repsiartion into account. Following Malhi 2012 Journal of Ecology, leaf dark respiration ~ 0.5 total autotrophic respiration. this multiplication by 2 is clearly discutable
 
 return NPP;
 }*/


/*############################################
 ############################################
 ############     Tree  class    ############
 ############################################
 ############################################*/

class Tree {
    
private:
#ifdef TREEFALL
    t_C;                    /* flexural force intensity */
#endif
    
public:
    int   t_site,           /* location */
    t_NPPneg;               /* diagnostic variable: number of consecutive timesteps with NPP<0 -- V.2.2 */
    float t_dbhmaxim,       /* dbh max */
#ifdef TREEFALL
    t_angle,                /* orientation of applied force */
#endif
    t_hmature,              /* reproductive size threshold */
    t_dbh,                  /* diameter at breast height (in m, but used in number of horizontal cells throughout all the code) */
    t_Tree_Height,          /* total tree height (in m, but used in number of vertical cells throughout all the code) */
    t_Crown_Depth,          /* crown depth (in m, but used in number of vertical cells throughout all the code) */
    t_Crown_Radius,         /* crown radius (in m, but used in number of horizontal cells throughout all the code)*/
#ifdef BASICTREEFALL
    t_Ct,                   /* flexural force threshold */
#endif
    t_GPP,                  /* tree gross primary productivity */
    t_NPP,                  /* tree net primary productivity (in gC/timestep) */
    t_Rday,                 /* leaf respiration during day */
    t_Rnight,               /* leaf respiration during night */
    t_Rstem,                /* stem respiration */
    /* t_NPPleaf,             Deprecated v.2.2.  tree net primary productivity at the leaf level (in gC/timestep) -- this variable is A (micromoleC/m^2/s), the assimilation rate of Farquhar (t_GPP minus the respiration term, see function Species::GPPleaf and Species NPPLeaf),
     times the surface area of leaf (in m2) which is dens*Pi*t_Crown_Radius^2 [*LH*LH ??]
     times 378.72*timestep (conversion from micromolC/s to gC/timestep) */
    t_PPFD,                 /* light intensity received by the tree (computed by Tree::Flux, depending of the LAI at the tree height) */
    t_D,                    /* VPD at tree height -- v.2.2 */
    t_T,                    /* Temperature at tree height -- v.2.2 */
    t_ddbh,                 /* increment of dbh per timestep */
    t_age,                  /* tree age */
    t_youngLA,              /* total young leaf area, in m2 -- v.2.2  */
    t_matureLA,             /* total mature leaf area, in m2 -- v.2.2  */
    t_oldLA,                /* total old leaf area, in m2 -- v.2.2  */
    t_leafarea,             /* total crown leaf area in m2 -- v.2.2  */
    t_dens,                 /* tree crown average leaf density in m2/m2 -- v.2.2  */
    t_litter;               /* tree litterfall at each timestep, in g (of dry mass) -- v.2.2  */
#ifdef NDD
    float *t_NDDfield;
#endif
    
    Species *t_s;
    
    unsigned short
    t_sp_lab,               /* species label */
    t_hurt;                 /* treefall index */
    
    Tree(){					/* constructor */
        t_sp_lab = 0;
        t_age = 0;
        t_hurt = 0;
        t_NPP=t_GPP=t_Rday=t_Rnight=t_Rstem=t_PPFD=t_D=t_T=0.0; /* new v.2.2 */
#ifdef TREEFALL
        t_C  = 0;
        t_angle = 0.0;
#endif
#ifdef BASICTREEFALL
        t_Ct = 0.0;
#endif
        t_dbh = t_Tree_Height = t_Crown_Radius = 0.0;
        
    };
    
    virtual ~Tree() {
#ifdef NDD
        delete [] t_NDDfield;
#endif
    };	/* destructor */
    
    
    void Birth(Species*,int,int);	/* tree birth */
    void Death();                   /* tree death */
    void Growth();                  /* tree growth */
    //float Flux();                 /* mean incident light flux (PPFD) reaching the top on the tree crown*/
    float Fluxh(int h);             /* compute mean light flux received by the tree crown layer at height h */
#ifdef TREEFALL
    int   Couple();                 /* force exerted by other trees */
#endif
    void DisperseSeed();            /* update Seed field */
#ifdef BASICTREEFALL
    void FallTree();                /* tree falling routine */
#endif
    void Update();                  /* tree evolution at each timestep */
    void Average();                 /* local computation of the averages */
    void CalcLAI();                 /* computation of LAI */
    //void TrunkLAI();                /* computation  of trunk LAI -- deprecated v.2.2 */
    void histdbh();                 /* computation of dbh histograms */
    
};


/*##############################################
 ####	            Tree birth              ####
 ####  called by BirthInit and UpdateTree   ####
 ##############################################*/

void Tree::Birth(Species *S, int nume, int site0) {
    
    t_site = site0;
    t_sp_lab = nume;            /* t_sp_lab is the species label of a site. Can be defined even if the site is empty
                                 (cf. persistence function defined in Chave, Am Nat. 2001) */
    t_NPPneg=0; /* v.2.2 */
    t_s = S+t_sp_lab;
    t_age = 1; /* in V.2.2 Isa replaced by  t_age = timestep; don't get it */
    t_hurt = 0;
    t_dbh = DBH0;
    t_ddbh=0.0;
    t_dbhmaxim = ((t_s->s_dmax)-t_dbh)*flor(1.0+log(genrand2())*0.01)+t_dbh;
    t_Tree_Height = H0;
    t_Crown_Radius  = ra0;
    t_Crown_Depth = de0;
    /* leaf ageing module -- new in v.2.2. */
    t_dens=dens;
    t_youngLA=t_dens*PI*t_Crown_Radius*LH*t_Crown_Radius*LH*t_Crown_Depth;
    /* initially, all stems have only young leaves -- LA stands for leaf area */
    t_matureLA=0; /* this is the amount of leaf area at maturity */
    t_oldLA=0; /* leaf area of senescing leaves */
    t_leafarea=t_youngLA;
    /* should be sum of LA young+mature+old, but the equation is correct initially */
    t_PPFD = Fluxh(int(t_Tree_Height));
    tempRday=0.0;
    
#ifdef DAILYLIGHT
    t_GPP=t_s->dailyGPPleaf(t_PPFD, t_D, t_T)*0.5*t_youngLA*189.3*timestep;
    /* The function Species::dailyGPPleaf returns the primary productivity at the leaf level for the plant in micromoles C/m^2/s.
     It is converted into gC per m^2 of leaf per timestep by "*189.3*timestep"
     accounting only for the ligt hours (12 hours instead of 24): 189.3=12*3600*365.25*12/1000000
     12 is the molar mass of C, and also the number of light hours in a day) */
#else
    t_GPP=t_s->GPPleaf(t_PPFD, t_D, t_T)*(0.5*t_youngLA)*189.3*timestep;
#endif
    
    /* Compute respiration terms -- new v.2.2 */ t_Rnight=t_s->s_Rdark*exp((tnight-25)*0.1*log(3.09-0.0215*(25+tnight)))*0.5*t_youngLA*189.3*timestep;                   // exp((tnight-25)*0.1*log(3.09-0.0215*(25+tnight))) is the temperature dependencies used by Atkin 2015 (equ1)
    t_Rday=tempRday*0.5*t_youngLA*189.3*timestep*0.40;
    tempRday=0.0;
    t_Rstem=39.6*0.5*t_dbh*PI*(t_dbh-0.5*t_dbh)*(t_Tree_Height-t_Crown_Depth)*exp((temp-25)/10*log(2))*378.7*timestep;          //this is the computation of Rstem for dbh<0.2 (see Tree::Growth); exp((temp-25)/10*log(2)) is the temperature dependency of Rstem, supposing a constant Q10=2, according to Ryan et al 1994 and Meir & Grace 2002
    t_NPP = 0.75*(t_GPP- 1.5*(t_Rday+t_Rnight)- 1.50*t_Rstem);                // 0.75 is to account for growth respiration (Thornley & canadell 2000,
    
    float hrealmax=t_s->s_hmax*t_dbhmaxim/(t_dbhmaxim+t_s->s_ah);
    t_hmature= (-11.47+0.90*hrealmax*LV)*NV*flor(1.0+log(genrand2())*0.01);  // this expression (hmature function of hmax, linear regression) is drawn from Wright et al 2005 JTE (from data on 22 species; 11 of which from BCI, Panama from Wright et al 2005; 11 of which from Thomas 1996 Oikos, Pasoh, Malaysia; all of which with hmax>20m)
#ifdef BASICTREEFALL
    //t_Ct = (t_s->s_hmax)*flor(1.0-vC*sqrt(-log(genrand2())));
    t_Ct = hrealmax*flor(1.0-vC*sqrt(-log(genrand2())));
    //t_Ct =t_s->s_hmax*t_s->s_dmax/(t_s->s_dmax+t_s->s_ah)*2*sqrt(-log(genrand2())*iPi);
    //t_Ct =flor(t_s->s_hmax*t_s->s_dmax/(t_s->s_dmax+t_s->s_ah)*(1 + 0.7*sqrt(-2*log(genrand2()))*cos(2*PI*genrand2())));
    //t_Ct = int(Cm*(t_s->s_hmax)*flor(1.0-sqrt(-log(genrand2()))*
    //		     vC*float(2*(genrand2i()%2)-1)));
#endif
    (t_s->s_nbind)++;
    nblivetrees++;
    
}


/*################################################
 #### Contribution of each tree in LAI field  ####
 ####          called by UpdateField          ####
 #################################################*/

void Tree::CalcLAI() {
    
    int haut0,h;
    if(t_age>0) {  // In V.2.2 -- but if(t_age) and if(t_age > 0) is equivalent
        dm=int(t_Crown_Radius);
        row0=t_site/cols;
        col0=t_site%cols;
        haut0 = int(t_Tree_Height-t_Crown_Depth);
        haut = int(t_Tree_Height)-1;                    // v.2.0 without the '-1' (ie haut=tree_height and not tree_height - 1), the tree contributes to LAI t_Crown_Depth+1 times, and contributes to shade above him. With this computation of each tree contribution in LAI, Fluxh can be compute with LAI at tree height (and not tree-height+1) without self shading.
        for(col=max(0,col0-dm);col<=min(cols-1,col0+dm);col++) {
            for(row=row0-dm;row<=row0+dm;row++) {                   /* loop over the tree crown */ /* no constraint on rows */
                xx=col0-col;
                yy=row0-row;
                if(xx*xx+yy*yy<=dm*dm)             /* is the voxel within crown? */
                    for(h=haut0;h <= haut;h++)   { /* loop over the crown depth */
                        /* Below if_else new in v.2.2. */
                        if (dm*(haut-haut0+1) <= 0) {
                            LAI3D[h][col+cols*row+SBORD] += t_dens*t_Crown_Depth;
                        }
                        else {
                            LAI3D[h][col+cols*row+SBORD] += t_dens;
                        }
                    }
            }
        }
    }
}


/* Effect of shading by the tree stems deprecated v.2.1 */

/*void Tree::TrunkLAI() {
 
 if(t_age)
 for(int h=0;h<int(t_Tree_Height);h++)
 LAI3D[h][t_site+SBORD]+=700;
 
 }*/


/*###################################################
 ####  Computation of PPFD right above the tree  ####
 ####    called by Tree::Birth and Growth        ####
 ####################################################*/

/* deprecated in v.2.1 */
/* Tree::Flux() computes the average light flux recieved by the tree crown top layer (ie at h=Tree_height) */

/*
 
 float Tree::Flux() {
 
 int count=0;
 float absorb,flux = 0.0;
 dm = int(t_Crown_Radius);
 if(dm == 0) {
 count=1;
 if (int(t_Tree_Height) >= HEIGHT) {                                                 // modification IM 01/2015 to avoid that trees (and especially top canopy trees) shade themselves
 absorb=0;
 }
 else {
 //absorb = LAI3D[int(t_Tree_Height+1)][t_site+SBORD];
 absorb = LAI3D[int(t_Tree_Height)][col+cols*row+SBORD];                //modified 22_04_2015, according to modification of Tree::CalcLAI
 }
 flux = exp(-absorb*klight);
 }
 else {
 row0=t_site/cols;
 col0=t_site%cols;
 for(col=max(0,col0-dm);col<min(cols,col0+dm+1);col++) {
 for(row=row0-dm;row<=row0+dm;row++) {                                           //loop over the tree crown
 xx=col0-col;
 yy=row0-row;
 if(xx*xx+yy*yy <= dm*dm) {                                                  //is the voxel within crown?
 count++;
 if (int(t_Tree_Height) >= HEIGHT) {
 absorb=0;
 }
 else {
 //absorb = LAI3D[int(t_Tree_Height+1)][col+cols*row+SBORD];
 absorb = LAI3D[int(t_Tree_Height)][col+cols*row+SBORD];                //modified 22_04_2015, according to modification of Tree::CalcLAI
 }
 flux += exp(-absorb*klight);
 
 }
 }
 }
 }
 flux*=Wtot/float(count);                                                                //this is the mean photosynthetic flux reaching the tree, computeted according to Beer-Lambert's law
 return flux;
 }
 
 */

/* Tree::Fluxh() computes the average light flux received by a tree crown layer at height h , and also the average VPD and T it is thriving in (modified 1/02/2016)*/


float Tree::Fluxh(int h) {
    
    int count=0;
    t_D=0;
    t_T=0;
    float absorb,flux = 0.0;
    dm = int(t_Crown_Radius);
    if(dm == 0) {
        count=1;
        if (h >= HEIGHT) {                                                 // modification IM 01/2015 to avoid that trees (and especially top canopy trees) shade themselves
            absorb=0;
        }
        else {
            //absorb = LAI3D[h+1][t_site+SBORD];
            absorb = LAI3D[h][t_site+SBORD];                    // modified 22-04-2015 according to modification in Tree::CalcLAI
        }
        flux = exp(-absorb*klight);
        t_D=0.25+sqrt(maxf(0.0 , 0.08035714*(7-absorb)));       /* new in v.2.2. 0.08035714=(1-0.25)^2/7, this expressions results from fit of observations of relationships between VPD and height within dense canopy (HOBO data on COPAS tower, Shuttleworth et al 1985; Camargo & Kapos 1995 journal of Tropical Ecology) */
        t_T=tmax-0.4285714*(minf(7,absorb));                    /* new in v.2.2. 0.4285714=3/7, assuming deltaT between the top canopy and dense understorey is constant = 3°C, could be refined.*/
    }
    else {
        row0=t_site/cols;
        col0=t_site%cols;
        for(col=max(0,col0-dm);col<min(cols,col0+dm+1);col++) {
            for(row=row0-dm;row<=row0+dm;row++) {                                           //loop over the tree crown
                xx=col0-col;
                yy=row0-row;
                if(xx*xx+yy*yy <= dm*dm) {                                                  //is the voxel within crown?
                    count++;
                    if (h >= HEIGHT) {
                        absorb=0;
                    }
                    else {
                        //absorb = LAI3D[h+1][col+cols*row+SBORD];
                        absorb = LAI3D[h][col+cols*row+SBORD];        // modified 22-04-2015 according to modification in Tree::CalcLAI
                    }
                    flux += exp(-absorb*klight);
                    t_D+=0.25+sqrt(maxf(0.0 , 0.08035714*(7-absorb)));
                    t_T+=tmax-0.4285714*(minf(7,absorb));
                }
            }
        }
    }
    flux*=Wmax/float(count);                                                                //this is the mean photosynthetic flux reaching the tree crown layer (at height h), computed according to Beer-Lambert's law
    t_D*=VPDmax/float(count);
    t_T*=1/float(count);
    return flux;
}


/*############################################
 ####            Tree growth              ####
 ####         called by UpdateTree        ####
 #############################################*/

void Tree::Growth() {
    
    t_GPP=0.0;
    t_NPP=0.0;
    t_PPFD=0.0;
    t_Rday=0.0;
    t_Rnight=0.0;
    t_Rstem=0.0;
    t_D=0.0;
    t_T=0.0;
    tempRday=0.0;
    t_ddbh=0.0;
    t_litter=0.0;
    float flush=0.0;
    float Dmean=0.0;
    float Tmean=0.0;
    float PPFDmean=0.0;
    
    t_age+= timestep; /* new v.2.2: increments are not 1 yr, but the duration of the timestep (usually 1 or <1, i.e. 1/12 if monthly, 1/365 if daily */
    
#ifdef CROWNGRADIENT
    int haut0, hauth, h;
    float f;
#endif
    
#ifdef CROWNGRADIENT
    //computation of t_GPP from the sum of GPP of each tree crown layer:
    if (int(t_Crown_Depth)<=1) {
#endif
        
        t_PPFD=Fluxh(int(t_Tree_Height));
        Dmean+=t_D;
        Tmean+=t_T;
        PPFDmean+=t_PPFD;
#ifdef DAILYLIGHT /* new v.2.2 */
        t_GPP+=t_s->dailyGPPleaf(t_PPFD, t_D, t_T);
#else
        t_GPP+=t_s->GPPleaf(t_PPFD, t_D, t_T);
#endif
        t_Rday+=tempRday;
        tempRday=0;
        
#ifdef CROWNGRADIENT
    }
    
    else {
        
        haut0=int(t_Tree_Height-t_Crown_Depth)+1;                   //modified 22-04-2015, according to modification of Tree::CalcLAI and Flux and Fluxh.
        hauth=int(t_Tree_Height);
        t_PPFD=Fluxh(int(t_Tree_Height));
        for(h=haut0; h<=hauth; h++) {
            f=Fluxh(h);
            Dmean+=t_D;  /* new v.2.2 */
            Tmean+=t_T; /* new v.2.2 */
            PPFDmean+=f; /* new v.2.2 */
            
#ifdef DAILYLIGHT
            t_GPP+=t_s->dailyGPPleaf(f, t_D, t_T);     /* new v.2.2  t_D and t_T are modified by Fluxh */
#else
            t_GPP+=t_s->GPPleaf(f, t_D, t_T);
#endif
            t_Rday+=tempRday;
            tempRday=0;
        }
        
        float inb_layer=1/float(hauth-haut0+1);
        t_GPP*=inb_layer;
        t_Rday*=inb_layer;
        Dmean*=inb_layer;
        Tmean*=inb_layer;
        PPFDmean*=inb_layer;
        
    }
    
#endif
    
    /*if (t_PPFD < 100)
     cout << "PPFD: " << t_PPFD << "  GPP: "<< t_GPP << "\n";*/
    
    t_GPP*=(0.5*t_youngLA+t_matureLA+0.5*t_oldLA)*189.3*timestep;
    /* new v.2.2: assumes an efficiency of 0.5 for young and old leaves vs. 1 for mature leaves in computation of GPP */
    /* The function Species::NPP returns the net primary productivity at the leaf level for the plant in micromoles C/m^2/s.
     It is converted into gC per m^2 of leaf per timestep by "*189.3*timestep"
     accounting only for the ligt hours (12 hours instead of 24): 189.3=12*3600*365.25*12/1000000
     12 is the molar mass of C, and also the number of light hours in a day) */
    if (t_GPP<0) cerr << "Warning: negative t_GPP!!";
    
    /* new v.2.2. sapwood thickness */
    float sapthick=0.04;
    if (t_dbh < 0.08) sapthick=0.5*t_dbh;
    
    /* new v.2.2 stem respiration */
    t_Rstem=39.6*sapthick*PI*(t_dbh-sapthick)*(t_Tree_Height-t_Crown_Depth)*378.7*timestep;
    t_Rstem*=exp((temp-25)/10*log(2));                                                    /* exp((temp-25)/10*log(2)) is the temperature dependency of Rstem, supposing a constant Q10=2, according to Ryan et al 1994 and Meir & Grace 2002 */
    t_Rday*=0.40*(0.5*t_youngLA+t_matureLA+0.5*t_oldLA)*189.3*timestep;
    t_Rnight=t_s->s_Rdark*exp((tnight-25)*0.1*log(3.09-0.0215*(25+tnight)))*(0.5*t_youngLA+t_matureLA+0.5*t_oldLA)*189.3*timestep;
    /* exp((tnight-25)*0.1*log(3.09-0.0215*(25+tnight))) is the temperature dependencies used by Atkin 2015 (equ1) */
    t_NPP = 0.75*(t_GPP- 1.5*(t_Rday+t_Rnight) - 1.50*t_Rstem);
    /* Rleaf=Rday+Rnight is multiplied by 1.5 to also account for fine root respiration (cf as in Fyllas et al 2014 and Malhi 2012); Rstem is multliplied by 1.5 to account for coarse root respiration (according to the shoot root biomass ratio of 0.2 - Jérôme's paper in prep- and also to branch respiration (Meir & Grace 2002, Cavaleri 2006, Asao 2005). */
    
    
    if(t_NPP<0){
        t_NPPneg++;
        t_NPP=t_GPP- 1.5*(t_Rday+t_Rnight) - 1.50*t_Rstem;
    }
    
    else {
        t_NPPneg=0;
        /**** NPP allocation to wood and tree size increment *****/
        float volume=2.0*t_NPP/(t_s->s_wsg) * fallocwood * 1.0e-6;
        /* volume in m^3: the first factor of 2 is to convert C into biomass.
         the 1/s_wsg to convert biomass into volume.
         the 1e-6 term converts cm^3 into m^3 (the sole metric unit in the model).
         fallocwood is the fraction of biomass allocated to aboveground wood (stem + branches) growth
         For the time being, we shall assume that a fixed proportion of NPP is allocated into AGB production.
         Currently, 0.20=%biomasse allocated to stem increment could be a global variable, even though this % allocation could in fact vary with resouce variation/co-limitation*/
        
        /*if (t_dbh<0.1) {
         volume*=5.56*t_dbh+0.44;
         }*/
        if (t_dbh>t_dbhmaxim) volume*=flor(3.0-2.0*t_dbh/t_dbhmaxim);
        /* taking into account wood elements recycling (ex. fallen branches etc...) */
        
        t_ddbh = flor( volume* 4.0/( 3.0*PI*t_dbh*LH*t_Tree_Height*LV ) )* NH;
        /* With V=pi*r^2*h, increment of volume = dV = 2*pi*r*h*dr + pi*r^2*dh,
         With isometric growth assumption (ddbh/dbh=dh/h)and dbh=2*r:
         dV=3/4*pi*dbh*h*ddbh, ddbh in m   */
        
        t_dbh += t_ddbh;
        
        /* update of tree height */
        t_Tree_Height = (t_s->s_hmax) * t_dbh/(t_dbh + (t_s->s_ah));
        
        /* update of tree crown depth */
        //t_Crown_Depth = de0+de1*(t_Tree_Height-H0);
        //if (t_Tree_Height<10.0) t_Crown_Depth = 0.5+0.163*t_Tree_Height;
        //else t_Crown_Depth = -0.47+0.26*t_Tree_Height;
        if (t_Tree_Height<5.0) t_Crown_Depth = 0.133+0.168*t_Tree_Height;
        else t_Crown_Depth = -0.48+0.26*t_Tree_Height;
        //if (t_Tree_Height>7.0) t_Crown_Depth = 3.09-0.14*t_Tree_Height+0.01*t_Tree_Height*t_Tree_Height;              // allometry deduced from Piste Saint-Elie dataset
        //else t_Crown_Depth = 0.15+0.35*t_Tree_Height;
        //t_Crown_Depth=exp(-1.169+1.098*log(t_Tree_Height));                         // 29/04/15: try with allometry used in Fyllas et al 2014 (see SI, from Poorter et al 2006)
        
        t_Crown_Radius  = 0.80+10.47*t_dbh-3.33*t_dbh*t_dbh;                                               /* update of tree crown radius */
        // allometry deduced from Piste Saint-Elie dataset
        //t_Crown_Radius  = ra0+ra1*t_dbh;
        //t_Crown_Radius=sqrt(iPi*exp(-1.853+1.888*log(t_Tree_Height)));               // 29/04/15: try with allometry used in Fyllas et al 2014 (see SI, from Poorter et al 2006)
        
        /**** NPP allocation to leaves *****/                                       /* rk: in this current scheme of leaf demography and phenology in three leaf age classes: only the old leaves generate treefall, and the dynamic of leaves cycle is generated by the dynamic of NPP, with a total leaf biomass varying - as opposed to De Weirdt et al 2012 in ORCHIDEE, but as in Wu et al 2016 but importantly without prescribing litterfall- */
        
        flush=2.0*t_NPP*falloccanopy*0.68/(t_s->s_LMA);                           /* this is to convert the NPP allocated to leaves (falloccanopy is the fraction of biomass assumed to be alloacted to canopy (leaves+reproductive organs+twigs) at each timestep - Malhi et al 2011-, 68% of which is alloctaed to leaves - chave et al 2008, Chave et al 2010-), in new leaf area (2 is to convert carbon mass in biomass and LMA to convert leaf biomass into leaf area).*/
        /*if (t_dbh<0.1) {
         flush*=5.56*t_dbh+0.44;
         }*/
    }
    
    t_litter=t_oldLA/(t_s->s_time_old);
    float new_mature=t_youngLA/(t_s->s_time_young);
    float new_old=t_matureLA/(t_s->s_time_mature);
    t_youngLA+= flush - new_mature;
    t_matureLA+= new_mature - new_old;
    t_oldLA+= new_old - t_litter;
    t_litter*=t_s->s_LMA;
    t_leafarea=t_youngLA + t_matureLA + t_oldLA;
    float crownvolume=PI*t_Crown_Radius*LH*t_Crown_Radius*LH*t_Crown_Depth*LV;
    t_dens=t_leafarea/crownvolume;
    
    if(isnan(t_dens)|isnan(t_leafarea)|isnan(flush)) {
        cout<<"problem in growth: nan";
    }
    
    if (iter == 2)  {
        sor[151] <<  t_sp_lab << "\t" << t_Tree_Height << "\t" << t_dbh << "\t"  << t_ddbh << "\t" << flush << "\t" << t_litter << "\t" << t_age << "\t" << t_leafarea << "\t" << t_youngLA<< "\t" << t_matureLA << "\t" << t_oldLA << "\t" << crownvolume << "\t" << t_Crown_Radius << "\t" << t_Crown_Depth << "\t" << t_dens  <<  "\t" << t_PPFD  <<"\t" << t_GPP  <<"\t" << t_NPP <<"\t" << t_Rstem <<"\t" << t_Rday  <<"\t" << t_Rnight << "\t" << t_site << "\t" << Dmean << "\t" << Tmean << "\t" << LAI3D[int(t_Tree_Height)][t_site+SBORD] << "\t" << LAI3D[int(t_Tree_Height-t_Crown_Depth)+1][t_site+SBORD] << "\t" <<  PPFDmean <<  "\n";
    }
    if (iter == int(nbiter/2))  {
        sor[152] <<  t_sp_lab << "\t" << t_Tree_Height << "\t" << t_dbh << "\t"  << t_ddbh << "\t" << flush << "\t" << t_litter << "\t" << t_age << "\t" << t_leafarea << "\t" << t_youngLA<< "\t" << t_matureLA << "\t" << t_oldLA << "\t" << crownvolume << "\t" << t_Crown_Radius << "\t" << t_Crown_Depth << "\t" << t_dens  <<  "\t" << t_PPFD  <<"\t" << t_GPP  <<"\t" << t_NPP <<"\t" << t_Rstem <<"\t" << t_Rday  <<"\t" << t_Rnight << "\t" << t_site << "\t" << Dmean << "\t" << Tmean << "\t" << LAI3D[int(t_Tree_Height)][t_site+SBORD] << "\t" << LAI3D[int(t_Tree_Height-t_Crown_Depth)+1][t_site+SBORD]<< "\t"  << PPFDmean << "\n";      }
    if (iter == int(nbiter-1))  {
        sor[153] <<  t_sp_lab << "\t" << t_Tree_Height << "\t" << t_dbh << "\t"  << t_ddbh << "\t" << flush << "\t" << t_litter << "\t" << t_age << "\t" << t_leafarea << "\t" << t_youngLA<< "\t" << t_matureLA << "\t" << t_oldLA << "\t" << crownvolume << "\t" << t_Crown_Radius << "\t" << t_Crown_Depth << "\t" << t_dens  <<  "\t" << t_PPFD  <<"\t" << t_GPP  <<"\t" << t_NPP <<"\t" << t_Rstem <<"\t" << t_Rday  <<"\t" << t_Rnight  << "\t" << t_site <<"\t" << Dmean << "\t" << Tmean << "\t" << LAI3D[int(t_Tree_Height)][t_site+SBORD] << "\t" << LAI3D[int(t_Tree_Height-t_Crown_Depth)+1][t_site+SBORD] << "\t" << PPFDmean << "\n";       }
    
    
    
    if (t_site==2500)  {
        sor[11]<< iter << "\t" << t_site << "\t" << t_sp_lab << "\t" << t_Tree_Height << "\t" << t_dbh << "\t"  << t_ddbh << "\t" << flush << "\t" << t_litter << "\t" << t_age << "\t" << t_leafarea << "\t" << t_youngLA<< "\t" << t_matureLA << "\t" << t_oldLA << "\t" << t_Crown_Radius << "\t" << t_Crown_Depth << "\t" << t_dens  <<  "\t" << t_PPFD  <<"\t" << t_GPP  <<"\t" << t_NPP <<"\t" << t_Rstem <<"\t" << t_Rday  <<"\t" << t_Rnight  << "\t" << Dmean << "\t" << Tmean << "\t" << LAI3D[int(t_Tree_Height)][t_site+SBORD] << "\t" << LAI3D[int(t_Tree_Height-t_Crown_Depth)+1][t_site+SBORD] << "\t" << PPFDmean <<"\n";
    }
    
    if (t_site==10380)  {
        sor[12]<< iter << "\t" << t_site << "\t" << t_sp_lab << "\t" << t_Tree_Height << "\t" << t_dbh << "\t"  << t_ddbh << "\t" << flush << "\t" << t_litter << "\t" << t_age << "\t" << t_leafarea << "\t" << t_youngLA<< "\t" << t_matureLA << "\t" << t_oldLA << "\t" << t_Crown_Radius << "\t" << t_Crown_Depth << "\t" << t_dens  <<  "\t" << t_PPFD  <<"\t" << t_GPP  <<"\t" << t_NPP <<"\t" << t_Rstem <<"\t" << t_Rday  <<"\t" << t_Rnight  << "\t" << Dmean << "\t" << Tmean << "\t" << LAI3D[int(t_Tree_Height)][t_site+SBORD] << "\t" << LAI3D[int(t_Tree_Height-t_Crown_Depth)][t_site+SBORD] <<"\t" << PPFDmean <<"\n";
    }
    
    if (t_site==100950)  {
        sor[13]<< iter << "\t" << t_site << "\t" << t_sp_lab << "\t" << t_Tree_Height << "\t" << t_dbh << "\t"  << t_ddbh << "\t" << flush << "\t" << t_litter << "\t" << t_age << "\t" << t_leafarea << "\t" << t_youngLA<< "\t" << t_matureLA << "\t" << t_oldLA << "\t" << t_Crown_Radius << "\t" << t_Crown_Depth << "\t" << t_dens  <<  "\t" << t_PPFD  <<"\t" << t_GPP  <<"\t" << t_NPP <<"\t" << t_Rstem <<"\t" << t_Rday  <<"\t" << t_Rnight  << "\t" << Dmean << "\t" << Tmean <<"\t" << LAI3D[int(t_Tree_Height)][t_site+SBORD] << "\t" << LAI3D[int(t_Tree_Height-t_Crown_Depth)][t_site+SBORD] << "\t" << PPFDmean <<"\n";
    }
    
    if (t_site==12090)  {
        sor[155]<< iter << "\t" << t_site << "\t" << t_sp_lab << "\t" << t_Tree_Height << "\t" << t_dbh << "\t"  << t_ddbh << "\t" << flush << "\t" << t_litter << "\t" << t_age << "\t" << t_leafarea << "\t" << t_youngLA<< "\t" << t_matureLA << "\t" << t_oldLA << "\t" << t_Crown_Radius << "\t" << t_Crown_Depth << "\t" << t_dens  <<  "\t" << t_PPFD  <<"\t" << t_GPP  <<"\t" << t_NPP <<"\t" << t_Rstem <<"\t" << t_Rday  <<"\t" << t_Rnight  << "\t" << Dmean << "\t" << Tmean <<"\t" << LAI3D[int(t_Tree_Height)][t_site+SBORD] << "\t" << LAI3D[int(t_Tree_Height-t_Crown_Depth)][t_site+SBORD] << "\t" << PPFDmean <<"\n";
    }
    
    if (t_site==120090)  {
        sor[156]<< iter << "\t" << t_site << "\t" << t_sp_lab << "\t" << t_Tree_Height << "\t" << t_dbh << "\t"  << t_ddbh << "\t" << flush << "\t" << t_litter << "\t" << t_age << "\t" << t_leafarea << "\t" << t_youngLA<< "\t" << t_matureLA << "\t" << t_oldLA << "\t" << t_Crown_Radius << "\t" << t_Crown_Depth << "\t" << t_dens  <<  "\t" << t_PPFD  <<"\t" << t_GPP  <<"\t" << t_NPP <<"\t" << t_Rstem <<"\t" << t_Rday  <<"\t" << t_Rnight  << "\t" << Dmean << "\t" << Tmean <<"\t" << LAI3D[int(t_Tree_Height)][t_site+SBORD] << "\t" << LAI3D[int(t_Tree_Height-t_Crown_Depth)][t_site+SBORD] << "\t" << PPFDmean <<"\n";
    }
    
    if (t_site==150667)  {
        sor[157]<< iter << "\t" << t_site << "\t" << t_sp_lab << "\t" << t_Tree_Height << "\t" << t_dbh << "\t"  << t_ddbh << "\t" << flush << "\t" << t_litter << "\t" << t_age << "\t" << t_leafarea << "\t" << t_youngLA<< "\t" << t_matureLA << "\t" << t_oldLA << "\t" << t_Crown_Radius << "\t" << t_Crown_Depth << "\t" << t_dens  <<  "\t" << t_PPFD  <<"\t" << t_GPP  <<"\t" << t_NPP <<"\t" << t_Rstem <<"\t" << t_Rday  <<"\t" << t_Rnight  << "\t" << Dmean << "\t" << Tmean <<"\t" << LAI3D[int(t_Tree_Height)][t_site+SBORD] << "\t" << LAI3D[int(t_Tree_Height-t_Crown_Depth)][t_site+SBORD] << "\t" << PPFDmean <<"\n";
    }
    
    /*if (t_leafarea==0.0) {
     cout << "leafarea=0 !" << "\t" << t_site << "\t" << t_PPFD  << "\t" << t_s->s_LCP  << "\t" << t_NPPneg  << "\t" << t_s->s_name  << "\t" << t_Tree_Height <<  "\n";
     }*/
    
}


/*####################################################
 ####           Death of the tree                ####
 ####         called by Tree::Update             ####
 ####################################################*/

void Tree::Death() {
    
    t_age=0;
    t_dbh = t_Tree_Height = t_Crown_Radius = t_Crown_Depth= 0.0;
    t_hurt = 0;
#ifdef TREEFALL
    t_angle = 0.;
    t_C = t_Ct = 0;
#endif
    if ((t_s->s_nbind)>0) (t_s->s_nbind)--;
    nblivetrees--;
    t_s = NULL;
    
}


/*#################################
 ####      Seed dispersal      ####
 ####  called by UpdateField   ####
 #################################*/

void Tree::DisperseSeed(){
    //if((t_Tree_Height>=t_hmature)&&(t_PPFD>2.0*(t_s->s_LCP))&&(fmod(t_age,(t_s->s_freqgerm))==0)) {         /* Test maturity germination frequency */
    /* New v.2.0 reproduction can only occur for trees that receive enough
     light (twice the LCP) */
    /* New v.2.1 threshold of maturity is defined as a size threshold
     (and not age as before), following Wright et al 2005 JTE */
    if((t_Tree_Height>=t_hmature)&&(t_PPFD>2.0*(t_s->s_LCP))) {
        
        float rho,theta;
        int nbs=0;
#ifdef SEEDTRADEOFF
        nbs=int(t_NPP*2*falloccanopy*0.08*0.5*(t_s->s_iseedmass));    /* nbs is the number of seeds produced at this time step; it is computed from the NPP (in g) allocated to reproductive organs -fruits and seeds-, *2 is to convert in biomass,  * 0.40 is to obtained the NPP allocated to canopy (often measured as litterfall), drawn from Malhi et al 2011 Phil. trans roy. Soc. and 0.08 is the part of litterfall corresponding the fruits+seeds, drawn from Chave et al 2008 JTE; assumed to be twice the biomass dedicated to seeds only (ie without fruits), and then divided by the mass of a seed to obtain the number of seeds */
        nbs*=0.2;
        //nbs=(int)nbs;
#else
        nbs=nbs0;
#endif
        
        for(int ii=0;ii<nbs;ii++){                                                 /* Loop over number of produced seeds */
            
            rho = 2*((t_s->s_ds)+t_Crown_Radius)*float(sqrt(fabs(log(genrand2()*iPi))));    /* Dispersal distance rho: P(rho) = rho*exp(-rho^2) */
            theta = float(deuPi*genrand2());                                                /* Dispersal angle theta */
            t_s->FillSeed(flor(int(rho*cos(theta))+t_site%cols), /* column */               /* Update du Field s_Seed */
                          int(rho*sin(theta))+t_site/cols);      /* line */
        }
    }
    
}


/*##################################
 ####   Tree death and growth   ####
 ####   called by UpdateTree    ####
 ##################################*/

void Tree::Update() {
    
    if(t_age) {
        
#ifdef NDD
        if(int(genrand2() + t_s->DeathRate(t_PPFD, t_dbh, t_NDDfield[t_sp_lab]))){
#else
            if(int(genrand2() + t_s->DeathRate(t_PPFD, t_dbh, t_NPPneg))){
                /* Natural death caused by unsustained photosynthesis and density dependance */
#endif
                nbm_n1++;
                /* Increments number of deaths */
                if(t_dbh*LH>0.1) nbm_n10++;
                /* same but only for trees > 10cm */
                /* New v.2.2. new outputs */
                if(iter == 2) sor[146] << "N\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
                if(iter == int(nbiter/2)) sor[147]<< "N\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
                if(iter == int(nbiter-1)) sor[148]<< "N\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
                Death();
            }
#ifdef BASICTREEFALL
            else if(t_Tree_Height<2*t_hurt*genrand2()) {
                /* Death caused by a treefall (t_hurt) */
                nbm_c1++;   /* Increments number of treefalls */
                if(t_dbh*LH>0.1) nbm_c10++;/* same but only for trees > 10cm */
                if(iter == 2) sor[146] << "OT\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
                if(iter == int(nbiter/2)) sor[147]<< "OT\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
                if(iter == int(nbiter-1)) sor[148]<< "OT\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
                Death();
#ifdef NDD
            }
#else
        }
#endif
        
#endif
        else {/* If no death, then growth */
            t_hurt = 0;
            Growth();
        }
    }
}


/*##################################
 ####           Treefall         ####
 #### called by UpdateTreefall   ####
 ####################################*/

#ifdef BASICTREEFALL
void Tree::FallTree() {
#ifdef TREEFALL           /* above a given stress threshold the tree falls */
    if(Couple()>t_Ct) {
#else
        if(genrand2()*t_Tree_Height > t_Ct){
            float t_angle = float(deuPi*genrand2());
#endif
            int h;
            float hvrai = t_Tree_Height*LV;
            nbTreefall1++;
            if(t_dbh*LH>0.1) nbTreefall10++;
            Thurt[0][t_site+sites] = int(t_Tree_Height);
            row0=t_site/cols;       /* fallen stem destructs other trees */
            col0=t_site%cols;
            haut = int(hvrai*NH);
            for(h=1;h<haut;h++) {   /* loop on the fallen stem (horizontally) */
                xx=int(flor(col0+h*cos(t_angle)));
                if(xx<cols){
                    yy=   int(row0+h*sin(t_angle));
                    Thurt[0][xx+(yy+rows)*cols] = int(t_Tree_Height);
                }
            }
            xx=col0+int((hvrai*NH-t_Crown_Radius)*cos(t_angle));
            /* fallen crown destructs other trees */
            yy=row0+int((hvrai*NH-t_Crown_Radius)*sin(t_angle));
            h = int(t_Crown_Radius);
            for(col=max(0,xx-h);col<min(cols,xx+h+1);col++) {
                /* loop on the fallen crown (horizontally) */
                for(row=yy-h;row<yy+h+1;row++) {
                    if((col-xx)*(col-xx)+(row-yy)*(row-yy)<h*h) Thurt[0][col+(row+rows)*cols] = int((t_Tree_Height-t_Crown_Radius*NV*LH)*0.5);
                }
            }
            if(iter == 2) sor[146] << "T\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
            if(iter == int(nbiter/2)) sor[147]<< "T\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
            if(iter == int(nbiter-1)) sor[148]<< "T\t" << t_sp_lab << "\t" << t_dbh << "\t" << t_age << "\t" << t_Tree_Height <<  "\n";
            Death();
#ifdef TREEFALL
        }
#else
    }
#endif
}

#ifdef TREEFALL
int Tree::Couple() {
    int site2,quadist,h,haut0;
    float fx, fy,temp,lai;
    dm = int(t_Crown_Radius);
    haut = int(t_Tree_Height);
    haut0 = int(t_Tree_Height-t_Crown_Depth);
    if(dm){
        row0=t_site/cols;
        col0=t_site%cols;
        fx = fy = 0.0;
        for(col=max(0,col0-dm);col<min(cols,col0+dm+1);col++) {
            for(row=row0-dm;row<=row0+dm;row++) {
                xx=col0-col;
                yy=row0-row;
                quadist = xx*xx+yy*yy;
                if((quadist<=dm*dm)&&quadist) {
                    //site2 = col+cols*(row+RMAX); //modif 23/03/2011
                    site2 = col+cols*row+SBORD;
                    for(h=haut0;h<=haut;h++) {
                        if(haut<HEIGHT) lai = LAI3D[haut][site2]-LAI3D[haut+1][site2];
                        else  lai = LAI3D[haut][site2];
                        if(lai>dens) { // needs to be changed when TREEFALL is revised
                            temp = 1.0/sqrt(float(quadist));
                            if(temp>0.0) {
                                fx += float(xx)*temp;
                                fy += float(yy)*temp;
                            }
                        }
                    }
                }
            }
        }
        t_C = int(sqrt(fx*fx+fy*fy)*t_Tree_Height);
        if(fx!=0.0) t_angle=atan2(fy,fx);
        else t_angle = Pis2*sgn(fy);
    }
    else{t_C = 0; t_angle = 0.0; }
    return t_C;
}
#endif
#endif


/*#####################################################
 ####      For Average and OutputField             ####
 ######################################################*/

void Tree::Average() {
    if(t_age>0) {
        if(t_dbh*LH >= 0.1) {
            (t_s->s_chsor[1])++;
            t_s->s_chsor[6] += t_dbh*LH*t_dbh*LH;
        }
        if(t_dbh*LH >= 0.3) (t_s->s_chsor[2])++;
        t_s->s_chsor[3] += t_dbh*LH*t_dbh*LH;
        t_s->s_chsor[4] += t_NPP*1.0e-6;
        t_s->s_chsor[5] += t_GPP*1.0e-6;
        t_s->s_chsor[7] += 0.0673*pow(t_s->s_wsg*t_Tree_Height*LV*t_dbh*t_dbh*LH*LH*10000, 0.976);  // this is the allometrtic equ 4 in Chave et al. 2014 Global Change Biology to compute above ground biomass
        t_s->s_chsor[8] += t_Rday*1.0e-6;
        t_s->s_chsor[9] += t_Rnight*1.0e-6;
        t_s->s_chsor[10] += t_Rstem*1.0e-6;
        t_s->s_chsor[11] += t_litter*1.0e-6;
    }
}
void Tree::histdbh() {
    if(t_age&&(t_age<500)) nbdbh[int(100.*t_dbh*LH)]++;
    // compute the diameter distribution density
    // where dbh is in cm (it is in number of horizontal cells throughout the code)
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

//#define NIV_PERS 20

/* Class objects */

Species *S=NULL;
Tree *T=NULL;

//int *distr;


/*############################################
 ############################################
 ############     MAIN PROGRAM    ###########
 ############################################
 ############################################*/

int main(int argc,char *argv[]) {
    
    /***********************/
    /*** Initializations ***/
    /***********************/
#ifdef MPI   /* Lookup processor number / total number of processors */
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&p_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
#else
    mpi_rank = 0;
    mpi_size = 1;
#endif
#ifdef easyMPI
    MPI::Init();
    easympi_rank = MPI::COMM_WORLD.Get_rank();
    easympi_size = MPI::COMM_WORLD.Get_size();
    //MPI_Get_library_version(version, &len);
    //std::cout << "Hello, world!  I am " << rank << " of " << size << "(" << version << ", " << len << ")" << std::endl;
    //cout<< "Hello, world!  I am " << easympi_rank << " of " << easympi_size<<endl;
#else
    easympi_rank = 0;
    easympi_size =1;
#endif
    
    for(int argn=1;argn<argc;argn++){ /* Arguments of the input and output files */
        if(*argv[argn] == '-'){
            switch(*(argv[argn]+1)){
                case 'i':
                    bufi = argv[argn]+2;
                    break;
                case 'o':
                    buf = argv[argn]+2;
                    break;
                case 'n':
                    easympi_rank=atoi(argv[argn]+2); /* new v.2.2. */
            }
        }
    }
    cout<< easympi_rank;
    int t=(int) time(NULL);
    int seed =3*t+2*(easympi_rank+1)+1;
    sgenrand2i(seed);
    sgenrand2(seed);
    cout<< "On proc #" << easympi_rank << " seed: " << seed <<  " rng: "<< genrand2() << endl;
    cout << genrand2()<< endl;
    cout << genrand2()<< endl;
    sprintf(nomfi,"%s",bufi);
    //in this version of the code, the code user is expected to use the input txt file name with its extension name, as eg. "input.txt", when setting the program's argument
    //sprintf(nomfi,"%s.txt",bufi);
    
#ifdef outputMPI
    //nothing
#else
    sprintf(nomf0,"%s_%i_par.txt",buf, easympi_rank);                     /* Files with general output info */
    cout<< "On proc #" << easympi_rank << " seed: " << seed <<  " rng: "<< genrand2() <<  endl;
    out.open(nomf0, ios::out);
    if(!out) cerr<< "ERROR with par file"<< endl;
    sprintf(nomf0,"%s_%i_info.txt",buf, easympi_rank);
    out2.open(nomf0, ios::out);
    if(!out2) cerr<< "ERROR with info file"<< endl;
#endif
    Initialise();      /* Read global parameters */
    AllocMem();        /* Memory allocation */
    BirthInit();       /* Initial configuration of the forest */
    out.close();
    
    /***********************/
    /*** Evolution loop  ***/
    /***********************/
    
    double start_time,stop_time, duration=0.0;           /* for simulation duration */
    stop_time = clock();
    for(iter=0;iter<nbiter;iter++) {
        start_time = stop_time;
        
        /* set the iteration environment -- nb: the current structure of code suppose that environment is periodic (a period = a year), if one want to make climate vary, with interannual variation and climate change along the simulation, one just need to provide the full climate input of the whole simulation (ie number of columns=iter and not iterperyear) and change iterperyear by nbiter here. */
        temp=Temperature[iter%iterperyear];
        tmax=DailyMaxTemperature[iter%iterperyear];
        tnight=NightTemperature[iter%iterperyear];
        precip=Rainfall[iter%iterperyear];
        WS=WindSpeed[iter%iterperyear];
        Wmax=MaxIrradiance[iter%iterperyear]*1.678;       // 1.678 is to convert irradiance from W/m2 to micromol of PAR /s /m2, the unit used in the FvCB model of photosynthesis
        Wmean=MeanIrradiance[iter%iterperyear];            // still in W/m2
        e_s=SaturatedVapourPressure[iter%iterperyear];
        e_a=VapourPressure[iter%iterperyear];
        VPD=VapourPressureDeficit[iter%iterperyear];
        VPDday=DailyVapourPressureDeficit[iter%iterperyear];
        VPDmax=DailyMaxVapourPressureDeficit[iter%iterperyear];
        
        Evolution();
        stop_time = clock();
        duration +=flor(stop_time-start_time);
        
        /****** Final pattern v.2.2. ******/
        if (iter==nbiter-2) {
            for(int row=0;row<rows;row++)
                for(int col=0;col<cols;col++){
                    sor[10] << col << "\t" << row << "\t" << T[col + cols*row].t_age << "\t" << T[col + cols*row].t_dbh << "\t" << T[col + cols*row].t_Tree_Height << "\t" << T[col + cols*row].t_Crown_Radius << "\t" << T[col + cols*row].t_Crown_Depth << "\t" << T[col + cols*row].t_sp_lab << endl;
                }
        }
        
        //        if(nbiter>100 && iter > (nbiter-101))
        //            for(int row=0;row<rows;row++)
        //                for(int col=0;col<cols;col++){
        //                    sor[17] << iter << "\t" << col+cols*row << "\t" << col << "\t" << row << "\t" << T[col + cols*row].t_age << "\t" << T[col + cols*row].t_sp_lab << "\t" << T[col + cols*row].t_dbh << "\t" << T[col + cols*row].t_Tree_Height << "\t" << T[col + cols*row].t_Crown_Radius << "\t" << T[col + cols*row].t_Crown_Depth << "\t" << T[col + cols*row].t_leafarea << "\t" << T[col + cols*row].t_dens << "\t" << T[col + cols*row].t_GPP << "\t" << T[col + cols*row].t_hurt << endl;
        //                }
    }
    
    for(int sp=1;sp<=numesp;sp++) sor[15] << S[sp].s_name << "\t" << S[sp].s_Nmass << "\t" << S[sp].s_Pmass << "\t" << S[sp].s_LMA << "\t" << S[sp].s_Vcmax << "\t" << S[sp].s_Jmax << "\t" << S[sp].s_Rdark << "\t" << S[sp].s_LCP << "\n";
    
    /*************************/
    /*** End of simulation ***/
    /*************************/
    
    float durf = duration/double(CLOCKS_PER_SEC);        /* output of the effective CPU time */
#ifdef MPI
    MPI_Reduce(&durf,&durf,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
#endif
    if(!mpi_rank) {
        cout << "\n";
#ifdef MPI
        out2 << "Number of processors : "<< mpi_size << "\n";
#endif
        out2 << "Average computation time : "<< durf/float(mpi_size) << " seconds.\n";
        out2 << "End of simulation.\n";
        cout << "\nNumber of processors : "<< mpi_size << "\n";
        cout << "Average computation time : "<< durf/float(mpi_size) << " seconds.\n";
        cout << "End of simulation.\n";
    }
    out2.close();
    
    /*Free dynamic memory */ /* added in oct2013 */
    FreeMem();
    
#ifdef easyMPI
    MPI::Finalize();
#endif
    
    exit(1);
}


/*###########################################
 ############################################
 #######  Initialization routines    ########
 ############################################
 ############################################*/

void Initialise() {
    
    int ligne;
    
    fstream In(nomfi, ios::in);
    
    /*** Initialization of the simulation parametres ***/
    /***************************************************/
    
    if (In) {
        for(ligne=0;ligne<4;ligne++) In.getline(buffer,128,'\n');
        
        /*General parameters */
        In >> cols; In.getline(buffer,128,'\n');
        In >> rows; In.getline(buffer,128,'\n');
        sites = rows*cols;
        In >> nbiter; In.getline(buffer,128,'\n');
        In >> iterperyear; In.getline(buffer,128,'\n');
        timestep=1/float(iterperyear);
        In >> NV; In.getline(buffer,128,'\n');
        In >> NH; In.getline(buffer,128,'\n');
        LV = 1.0/NV;
        LH = 1.0/NH;
        In >> nbout; In.getline(buffer,128,'\n');
        if(nbout) freqout = nbiter/nbout;
        In >> numesp; In.getline(buffer,128,'\n');
        In >> p; In.getline(buffer,128,'\n');
        In >> disturb_iter; In.getline(buffer,128,'\n');
        In >> disturb_intensity; In.getline(buffer,128,'\n');
        for (int i=0; i<=23; i++) In >> daily_light[i];
        In.getline(buffer,128,'\n');
        for (int i=0; i<=23; i++) In >> daily_vpd[i];
        In.getline(buffer,128,'\n');
        for (int i=0; i<=23; i++) In >> daily_T[i];
        In.getline(buffer,128,'\n');
        In.getline(buffer,128,'\n');
        
        /*Characters shared by species */
        In >> klight; In.getline(buffer,128,'\n');
        In >> phi; In.getline(buffer,128,'\n');
        In >> vC; In.getline(buffer,128,'\n');
        In >> DBH0; In.getline(buffer,128,'\n');
        In >> H0; In.getline(buffer,128,'\n');
        In >> ra0; In.getline(buffer,128,'\n');
        In >> ra1; In.getline(buffer,128,'\n');
        In >> de0; In.getline(buffer,128,'\n');
        In >> de1; In.getline(buffer,128,'\n');
        In >> dens; In.getline(buffer,128,'\n');
        In >> fallocwood; In.getline(buffer,128,'\n');
        In >> falloccanopy; In.getline(buffer,128,'\n');
        In >> Cseedrain; In.getline(buffer,128,'\n');
        In >> nbs0; In.getline(buffer,128,'\n');
        In >> m; In.getline(buffer,128,'\n');
        In >> m1; In.getline(buffer,128,'\n');
        In >> Cair; In.getline(buffer,128,'\n');
        
        /* tests for v.2.2 */
        //klight+=0.5*genrand2();
        //phi+=0.05*genrand2();
        //fallocwood+=0.30*genrand2();
        //falloccanopy+=0.30*genrand2();
        //m+=0.09*genrand2();
        //m1+=0.09*genrand2();
        //Cseedrain+=490*genrand2();
        
        DBH0 *= NH;
        H0 *= NV;
        ra0 *= NH;
        de0 *= NV;
        alpha=4*phi;
        /* apparent quantum yield to electron transport in mol e-/mol photons */
        // see Mercado et al 2009 , the conversion of the apparent quantum yield in micromolCO2/micromol quantum into micromol e-/micxromol quantum is done by multipliyng by 4, since four electrons are needed to regenerate RuBP.
        /* alpha is fixed at 0.3 mol e-/mol photons in Medlyn et al 2002
         but see equ8 and Appendix 1 in Farquahr et al 1980: it seems that alpha should vary with leaf thickness: there is a fraction of incident light which is lost by absorption by other leaf parts than the chloroplast lamellae, and this fraction f may increase with leaf thickness.
         With the values of the paper: alpha= 0.5*(1-f)=0.5*(1-0.23)=0.385, but this is a theoretical value and observations often report lower values (see ex discussion in medlyn et al 2005 Tree phsyiology, Lore Veeryckt values, Mercado et al 2009 Table 10, Domingues et al. 2014)*/
        
        /* deprecated section below v.2.2 */
#ifdef DAILYLIGHT
        //for (int i=0; i<=23; i++) cout<< "Daily light" << i << "\t"  << daily_light[i] << endl;
        //for (int i=0; i<=23; i++) cout<< "Daily VPD" << i << "\t"  << daily_vpd[i] << endl;
        //for (int i=0; i<=23; i++) cout<< "Daily T" << i << "\t"  << daily_T[i] << endl;
#endif
    }
    
    else {
        cout<< "ERROR with the input file of parameters" << endl;
    }
    
    
    /*** Information in file info ***/
    /***********************************/
    
    if(!mpi_rank){
        out2 << "\nTROLL simulator\n\n";
        out2 << "\n   2D discrete network: horizontal step = " << LH
        << " m, one tree per "<< LH*LH << " m^2 \n\n";
        out2 << "\n   Tree : (t_dbh,t_Tree_Height,t_Crown_Radius,t_Crown_Depth) \n\n";
        out2 << "\n            + one species label \n\n";
        out2 << " Number of sites      : "<<rows<<"x"<<cols<<"\n";
        out2 << " Number of iterations : "<<nbiter<<"\n";
        out2 << " Duration of timestep : "<<timestep<<" years\n";
        out2 << " Number of Species    : "<<numesp << "\n\n";
        out2.flush();
    }
    
    
    /*** Initialization of trees ***/
    /*******************************/
    
    if(NULL==(T=new Tree[sites])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Tree" << endl;
    }
    
#ifdef NDD
    for(int site=0;site<sites;site++) {
        if (NULL==(T[site].t_NDDfield = new float[numesp+1])) cerr<<"!!! Mem_Alloc\n";
        for(int ii=0;ii<(numesp+1);ii++) T[site].t_NDDfield[ii]=0;
    }
#endif
    
    
    /*** Initialization of species ***/
    /*********************************/
    
    int sp;
    if(NULL==(S=new Species[numesp+1])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Species" << endl;
    }
    
    for(ligne=0;ligne<3;ligne++) In.getline(buffer,128,'\n');                           /* Read species parameters (ifstream In) */
    for(sp=1;sp<=numesp;sp++) {
        S[sp].Init(sp,In);
    }
    
    /*** Initialization of environmental variables ***/
    /*************************************************/
    
    In.getline(buffer,128,'\n');
    In.getline(buffer,128,'\n');
    In.getline(buffer,128,'\n');
    
    if(NULL==(Temperature=new float[iterperyear])) {                                // rk, the current structure of code suppose that environment is periodic (a period = a year), if one want to make climate vary, with interannual variation and climate change along the simulation, one just need to provide the full climate input of the whole simulation (ie number of columns=iter and not iterperyear) and change iterperyear by nbiter here.
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Temperature" << endl;
    }
    
    for (int i=0; i<iterperyear; i++) In >> Temperature[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "Temperature" << i << "\t"  << Temperature[i] <<  "\t";
    //cout<<endl;
    
    if(NULL==(DailyMaxTemperature=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc DailyMaxTemperature" << endl;
    }
    for (int i=0; i<iterperyear; i++) In >> DailyMaxTemperature[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "DailyMaxTemperature" << i << "\t"  << DailyMaxTemperature[i] << "\t";
    //cout<<endl;
    
    if(NULL==(NightTemperature=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc NightTemperature" << endl;
    }
    for (int i=0; i<iterperyear; i++) In >> NightTemperature[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "NightTemperature" << i << "\t"  << NightTemperature[i] << "\t";
    //cout<<endl;
    
    if(NULL==(Rainfall=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Rainfall" << endl;
    }
    for (int i=0; i<iterperyear; i++) In >> Rainfall[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "Rainfall" << i << "\t"  << Rainfall[i] << "\t";
    //cout<<endl;
    
    if(NULL==(WindSpeed=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc WindSpeed" << endl;
    }
    for (int i=0; i<iterperyear; i++) In >> WindSpeed[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "WindSpeed" << i << "\t"  << WindSpeed[i] << "\t";
    //cout<<endl;
    
    if(NULL==(MaxIrradiance=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Irradiance" << endl;
    }
    for (int i=0; i<iterperyear; i++) In >> MaxIrradiance[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "MaxIrradiance" << i << "\t"  << MaxIrradiance[i] << "\t";
    //cout<<endl;
    
    if(NULL==(MeanIrradiance=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc Irradiance" << endl;
    }
    for (int i=0; i<iterperyear; i++) In >> MeanIrradiance[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "MeanIrradiance" << i << "\t"  << MeanIrradiance[i] << "\t";
    //cout<<endl;
    
    if(NULL==(SaturatedVapourPressure=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc SaturatedVapourPressure" << endl;
    }
    for (int i=0; i<iterperyear; i++) In >> SaturatedVapourPressure[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "SaturatedVapourPressure" << i << "\t"  << SaturatedVapourPressure[i] << "\t";
    //cout<<endl;
    
    if(NULL==(VapourPressure=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc VapourPressure" << endl;
    }
    for (int i=0; i<iterperyear; i++) In >> VapourPressure[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "VapourPressure" << i << "\t"  << VapourPressure[i] << "\t";
    //cout<<endl;
    
    if(NULL==(VapourPressureDeficit=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc VapourPressureDeficit" << endl;
    }
    for (int i=0; i<iterperyear; i++) In >> VapourPressureDeficit[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "VapourPressureDeficit" << i << "\t"  << VapourPressureDeficit[i] << "\t";
    //cout<<endl;
    
    if(NULL==(DailyVapourPressureDeficit=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc DailyVapourPressureDeficit" << endl;
    }
    for (int i=0; i<iterperyear; i++) In >> DailyVapourPressureDeficit[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "DailyVapourPressureDeficit" << i << "\t"  << DailyVapourPressureDeficit[i] << "\t";
    //cout<<endl;
    
    if(NULL==(DailyMaxVapourPressureDeficit=new float[iterperyear])) {
        cerr<<"!!! Mem_Alloc\n";
        cout<<"!!! Mem_Alloc DailyMaxVapourPressureDeficit" << endl;
    }
    for (int i=0; i<iterperyear; i++) In >> DailyMaxVapourPressureDeficit[i];
    In.getline(buffer,128,'\n');
    
    //for (int i=0; i<iterperyear; i++) cout<< "DailyMaxVapourPressureDeficit" << i << "\t"  << DailyMaxVapourPressureDeficit[i] << "\t";
    //cout<<endl;
    
    temp=Temperature[iter%iterperyear];
    tmax=DailyMaxTemperature[iter%iterperyear];
    tnight=NightTemperature[iter%iterperyear];
    precip=Rainfall[iter%iterperyear];
    WS=WindSpeed[iter%iterperyear];
    Wmax=MaxIrradiance[iter%iterperyear]*1.678;       // 1.678 is to convert irradiance from W/m2 to micromol of PAR /s /m2, the unit used in the FvCB model of photosynthesis
    Wmean=MeanIrradiance[iter%iterperyear];            // still in W/m2
    e_s=SaturatedVapourPressure[iter%iterperyear];
    e_a=VapourPressure[iter%iterperyear];
    VPD=VapourPressureDeficit[iter%iterperyear];
    VPDday=DailyVapourPressureDeficit[iter%iterperyear];
    VPDmax=DailyMaxVapourPressureDeficit[iter%iterperyear];
    
    In.close();
    /* Close ifstream In */
    
    In.open(nomfi, ios::in);
    if(!mpi_rank) {
        do{
            In.getline(buffer,256,'\n');
            out << buffer <<endl;
        }while (!In.eof()) ;
    }
    
    In.close();
    /* Close ifstream In */
    
    
    
    /*** Initialization of output streams ***/
    /****************************************/
    
    char nnn[100];
    if(!mpi_rank) {
        
#ifdef outputMPI
        sprintf(nnn,"%s_%i_outputs.txt",buf, easympi_rank);
        sor[155].open(nnn, ios::out);
#else
        sprintf(nnn,"%s_%i_abund.txt",buf, easympi_rank);
        sor[0].open(nnn, ios::out);
        if (!sor[0]) {
            cout<< "ERROR with abund file"<< endl;
        }
        sprintf(nnn,"%s_%i_abu10.txt",buf, easympi_rank);
        sor[1].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_abu30.txt",buf, easympi_rank);
        sor[2].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_ba.txt",buf, easympi_rank);
        sor[3].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_npp.txt",buf, easympi_rank);
        sor[4].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_gpp.txt",buf, easympi_rank);
        sor[5].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_ba10.txt",buf, easympi_rank);
        sor[6].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_ppfd0.txt",buf, easympi_rank);
        sor[7].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_death.txt",buf, easympi_rank);
        sor[8].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_state.txt",buf, easympi_rank);
        sor[9].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_final_pattern.txt",buf, easympi_rank);
        sor[10].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_site1.txt",buf, easympi_rank);
        // used to be added to explore potential link with Belassen curve, could be suppressed, but maybe useful to have an idea of the magnitude and distribution of increment of dbh
        sor[11].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_site2.txt",buf, easympi_rank);
        sor[12].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_site3.txt",buf, easympi_rank);
        sor[13].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_site4.txt",buf, easympi_rank);
        sor[155].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_site5.txt",buf, easympi_rank);
        sor[156].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_site6.txt",buf, easympi_rank);
        sor[157].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_cica.txt",buf, easympi_rank);
        sor[14].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_sp_par.txt",buf, easympi_rank);
        sor[15].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_agb.txt",buf, easympi_rank);
        sor[16].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_100yearsofsolitude.txt",buf, easympi_rank);
        sor[17].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_NDDfield.txt",buf, easympi_rank);
        sor[142].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_Rday.txt",buf, easympi_rank);
        sor[143].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_Rnight.txt",buf, easympi_rank);
        sor[144].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_Rstem.txt",buf, easympi_rank);
        sor[145].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_death1.txt",buf, easympi_rank);
        sor[146].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_death2.txt",buf, easympi_rank);
        sor[147].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_death3.txt",buf, easympi_rank);
        sor[148].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_deathrate.txt",buf, easympi_rank);
        sor[149].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_litterfall.txt",buf, easympi_rank);
        sor[150].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_leafdens1.txt",buf, easympi_rank);
        sor[151].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_leafdens2.txt",buf, easympi_rank);
        sor[152].open(nnn, ios::out);
        sprintf(nnn,"%s_%i_leafdens3.txt",buf, easympi_rank);
        sor[153].open(nnn, ios::out);
        
        sprintf(nnn,"%s_%i_dbh.txt",buf, easympi_rank);
        sch[1].open(nnn,ios::out);
        sprintf(nnn,"%s_%i_vertd.txt",buf, easympi_rank);
        sch[2].open(nnn,ios::out);
        
#endif
        
        sprintf(nnn,"%s_%i_paramspace.txt",buf, easympi_rank);
        sor[154].open(nnn, ios::out);
        
        
        sor[154] << "proc"  << "\t" << easympi_rank << endl;
        sor[154] << "phi"  << "\t" << phi << endl;
        sor[154] << "k"  << "\t" <<klight << endl;
        sor[154] << "fallocwood"  << "\t" << fallocwood << endl;
        sor[154] << "falloccanopy"  << "\t" << falloccanopy << endl;
        sor[154] << "m"  << "\t" << m << endl;
        sor[154] << "m1"  << "\t" << m1 << endl;
    }
}


/***************************************
 *** Field dynamic memory allocation ***
 ***************************************/

void AllocMem() {
    
    int spp,haut,i;
    
    /*HEIGHT = 0;
     Maximal geometric constants
     for(spp=1;spp<=numesp;spp++) {
     HEIGHT = max(HEIGHT,int(S[spp].s_hmax*S[spp].s_dmax/(S[spp].s_dmax+S[spp].s_ah)));    in number of vertical cells
     }*/
    HEIGHT = 80;
    // try 12/05/2015, in order to enable big trees still to grow
    
    float d,r;
    d = 0.0;
    r = 0.0;
    for(spp=1;spp<=numesp;spp++){
        d = maxf(d,S[spp].s_dmax);
        /* in number of horizontal cells */
        r = maxf(r,ra0+S[spp].s_dmax*ra1);
        /* in number of horizontal cells */
    }
    
    RMAX = int(r+p*NH*LV*HEIGHT);
    //  RMAX = int(r);
    SBORD = cols*RMAX;
    dbhmaxincm = int(100.*d);
    
    if(!mpi_rank) {
        /*cout << "HEIGHT : " << HEIGHT
         << " RMAX : " << RMAX << " DBH : " << DBH <<"\n";
         cout.flush();*/
        
        if(RMAX>rows){
            /* Consistency tests */
            cerr << "Erreur : RMAX > rows \n";
            exit(-1);
        }
        if(HEIGHT > rows){
            cerr << "Erreur : HEIGHT > rows \n";
            exit(-1);
        }
    }
    
    
    /*** Initialization of dynamic Fields ***/
    /****************************************/
    
    if (NULL==(nbdbh=new int[dbhmaxincm])) cerr<<"!!! Mem_Alloc\n";                         /* Field for DBH histogram */
    if (NULL==(layer=new float[HEIGHT+1])) cerr<<"!!! Mem_Alloc\n";                         /* Field for variables averaged by vertical layer */
    if (NULL==(ESP_N=new int[numesp+1])) cerr<<"!!! Mem_Alloc\n";                           /* Field for democratic seed germination */
#ifdef SEEDTRADEOFF
    if (NULL==(PROB_S=new float[numesp+1])) cerr<<"!!! Mem_Alloc\n";
#endif
    //  if (NULL==(persist=new long int[nbiter])) cerr<<"!!! Mem_Alloc\n";                  /* Field for persistence */
    //  if (NULL==(distr=new int[cols])) cerr<<"!!! Mem_Alloc\n";
    
    if (NULL==(LAI3D=new float*[HEIGHT+1]))                                                   /* Field 3D */
        cerr<<"!!! Mem_Alloc\n";                                                            /* Trees at the border of the simulated forest need to know the canopy occupancy by trees in the neighboring processor.*/
    for(haut=0;haut<(HEIGHT+1);haut++)                                                          /* For each processor, we define a stripe above (labelled 0) and a stripe below (1). Each stripe is SBORD in width.*/
        if (NULL==(LAI3D[haut]=new float[sites+2*SBORD]))                                   /* ALL the sites need to be updated.*/
            cerr<<"!!! Mem_Alloc\n";
    for(haut=0;haut<(HEIGHT+1);haut++)
        for(int site=0;site<sites+2*SBORD;site++)
            LAI3D[haut][site] = 0.0;
    
    if (NULL==(Thurt[0]=new unsigned short[3*sites]))                                       /* Field for treefall impacts */
        cerr<<"!!! Mem_Alloc\n";
    for(i=1;i<3;i++)
        if (NULL==(Thurt[i]=new unsigned short[sites]))
            cerr<<"!!! Mem_Alloc\n";
    
#ifdef MPI                                                                                      /* Fields for MPI operations */
    for(i=0;i<2;i++){                                                                       /*  Two fields: one for the CL north (0) one for the CL south (1) */
        if (NULL==(LAIc[i]=new unsigned short*[HEIGHT+1]))                                    /* These fields contain the light info in the neighboring procs (2*SBORD in width, not SBORD !). They are used to update local fields */
            cerr<<"!!! Mem_Alloc\n";
        for(haut=0;haut<(HEIGHT+1);haut++)
            if (NULL==(LAIc[i][haut]=new unsigned short[2*SBORD]))
                cerr<<"!!! Mem_Alloc\n";
    }
#endif
}


/***************************************
 **** Initial non-local germination ****
 ***************************************/

void BirthInit() {
    
    //int site,temp2;
    int temp1;
    
#if INHOMOGENE
    if(!p_rank)
#endif
        
        nblivetrees=0;
    
    cout << "Initial germination" << endl;
    
    for(int spp=1;spp<=numesp;spp++) {                              /* Initial condition: intial germination due to external seed rain (of s_nbext seeds for each species) */
        
        for(int ii=0;ii<S[spp].s_nbext;ii++){
            temp1 = genrand2i()%sites;
            if(T[temp1].t_age==0) T[temp1].Birth(S,spp,temp1);
        }
        
        cout << S[spp].s_name << "\t" << S[spp].s_nbind <<  endl;
        
    }
    
    cout << "\tNBTrees\tinitial\tgermination" << nblivetrees << endl;
    
    
    /* for(int site=0;site<sites;site++)                             Initial condition: seeds (of species 1 or 2 only, randomly) everywhere (loop over all sites)
     T[site].Birth(S,genrand2i()%2+1,site);
     
     for(int spp=1;spp<=numesp;spp++) {
     cout << S[spp].s_name << "\t" << S[spp].s_nbind << endl;
     }
     
     cout << "\tNBTrees\t" << nblivetrees << endl;
     
     */
    
    /* for(int row=0;row<rows;row++)                                  Initial condition: two zones (one side filled with species 1, the other with species 2)
     for(int col=0;col<cols;col++){
     site = col+cols*row;
     if(col<cols/2)
     T[site].Birth(S,1,site);
     else
     T[site].Birth(S,2,site);
     }*/
    
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    cout<<endl;
}


/*###############################################
 ###############################################
 #######  Evolution at each timestep    ########
 ###############################################
 ###############################################*/

void Evolution() {
    
#ifdef DISTURBANCE
    Disturbance();                      /* Create a disturbance for a given iteration */
#endif
    UpdateField();                      /* Update light fields and seed banks */
    UpdateTree();                       /* Update trees */
#ifdef BASICTREEFALL
    UpdateTreefall();                   /* Compute and distribute Treefall events */
#endif
    Average();                          /* Compute averages for outputs */
    OutputField();                      /* Output the statistics */
}


/*##################################
 ####    Compute a disturbance   ###
 ##################################*/

void Disturbance() {

	// Two issues:
	// 1 - Can the function go two times on the same site and count two time the same dbh ?
	// 2 - If the last tree selected dbh explode the disturbance intensity (need to add real disturbance intensity output ?)

    int site;
    float dbh=0.0, disturb_dbh=0.0;

    for(site=1;site<=sites;site++)
        dbh += T[site].t_dbh;  

    if(iter == disturb_iter) {
        cout << "Disturbance of " << disturb_intensity * 100 << "% of BA." << endl;

        while (disturb_dbh/dbh < disturb_intensity) {
            int site=floor(genrand2()*sites);
            disturb_dbh += T[site].t_dbh;
            T[site].t_hurt += 10*T[site].t_Tree_Height;
            // cout << "Site " << site << " loose is tree of dbh " << T[site].t_dbh << " disturbed dbh is now " << disturb_dbh << endl;
        }
    }

}


/*##################################
 ####    Compute field LAI 3D    ###
 ####    Compute field Seed     ####
 ##################################*/

void UpdateField() {
    
    int site,haut,spp;
    
    /***  Compute Field LAI3D  ***/
    /*****************************/
    
    
#ifdef MPI
    /* Reinitialize field LAI3D */
    for(i=0;i<2;i++)
        for(haut=0;haut<(HEIGHT+1);haut++)
            for(site=0;site<2*SBORD;site++)
                LAIc[i][haut][site] = 0;
#endif
    
    for(haut=0;haut<(HEIGHT+1);haut++)
        for(site=0;site<sites+2*SBORD;site++)
            LAI3D[haut][site] = 0.0;
    
    for(site=0;site<sites;site++)                                    /* Each tree contribues to LAI3D */
        T[site].CalcLAI();
    
    int sbsite;
    //int cinde;
    for(haut=HEIGHT;haut>0;haut--)                                 /* LAI is computed by summing LAI from the canopy top to the ground */
        for(site=0;site<sites;site++){
            sbsite=site+SBORD;
#ifdef ONLYVERTICAL                                                      /* Several choices of light incidence */
            LAI3D[haut-1][sbsite] += LAI3D[haut][sbsite];
#else
            row0=site/cols;
            col0=site%cols;
            temp=0.0;
            cinde=0;
            if(col0>0) { cinde++;temp+=LAI3D[haut][sbsite-1]; }
            if(col0<cols-1) { cinde++;temp+=LAI3D[haut][sbsite+1]; }
            if(row0>0) { cinde++;temp+=LAI3D[haut][sbsite-cols]; }
            if(row0<rows-1) { cinde++;temp+=LAI3D[haut][sbsite+cols]; }
            //site2 = col+cols*row+SBORD;
            //temp=LAI3D[haut][sbsite-1]+LAI3D[haut][sbsite-cols]+LAI3D[haut][sbsite+1]+LAI3D[haut][sbsite+cols];
            //if((site<SBORD) || (site >= sites+SBORD)) temp =0;
            LAI3D[haut-1][sbsite] += p1*LAI3D[haut][sbsite]+p2*temp/float(cinde);
#endif
        }
    
    /*  for(site=0;site<sites;site++)       Trunks shade the cells in which they are
     T[site].TrunkLAI();*/
    
    
#ifdef MPI
    /* Communicate border of field */
    /*MPI_ShareField(LAI3D,LAIc,2*SBORD);
     This MPI command no longer exists in openMPI
     Action 20/01/2016 TODO: FIX THIS */
    MPI_ShareField(LAI3D,LAIc,2*SBORD);
    for(haut=0;haut<(HEIGHT+1);haut++){
        /* Add border effects in local fields */
        if(mpi_rank)
            for(site=0;site<2*SBORD;site++)
                LAI3D[haut][site] += LAIc[0][haut][site];
        if(mpi_rank<mpi_size-1)
            for(site=0;site<2*SBORD;site++)
                LAI3D[haut][site+sites] += LAIc[1][haut][site];
    }
#endif
    
    /*** Evolution of the field Seed **/
    /*********************************/
    
    /* Pass seeds across processors => two more fields to be communicated between n.n. (nearest neighbor) processors.
     NB: dispersal distance is bounded by the value of 'rows'.
     At least 99 % of the seeds should be dispersed within the stripe or on the n.n. stripe.
     Hence rows > 4.7*max(dist_moy_dissemination),for an exponential dispersal kernel.*/
    
    for(site=0;site<sites;site++)                                   /* Seeds produced by mature trees */
        if(T[site].t_age)
            T[site].DisperseSeed();
    
#ifdef MPI                                                              /* Seeds passed across processors */
    for(spp=1;spp<=numesp;spp++) {
        MPI_ShareSeed(S[spp].s_Gc,sites);
        S[spp].AddSeed();
    }
#endif
    
#if INHOMOGENE
    if(!p_rank) {
#endif
        for(spp=1;spp<=numesp;spp++) {                              /* External seed rain: constant flux from the metacommunity */
            for(int ii=0;ii<S[spp].s_nbext;ii++){
                site = genrand2i()%sites;
#ifdef SEEDTRADEOFF
                S[spp].s_Seed[site]++;
#else
                if(S[spp].s_Seed[site]!=1) S[spp].s_Seed[site] = 1;
#endif
            }
        }
#if INHOMOGENE
    }
    else{
        for(spp=1;spp<=numesp;spp++) {
            if(S[spp].s_nbind*sites > 50000)
                for(int ii=0;ii<S[spp].s_nbext;ii++){
                    site = genrand2i()%sites;
#ifdef SEEDTRADEOFF
                    S[spp].s_Seed[site]++;
#else
                    if(S[spp].s_Seed[site]!=1) S[spp].s_Seed[site] = 1;
#endif
                }
        }
    }
#endif
    
#ifdef MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    
#ifdef NDD
    /*** Evolution of the field NDDfield **/
    /**************************************/
    
    float R=15;
    for(site=0;site<sites;site++) {
        for(spp=1;spp<=numesp;spp++) {
            if (iter == int(nbiter/10))  {
                sor[142]<< T[site].t_NDDfield[spp] << "\t" ;
            }
            T[site].t_NDDfield[spp]=0;
        }
        sor[142]<< "\n";
        
        row0=T[site].t_site/cols;
        col0=T[site].t_site%cols;
        for(col=max(0,col0-R);col<=min(cols-1,col0+R);col++) {
            for(row=max(0,row0-R);row<=min(rows-1,row0+R);row++) {                   /* loop over the neighbourhood */
                xx=col0-col;
                yy=row0-row;
                float d=sqrt(xx*xx+yy*yy);
                if((d<=R)&&(d>0))  {                                             /* is the voxel within the neighbourhood?  */
                    int j=cols*row+col;
                    if(T[j].t_age) T[site].t_NDDfield[T[j].t_sp_lab]+= PI* T[j].t_dbh * T[j].t_dbh*0.25/d;
                }
            }
        }
    }
#endif
}


/*##############################
 ####	Birth, Growth       ####
 ####	and death of trees  ####
 ##############################*/

void UpdateTree() {
    
    int site,spp,iii;
    float flux;
    
    for(site=0;site<sites;site++) {                                     /***** Local germination *****/
        if(T[site].t_age == 0) {
            iii=0;
#ifdef SEEDTRADEOFF
            float tot=0;
#endif
            for(spp=1;spp<=numesp;spp++){                               /* lists all the species with a seed present at given site... */
                if(S[spp].s_Seed[site]) {
                    ESP_N[iii]=spp;
                    
#ifdef SEEDTRADEOFF
#ifdef NDD
                    float p=S[spp].s_Seed[site]*S[spp].s_seedmass/(T[site].t_NDDfield[spp]*10000+1);
                    PROB_S[iii]=tot+ p;
                    tot+=p;
#else
                    PROB_S[iii]=tot+ S[spp].s_Seed[site]*S[spp].s_seedmass;
                    tot+=S[spp].s_Seed[site]*S[spp].s_seedmass;
#endif
#endif
                    iii++;
                }
            }
            if(iii) {                                                   /* ... and then randomly select one of these species */
#ifdef SEEDTRADEOFF
                double p=genrand2();                                    /* if SEEDTRADEOFF is define, probability of species recruit are proportional to the number of seeds time the seed mass, if NDD is also defnes the probablility is also inversly proportional to the species NDDfield term at that site */
                float itot=1/tot;
                int s=0;
                while (p>PROB_S[s]*itot) {s++;}
                spp=ESP_N[s];
#else
                spp = ESP_N[rand()%iii];
                /* otherwise species are equiprobable */
#endif
                flux = Wmax*exp(-flor(LAI3D[0][site+SBORD])*klight);
                if(flux>(S[spp].s_LCP)){
                    /* If enough light, germination, initialization of NPP (LCP is the species light compensation point*/
                    /* here, light is the sole environmental resources tested as a limiting factor for germination, but we should think about adding nutrients (N,P) and water conditions... */
                    T[site].Birth(S,spp,site);
                }
            }
        }
    }
    
    nbm_n1=nbm_n10=0;
#ifdef BASICTREEFALL
    nbm_c1=nbm_c10=0;
#endif
    for(site=0;site<sites;site++) {
        /***** Tree evolution: Growth or death *****/
        T[site].Update();
        /***** Update the Seed field *****/
        for(spp=1;spp<=numesp;spp++)
            S[spp].UpdateSeed(T[site].t_age,site);
    }
}


/******************************
 *** Treefall gap formation ***
 ******************************/

#ifdef BASICTREEFALL
void UpdateTreefall(){
    int site;
    
    for(site=0;site<sites;site++){
        Thurt[0][site] = Thurt[0][site+2*sites] = 0;
        Thurt[0][site+sites] = 0;
    }
    
    nbTreefall1 = 0;
    nbTreefall10=0;
    for(site=0;site<sites;site++)
        if(T[site].t_age) {
            T[site].FallTree();
        }
    
    
#ifdef MPI
    /* Treefall field passed to the n.n. procs */
    MPI_ShareTreefall(Thurt, sites);
#endif
    
    for(site=0;site<sites;site++)
    /* Update of Field hurt */
        if(T[site].t_age) {
            T[site].t_hurt = Thurt[0][site+sites];
#ifdef MPI
            if(mpi_rank) T[site].t_hurt = max(T[site].t_hurt,Thurt[1][site]);
            if(mpi_rank<mpi_size-1) T[site].t_hurt = max(T[site].t_hurt,Thurt[2][site]);
#endif
        }
}
#endif


/*###############################################
 ###############################################
 #######        Output routines         ########
 ###############################################
 ###############################################*/

/*********************************************************
 *** Calculation of the global averages every timestep ***
 *********************************************************/

void Average(void){
    
    int site,spp,i;
    float sum1=0,sum10=0.0,sum30=0.0,ba=0.0,npp=0.0,gpp=0.0, ba10=0.0, agb=0.0, rday=0.0, rnight=0.0, rstem=0.0, litterfall=0.0;
    
    if(!mpi_rank) {
        for(spp=1;spp<=numesp;spp++)
            for(i=0;i<12;i++)
                S[spp].s_chsor[i]=0;
        
        float inbcells = 1.0/float(sites*mpi_size);
        float inbhectares = inbcells*NH*NH*10000.0;
        
#ifdef outputMPI
        sor[155] << iter << "\t";
        for(spp=1;spp<=numesp;spp++)
            sor[155] << float(S[spp].s_nbind)*inbhectares << "\t";
#else
        for(i=0;i<7;i++) sor[i] << iter << "\t";
        for(i=143;i<146;i++) sor[i] << iter << "\t";
        for(spp=1;spp<=numesp;spp++) sor[0] << float(S[spp].s_nbind)*inbhectares << "\t";
#endif
        for(site=0;site<sites;site++) T[site].Average();
        
        cout << iter << "\tNBtrees\t"<<nblivetrees << endl;
        
        for(spp=1;spp<=numesp;spp++) {
            S[spp].s_chsor[1] *= inbhectares;
            //species number of trees with dbh>10cm
            S[spp].s_chsor[2] *= inbhectares;
            //species number of trees with dbh>30cm
            S[spp].s_chsor[3] *= 3.1415*0.25*inbhectares;
            //species basal area
            S[spp].s_chsor[4] *= inbhectares;
            //species total NPP (sum of t_NPP) in MgC (per timestep)
            S[spp].s_chsor[5] *= inbhectares;
            //species total GPP (sum of t_GPP) in MgC (per timestep)
            S[spp].s_chsor[6] *= 3.1415*0.25*inbhectares;
            //species basal area; with only trees with dbh>10cm
            S[spp].s_chsor[7] *= inbhectares;
            //species aboveground biomass
            S[spp].s_chsor[8] *= inbhectares;
            /* species leaf Rday in MgC (per timestep) */
            S[spp].s_chsor[9] *= inbhectares;
            //species leaf Rnight in MgC (per timestep)
            S[spp].s_chsor[10] *= inbhectares;
            //species Rstem  in MgC (per timestep)
            S[spp].s_chsor[11] *= inbhectares;
            //species litterfall  in MgC (per timestep)
            sum1+= float(S[spp].s_nbind)*inbhectares;
            sum10 += S[spp].s_chsor[1];
            sum30 += S[spp].s_chsor[2];
            ba += S[spp].s_chsor[3];
            npp += S[spp].s_chsor[4];
            gpp += S[spp].s_chsor[5];
            ba10 += S[spp].s_chsor[6];
            agb += S[spp].s_chsor[7];
            rday += S[spp].s_chsor[8];
            rnight += S[spp].s_chsor[9];
            rstem += S[spp].s_chsor[10];
            litterfall += S[spp].s_chsor[11];
#ifdef outputMPI
            //nothing
#else
            for(i=1;i<7;i++) {
                sor[i] << S[spp].s_chsor[i] << "\t";
            }
            sor[16] << S[spp].s_chsor[7] << "\t";
            sor[143] << S[spp].s_chsor[8] << "\t";
            sor[144] << S[spp].s_chsor[9] << "\t";
            sor[145] << S[spp].s_chsor[10] << "\t";
            sor[150] << S[spp].s_chsor[11] << "\t";
#endif
        }
        
#ifdef outputMPI
        sor[155] << sum1 << "\t";
        //total number of trees (dbh>1cm=DBH0)
        for(spp=1;spp<=numesp;spp++) {
            sor[155] << S[spp].s_chsor[1] << "\t";
        }
        sor[155] << sum10 << "\t";
        //total number of trees (dbh>10cm=DBH0)
        
        for(spp=1;spp<=numesp;spp++) {
            sor[155] << S[spp].s_chsor[2] << "\t";
        }
        sor[155] << sum30 << "\t";
        //total number of trees (dbh>30cm=DBH0)
        for(spp=1;spp<=numesp;spp++) {
            sor[155] << S[spp].s_chsor[6] << "\t";
        }
        sor[155] << ba10 << "\t";
        for(spp=1;spp<=numesp;spp++) {
            sor[155] << S[spp].s_chsor[4] << "\t";
        }
        sor[155] << npp << "\t";
        for(spp=1;spp<=numesp;spp++) {
            sor[155] << S[spp].s_chsor[5] << "\t";
        }
        sor[155] << gpp << "\t";
        for(spp=1;spp<=numesp;spp++) {
            sor[155] << S[spp].s_chsor[7] << "\t";
        }
        sor[155] << agb << endl;
#else
        
        sor[0] << sum1 << endl;
        //total number of trees (dbh>1cm=DBH0)
        sor[1] << sum10 << endl;
        //total number of trees with dbh>10cm
        sor[2] << sum30 << endl;
        //total number of trees with dbh>30cm
        sor[3] << ba << endl;
        //total basal area
        sor[4] << npp << endl;
        //total NPP in MgC per ha (per time step)
        sor[5] << gpp << endl;
        //total GPP in MgC par ha (per time step)
        sor[6] << ba10 << endl;
        //total basal area with only trees with dbh>10cm
        sor[16] << agb << endl;
        //total above ground biomass
        sor[143] << rday << endl;
        //total leaf day respiration
        sor[144] << rnight << endl;
        //total leaf night respiration
        sor[145] << rstem << endl;
        //total stem respiration
        sor[150] << litterfall << endl;
        //total stem respiration
        
        float tototest=0.0, tototest2=0.0, flux;
        for(site=0;site<sites;site++) {
            flux = Wmax*exp(-flor(LAI3D[0][site+SBORD])*klight);
            tototest += flux;
            tototest2 += flux*flux;
        }
        tototest /=float(sites*LH*LH);                              // Average light flux (PPFD) on the ground
        tototest2 /=float(sites*LH*LH);
        if(iter)
            sor[7] << iter<< "\tMean PPFDground\t"<< tototest << "\t" << sqrt(tototest2-tototest*tototest) << "\n";
        
        sor[8] << iter << "\t" << nbm_n1*inbhectares << "\t" << nbm_n10*inbhectares<< "\t"
#ifdef BASICTREEFALL
        << nbm_c1*inbhectares << "\t" << nbm_c10*inbhectares << "\t" << nbTreefall1*inbhectares << "\t" << nbTreefall10*inbhectares << "\t"
#endif
        << endl;
#endif
    }
    
#ifdef MPI
    MPI_Reduce(&(S[spp].s_nbind),&sind,1,
               MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(S[spp].s_chsor,S[spp].s_chsor,5,
               MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(Mortality,Mortality,4,
               MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&S[spp].s_chsor[6],&S[spp].s_chsor[6],5,
               MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
#endif
    
    if(!mpi_rank) {
        if(iter == 200){                                        // State at the 200th iteration (for all living trees: dbh, height, crown radius and depth and dbh increment)
            for(site=0;site<sites;site++) {
                if(T[site].t_dbh>0.0)
                    sor[9] << T[site].t_dbh*LH*100 << "\t" << T[site].t_Tree_Height
                    << "\t" << T[site].t_Crown_Radius*LH
                    << "\t" << T[site].t_Crown_Depth*LV
                    << "\t" << T[site].t_ddbh*LH*100 << "\n";
            }
        }
    }
    
#ifdef MPI
    /* This section corresponds to the parallel version of
     the reporting of the global diagnostic variables. Since much work has been done
     on routine Average over the past years, this would need a full rewrite
     !!!!Action 20/01/2016: rework the parallel version of function Average!!!!
     
     MPI_Reduce(&(S[spp].s_nbind),&sind,1,
     MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(S[spp].s_chsor,S[spp].s_chsor,5,
     MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(Mortality,Mortality,4,
     MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&S[spp].s_chsor[6],&S[spp].s_chsor[6],5,
     MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
     */
#endif
    cout.flush();
}


/*********************************
 **** Output de Fields Shares ****
 *********************************/

void OutputField(){
    
    /*sch[0] vertical occupation of the voxel field
     
     sch[1] dbh histogram
     sch[2] mean LAI by height class
     */
    
    int site,haut;
    
    if((nbout)&&((iter%freqout)==freqout-1)) {
        // output fields, nbout times during simulation (every freqout iterations)
        
        int d;
        for(d=0;d<dbhmaxincm;d++) nbdbh[d]=0;
        for(site=0;site<sites;site++) T[site].histdbh();
        
        for(haut=0;haut<(HEIGHT+1);haut++){
            layer[haut] = 0;
            for(site=0;site<sites;site++) layer[haut] += LAI3D[haut][site+SBORD];
        }
        
#ifdef MPI
        MPI_Status status;
        MPI_Reduce(nbdbh,nbdbh,dbhmaxincm,
                   MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
        
        MPI_Reduce(layer,layer,HEIGHT,
                   MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
#endif
        
        if(!mpi_rank) {
            for(d=1;d<dbhmaxincm;d++)                                   // output of the dbh histograms (sch[1])
                sch[1] << d << "\t" << nbdbh[d]  << "\n";
            sch[1] <<  "\n";
            
            float norm = 1.0/float(sites*LH*LH*mpi_size);                   // output of the mean LAI per height class (sch[2])
            for(haut=0;haut<(HEIGHT+1);haut++)
                sch[2] << haut*LV << "\t" << layer[haut]*norm << "\n";
            sch[2] <<  "\n";
        }
    }
}


void SorCh(unsigned short * sort, int rrows, int ccols, fstream &oflux) {
    
    int col,row;
    //int ssites = rrows*ccols;
    
    /* Output for processor 0 */
    if(!mpi_rank)
        for(row=0;row<rrows;row++) {
            for(col=0;col<ccols;col++)
                oflux << sort[col+row*ccols] << " ";
            oflux << "\n";
        }
#ifdef MPI
    /* Envoi du processeur ran sur le processeur 0, puis ecriture */
    MPI_Status status;
    MPI_Barrier(MPI_COMM_WORLD);
    for(int ran=1;ran<size;ran++) {
        if(p_rank==ran)
            MPI_Send(sort,ssites,MPI_UNSIGNED_SHORT,
                     0,0,MPI_COMM_WORLD);
        if(!p_rank) {
            MPI_Recv(sort,ssites,MPI_UNSIGNED_SHORT,
                     ran,0,MPI_COMM_WORLD,&status);
            for(row=0;row<rrows;row++) {
                for(col=0;col<ccols;col++)
                    oflux << sort[col+row*ccols] << " ";
                oflux << "\n";
            }
        }
    }
#endif
    if(!mpi_rank)   oflux <<  "\n";
    
}

#ifdef MPI

/* Communication of border fields in the parallel version of the code */
/* Only if the MPI option has been enabled */
void MPI_ShareSeed(unsigned char **c, int n) {
    
    MPI_Status status;
    
    if(p_rank == size-1)
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_CHAR,size-2,0,
                     c[3],n,MPI_UNSIGNED_CHAR,0,0,MPI_COMM_WORLD,&status);
    if(p_rank == 0)
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_CHAR,size-1,0,
                     c[3],n,MPI_UNSIGNED_CHAR,1,0,MPI_COMM_WORLD,&status);
    if((p_rank) && (p_rank < size-1))
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_CHAR,p_rank-1,0,
                     c[3],n,MPI_UNSIGNED_CHAR,p_rank+1,0,MPI_COMM_WORLD,&status);
    
    if(p_rank == 0)
        MPI_Sendrecv(c[1],n,MPI_UNSIGNED_CHAR,1,1,
                     c[2],n,MPI_UNSIGNED_CHAR,size-1,1,MPI_COMM_WORLD,&status);
    if(p_rank == size-1)
        MPI_Sendrecv(c[1],n,MPI_UNSIGNED_CHAR,0,1,
                     c[2],n,MPI_UNSIGNED_CHAR,size-2,1,MPI_COMM_WORLD,&status);
    if((p_rank) && (p_rank < size-1))
        MPI_Sendrecv(c[1],n,MPI_UNSIGNED_CHAR,p_rank+1,1,
                     c[2],n,MPI_UNSIGNED_CHAR,p_rank-1,1,MPI_COMM_WORLD,&status);
}

void MPI_ShareField(unsigned short **cl, unsigned short ***cp, int n) {
    
    MPI_Status status;
    for(int haut=0;haut<(HEIGHT+1);haut++) {
        if(p_rank == 0)
            MPI_Sendrecv(cl[haut],n,MPI_UNSIGNED_SHORT,size-1,haut,
                         cp[1][haut],n,MPI_UNSIGNED_SHORT,1,haut,
                         MPI_COMM_WORLD,&status);
        if(p_rank == size-1)
            MPI_Sendrecv(cl[haut],n,MPI_UNSIGNED_SHORT,size-2,haut,
                         cp[1][haut],n,MPI_UNSIGNED_SHORT,0,haut,
                         MPI_COMM_WORLD,&status);
        if((p_rank) && (p_rank < size-1))
            MPI_Sendrecv(cl[haut],n,MPI_UNSIGNED_SHORT,p_rank-1,haut,
                         cp[1][haut],n,MPI_UNSIGNED_SHORT,p_rank+1,haut,
                         MPI_COMM_WORLD,&status);
        
        if(p_rank == 0)
            MPI_Sendrecv(cl[haut]+sites,n,MPI_UNSIGNED_SHORT,1,haut+HEIGHT,
                         cp[0][haut],n,MPI_UNSIGNED_SHORT,size-1,haut+HEIGHT,
                         MPI_COMM_WORLD,&status);
        if(p_rank == size-1)
            MPI_Sendrecv(cl[haut]+sites,n,MPI_UNSIGNED_SHORT,0,haut+HEIGHT,
                         cp[0][haut],n,MPI_UNSIGNED_SHORT,size-2,haut+HEIGHT,
                         MPI_COMM_WORLD,&status);
        if((p_rank) && (p_rank < size-1))
            MPI_Sendrecv(cl[haut]+sites,n,MPI_UNSIGNED_SHORT,p_rank+1,haut+HEIGHT,
                         cp[0][haut],n,MPI_UNSIGNED_SHORT,p_rank-1,haut+HEIGHT,
                         MPI_COMM_WORLD,&status);
    }
}

void MPI_ShareTreefall(unsigned short **c, int n) {
    
    MPI_Status status;
    if(p_rank == 0)
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_SHORT,size-1,0,
                     c[2],n,MPI_UNSIGNED_SHORT,1,0,MPI_COMM_WORLD,&status);
    if(p_rank == size-1)
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_SHORT,size-2,0,
                     c[2],n,MPI_UNSIGNED_SHORT,0,0,MPI_COMM_WORLD,&status);
    if((p_rank) && (p_rank < size-1))
        MPI_Sendrecv(c[0],n,MPI_UNSIGNED_SHORT,p_rank-1,0,
                     c[2],n,MPI_UNSIGNED_SHORT,p_rank+1,0,MPI_COMM_WORLD,&status);
    
    if(p_rank == 0)
        MPI_Sendrecv(c[0]+2*n,n,MPI_UNSIGNED_SHORT,1,1,
                     c[1],n,MPI_UNSIGNED_SHORT,size-1,1,MPI_COMM_WORLD,&status);
    if(p_rank == size-1)
        MPI_Sendrecv(c[0]+2*n,n,MPI_UNSIGNED_SHORT,0,1,
                     c[1],n,MPI_UNSIGNED_SHORT,size-2,1,MPI_COMM_WORLD,&status);
    if((p_rank) && (p_rank < size-1))
        MPI_Sendrecv(c[0]+2*n,n,MPI_UNSIGNED_SHORT,p_rank+1,1,
                     c[1],n,MPI_UNSIGNED_SHORT,p_rank-1,1,MPI_COMM_WORLD,&status);
}

#endif


/******************************************
 ******************************************
 *******  Free dynamic memory  ************
 ******************************************
 ******************************************/

void FreeMem () {
    
    delete [] T;
    delete [] S;
    delete [] nbdbh;
    delete [] layer;
    delete [] ESP_N;
#ifdef SEEDTRADEOFF
    delete [] PROB_S;
#endif
    
    int haut;
    for (haut=0; haut<(HEIGHT+1); haut++) {
        delete [] LAI3D[haut];
    }
    delete [] LAI3D;
    
    int i;
    for (i=0; i<3; i++) {
        delete [] Thurt[i];
    }
}


/***********************************
 ***********************************
 ***** RANDOM NUMBER GENERATOR *****
 ***********************************
 ***********************************/

/* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
/* When you use this, send an email to: matumoto@math.keio.ac.jp   */
/* with an appropriate reference to your work.                     */

//#include<stdio.h>

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializing the array with a NONZERO seed */
void sgenrand2(unsigned long seed)
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

/* generating reals */
/* unsigned long */ /* for integer generation */
double genrand2()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    
    if (mti >= N) { /* generate N words at one time */
        int kk;
        
        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand2(4357); /* a default initial seed is used   */
        
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];
        
        mti = 0;
    }
    
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    
    return ( (double)y / (unsigned long)0xffffffff ); /* reals */
    /* return y; */ /* for integer generation */
}


/* initializing the array with a NONZERO seed */
void sgenrand2i(unsigned long seed)
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

unsigned long genrand2i()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    
    if (mti >= N) { /* generate N words at one time */
        int kk;
        
        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand2i(4357); /* a default initial seed is used   */
        
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];
        
        mti = 0;
    }
    
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    
    return y;
}


