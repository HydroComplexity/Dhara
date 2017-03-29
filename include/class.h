#ifndef PROJECTCLASS_H
#define PROJECTCLASS_H
    class ProjectClass
    {
        public:
        const char *name;
        const char *folderoutput;
        const char *mlcanoutput;
        const char *olfoutput;
        const char *ssfoutput;
        const char *ssf1doutput;
        const char *statoutput;
        
        int co2elevation;
        int numberofplants;
        int planttype;
        int savemlcan;
        int saveolf;
        int savessf;
        int savestat;
        int saveinterval;
        int printiterval;
        int num_steps;              // Number of simulation steps

        int verbose;                // Force to verbose mode
        int version;

        double co2concentration;
        double dtimehr;
        double dx_meter;
        double dy_meter;
        double dz_meter;
    };
#endif


#ifndef MPICLASS_H
#define MPICLASS_H
    class mpiClass
    {
        public:
        int *neighbors;             // four neighbors of each global MPI in 2D coordinate
        int *xdispls;               // displacements for send/recv
        int *ydispls;                     
        int *xcounts;               // counts for send/recv 
        int *ycounts;       
        int3 global_size;           // size of the global domain
        int3 domain_size;           // size of the local domains
        int2 offset;                // offset of the local domains
        int2 topology_size;         // size of topology nodes (mpi)
        int2 topology_index;        // indices of topology nodes (mpi)
        MPI_Datatype blocktype;
        MPI_Comm cartComm;

    };
#endif

#ifndef FILENAMECLASS_H
#define FILENAMECLASS_H
    class FileNameClass
    {
        public:
        const char *config;
        const char *topography;     // Topography filename
        const char *forcings;       // Forcings filename
        const char *plants;
        const char *ovl_output;
        const char *ssf_output;
        const char *mlcan_output;
    };
#endif


#ifndef TIMEFORCINGCLASS_H
#define TIMEFORCINGCLASS_H
    class TimeForcingClass
    {
        public:
        int *doy;
        int *years;
        double *rg;
        double *pa;
        double *lwdn;
        double *zen;
        double *u;
        double *ppt;
        double *ta;
        double *ea;
        double *lai;
        double *vpd;
        double *decdoy;
        double *hour;
    };
#endif


#ifndef OVERLANDFLOWCLASS_H
#define OVERLANDFLOWCLASS_H
    class OverlandFlowClass
    {
        public:
        double *a2d;
        double *hpoten;
        double *ke;
        double *kn;
        double *ks;
        double *kw;
        double *mann;
        double *ph;
        double *rhs2d;
        double *qcapa;
        double *u;
        double *v;
        double *waterdepth;
        double *waterelev;
        double *ztopo;

        double delta;
        double dx;
        double dy;
        double dt;        
        double hmin;
        double hcri;
        double K0;
        double hn;
        double hs;
        double hw;
        double he;
    };
#endif


#ifndef SUBSURFACEFLOWCLASS_H
#define SUBSURFACEFLOWCLASS_H
    class SubsurfaceFlowClass
    {
        public:
        /* BC types. */
        int *bcb;               // Bottom
        int *bce;               // East
        int *bcn;               // North
        int *bcs;               // South
        int *bct;               // Top
        int *bcw;               // West
        int *type;
        int *procmap;           // Map of MPI processes

        /* BC values */
        double *bcpsib;         // pressure bottom
        double *bcpsie;         // pressure east
        double *bcpsin;         // pressure north
        double *bcpsis;         // pressure south
        double *bcpsit;         // pressure top
        double *bcpsiw;         // pressure west
        double *bcqb;           // flux bottom
        double *bcqe;           // flux east
        double *bcqn;           // flux north
        double *bcqs;           // flux south
        double *bcqt;           // flux top
        double *bcqw;           // flux west

        double *a3d;            // Left hand side matrix A
        double *deltam;         // difference between 2 iters 
        double *cnp1m;          // specific moisture capacity at n+1,m
        double *knp1m;          // hydraulic conductivity at n+1,m
        double *ksat;           // saturated hydraulic conductivity
        double *psin;           // pressure head input
        double *psinp1m;        // pressure head at n+1,m
        double *psinp1mp1;      // pressure head at n+1,m+1
        double *psiout;         // pressure head output
        double *qss;            // infiltration flux q
        double *rhs3d;          // right hand side b
        double *thetan;         // theta input
        double *thetanp1m;      // theta at n+1,m
        double *thetanp1mp1;    // theta at n+1,m+1
        double *thetaout;       // theta output

        double *quflux;         // flux up [m/dtime] - go deeper the soil
        double *qdflux;         // flux down [m/dtime]
        double *qwflux;         // flux west [m/dtime]
        double *qeflux;         // flux east [m/dtime]
        double *qsflux;         // flux south [m/dtime]
        double *qnflux;         // flux north [m/dtime]
        double *dtheta;         // soil moisture difference
        double *transp;         // transpiration [m/dtime]
        double *evapo;          // evapotration [m/dtime]
        double *ssflux;         // flux storage [m/dtime]
        double *mb_subsurfaceW; // subsurface water balance

        double *TR;             // Transpiration entire domain
        double *TR_root;        // Transpiration gather MPI
        double *ppt_ground;     // Throughfall entire domain
        double *ppt_root;       // Throughfall gather MPI
        double *E_soil;         // Evapotration entire domain
        double *E_soil_root;    // Evapotration gather MPI
        double *rda;            // Root density from all process
        double *psi_col;        // column average
        double *theta_col;

        int maxiter;            // max iteration
        int picardmethod;       // 0 or 1
        double alpha;
        double am;
        double airdry;
        double dx;
        double dy;
        double dz;
        double dt;
        double poresize;
        double porosity;
        double psimin;
        double Ss;
        double theta_s;
        double theta_r;
        double tolerance_psi;   // stopping criteria
        double tolerance_theta;
    };
#endif


#ifndef SWITCH_H
#define SWITCH_H
    class SwitchClass
    {
        public:
        int PH_type;                  // Number of sub time steps
        int Turbulence;
        int HydraulicRedistribution;
        int RootConductivity;
        int PrintSubStep;             // Number of time step for printing info
        int RootHydrauConduct;
        int SoilHeat;
        int Plotting;
        int ElevatedCO2;
        int LWequation;               // Longwave equation 1: without atmospheric correction, 2: with
    };
#endif


#ifndef FORCINGCLASS_H
#define FORCINGCLASS_H
    class ForcingClass
    {
        public:
        int doy;
        int years;
        double decdoy;
        double rg;
        double pa;
        double lwdn;
        double zen;
        double u;
        double ppt;
        double ta;
        double ea;
        double ca;
    };
#endif


#ifndef CONSTANTS_H
#define CONSTANTS_H
    class ConstantClass
    {
        public:
        double umoltoWm2;
        double Wm2toumol;
        double mmH2OtoMPa;
        double R;
        double R_kJ;
        double Lv;
        double Lv_g;
        double cp_mol;
        double cp_JkgK;
        double boltz;
        double vonk;
        double rho_dry_air;
        double grav;
        double timestep;
        double dtime;
    };
#endif


#ifndef CANOPY_H
#define CANOPY_H
    class CanopyClass
    {
        public:
        const char *name;
        const char *dists;
        int nl_can;
        double LEfact;
        double Hfact;
        double LWfact;
        double R_kJ;
        double hcan;
        double hobs;
        double z0;
        double d0;
        double leaftype;
        double ld;
        double lw;
        double Smax;
        double Ffact;
        double pptintfact;
    };
#endif


#ifndef SOILCLASS_H
#define SOILCLASS_H
    class SoilClass
    {
        public:
        int nl_soil;
        double clay;
        double Cd_soil;
        double dzs;
        double depths;
        double z0;
        double smpmin;
        double wimp;
        double K_rad;
        double K_axs;
        double scalek;
        double HC_air;
        double HC_liq;
        double HC_ice;
        double rho_liq;
        double rho_ice;
        double sand;
        double TK_liq;
        double TK_ice;
        double TK_iair;
        double Tf;
        double alphCN;
        double kpar_ax;
    };
#endif


#ifndef RADIATION_H
#define RADIATION_H
    class RadiationClass
    {
        public:
        double transmiss;
        double epsv;
        double epss;
        double epsa;
        double xx;
        double clump;
        double Kdf;
        double absorp_PAR;
        double absorp_NIR;
        double refl_PAR;
        double refl_NIR;
        double refl_soil;
        double trans_PAR;
        double trans_NIR;
    };
#endif


#ifndef PHOTOSYNTHESIS_H
#define PHOTOSYNTHESIS_H
    class PhotosynthesisClass
    {
        public:
        double ph_type;
        double Vmax_C4;
        double Rd_C4;
        double Q10_C4;
        double kk_C4;
        double theta_C4;
        double beta_C4;
        double al_C4;
        double kn_canopy;
        double ap;
        double bp;
        double beta_ph_C3;
        double Vcmax25_C3;
        double Vcmax25_fact;
        double Jmax25_C3;
        double Rd25;
        double Oi;        
    };
#endif


#ifndef RESPIRATION_H
#define RESPIRATION_H
    class RespirationClass
    {
        public:
        double Ro;
        double Q10;
    };
#endif


#ifndef STOMACONDUCT_H
#define STOMACONDUCT_H
    class StomaConductClass
    {
        public:
        double mslope;
        double bint;
        double sf;
        double psif;
        double Rp;
    };
#endif

#ifndef MICROENVIRONMENT_H
#define MICROENVIRONMENT_H
    class MicroEnvironmentClass
    {
        public:
        double Cd;
        double alph;
    };
#endif





