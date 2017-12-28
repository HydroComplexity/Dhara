#ifndef VERTCANOPY_H
#define VERTCANOPY_H
    class VerticalCanopyClass
    {
        public:
        double dzc;
        double *radabs_tot;
        double *radlost;
        double *radremain;
        double *SWout;
        double *LWout;
        double *PARout;
        double *NIRout;
        double *Evap_can;
        double *Ch2o_can;
        double *Sh2o_can;
        double *Sh2o_can_prev;
        double *An_can;
        double *LE_can;
        double *H_can;
        double *TR_can;
        double *Rnrad_can;
        double *mbw_can;        
        double *zhc;
        double *znc;
        double *LAD;
        double *Tl_sun;
        double *gsv_sun;
        double *Ci_sun;
        double *Tl_shade;
        double *gsv_shade;
        double *Ci_shade;
        double *TR;
        double *Sh2o_prof;
        double *LADnorm;
        double *LADz;
        double *LAIz;
        double *TAz;
        double *CAz;
        double *EAz;
        double *PAz;
        double *Uz;
        double *Vz;
        double *Km;
        double *TR_sun;
        double *TR_shade;
        double *wetfrac;
        double *dryfrac;
        double *Smaxz;
        double *PARabs_sun;
        double *PARabs_shade;
        double *NIRabs_sun;
        double *NIRabs_shade;
        double *PARabs_sun_lai;
        double *PARabs_shade_lai;
        double *fsun;
        double *fshade;
        double *LAIsun;
        double *LAIshade;
        double *taud;
        double *diffdn;
        double *diffup;
        double *SWabs_sun;
        double *SWabs_shade;
        double *LWabs_can;
        double *LWabs_sun;
        double *LWabs_shade;
        double *LWcan_out;
        double *LWsun_out;
        double *LWshade_out;
        double *LWemit_can;
        double *LWemit_sun;
        double *LWemit_shade;
        double *Rabs_sun;
        double *Rabs_shade;
        double *gbv_sun;
        double *gbh_sun;
        double *gbv_shade;
        double *gbh_shade;
        double *An_sun;
        double *An_shade;
        double *psil_sun;
        double *psil_shade;
        double *fsv_sun;
        double *fsv_shade;
        double *Tl;
        double *Ph_sun;
        double *Ph_shade;
        double *LE_sun;
        double *LE_shade;
        double *H_sun;
        double *H_shade;
        double *Evap_sun;
        double *Evap_shade;
        double *Evap_prof;
        double *Ch2o_sun;
        double *Ch2o_shade;
        double *Ch2o_prof;
        double *Rnrad_sun;
        double *Rnrad_shade;
        VerticalCanopyClass (int size)
        {
            radabs_tot       = new double[1];
            radlost          = new double[1];
            radremain        = new double[1];
            SWout            = new double[1];
            LWout            = new double[1];
            PARout           = new double[1];
            NIRout           = new double[1];
            Evap_can         = new double[1];
            Ch2o_can         = new double[1];
            Sh2o_can         = new double[1];
            Sh2o_can_prev    = new double[1];
            NIRout           = new double[1];
            An_can           = new double[1];
            LE_can           = new double[1];
            H_can            = new double[1];
            TR_can           = new double[1];
            Rnrad_can        = new double[1];
            mbw_can          = new double[1];
            zhc              = new double[size];
            znc              = new double[size];
            LAD              = new double[size];
            Tl_sun           = new double[size];
            gsv_sun          = new double[size];
            Ci_sun           = new double[size];
            Tl_shade         = new double[size];
            gsv_shade        = new double[size];
            Ci_shade         = new double[size];
            TR               = new double[size];
            Sh2o_prof        = new double[size];
            LADnorm          = new double[size];
            LADz             = new double[size];
            LAIz             = new double[size];
            TAz              = new double[size];
            CAz              = new double[size];
            EAz              = new double[size];
            PAz              = new double[size];
            Uz               = new double[size];
            Vz               = new double[size];
            Km               = new double[size+1];
            TR_sun           = new double[size];
            TR_shade         = new double[size];
            wetfrac          = new double[size];
            dryfrac          = new double[size];
            Smaxz            = new double[size];
            PARabs_sun       = new double[size];
            PARabs_shade     = new double[size];
            NIRabs_sun       = new double[size];
            NIRabs_shade     = new double[size];
            PARabs_sun_lai   = new double[size];
            PARabs_shade_lai = new double[size];
            fsun             = new double[size];
            fshade           = new double[size];
            LAIsun           = new double[size];
            LAIshade         = new double[size];
            taud             = new double[size];
            diffdn           = new double[size];
            diffup           = new double[size];
            SWabs_sun        = new double[size];
            SWabs_shade      = new double[size];
            LWabs_can        = new double[size];
            LWabs_sun        = new double[size];
            LWabs_shade      = new double[size];
            LWcan_out        = new double[size];
            LWsun_out        = new double[size];
            LWshade_out      = new double[size];
            LWemit_can       = new double[size];
            LWemit_sun       = new double[size];
            LWemit_shade     = new double[size];
            Rabs_sun         = new double[size];
            Rabs_shade       = new double[size];
            gbv_sun          = new double[size];
            gbh_sun          = new double[size];
            gbv_shade        = new double[size];
            gbh_shade        = new double[size];
            An_sun           = new double[size];
            An_shade         = new double[size];
            psil_sun         = new double[size];
            psil_shade       = new double[size];
            fsv_sun          = new double[size];
            fsv_shade        = new double[size];
            Tl               = new double[size];
            Ph_sun           = new double[size];
            Ph_shade         = new double[size];
            LE_sun           = new double[size];
            LE_shade         = new double[size];
            H_sun            = new double[size];
            H_shade          = new double[size];
            Evap_sun         = new double[size];
            Evap_shade       = new double[size];
            Evap_prof        = new double[size];
            Ch2o_sun         = new double[size];
            Ch2o_shade       = new double[size];
            Ch2o_prof        = new double[size];
            Rnrad_sun        = new double[size];
            Rnrad_shade      = new double[size];
        };
    };
#endif



#ifndef VERTSOIL_H
#define VERTSOIL_H
    class VerticalSoilClass
    {
        public:
        double *rpp_wgt;
        double *ppt_ground;
        double *pptrate_ground;
        double *SWabs_soil;
        double *Rnrad_soil;
        double *Totabs_soil;
        double *LWabs_soil;
        double *LWemit_soil;
        double *Fc_soil;
        double *LE_soil;
        double *H_soil;
        double *G;
        double *E_soil;
        double *T_surf;
        double *qss;
        double *PH;
        double *htb;
        double *diff;
        int *tbc;
        int *cases;
        double *zhs;
        double *dzs;
        double *zns;
        double *rootfr;
        double *sand;
        double *clay;
        double *HC_sol;
        double *porosity;
        double *psi0;
        double *bsw;
        double *TK_sol;
        double *TK_dry;
        double *HKsat;
        double *theta_dry;
        double *volliq;
        double *smp;
        double *dsmp;
        double *dwat;
        double *Ts;
        double *krad;
        double *kax;
        double *rpp;
        VerticalSoilClass (int size)
        {
            rpp_wgt        = new double[1];
            ppt_ground     = new double[1];
            pptrate_ground = new double[1];
            SWabs_soil     = new double[1];
            Rnrad_soil     = new double[1];
            Totabs_soil    = new double[1];
            LWabs_soil     = new double[1];
            LWemit_soil    = new double[1];
            Fc_soil        = new double[1];
            LE_soil        = new double[1];
            H_soil         = new double[1];
            T_surf         = new double[1];
            G              = new double[1];
            E_soil         = new double[1];
            qss            = new double[1];
            PH             = new double[1];
            htb            = new double[1];
            diff           = new double[1];
            tbc            = new int[1];
            cases          = new int[1];
            zhs            = new double[size];
            dzs            = new double[size];
            zns            = new double[size];
            rootfr         = new double[size];
            sand           = new double[size];
            clay           = new double[size];
            HC_sol         = new double[size];
            porosity       = new double[size];
            psi0           = new double[size];
            bsw            = new double[size];
            TK_sol         = new double[size];
            TK_dry         = new double[size];
            HKsat          = new double[size];
            theta_dry      = new double[size];
            volliq         = new double[size];
            smp            = new double[size];
            dsmp           = new double[size];
            dwat           = new double[size];
            Ts             = new double[size];
            krad           = new double[size];
            kax            = new double[size];
            rpp            = new double[size];
        };
    };
#endif


#ifndef EIGENCANOPYCLASS_H
#define EIGENCANOPYCLASS_H
    class EigenCanopyClass
    {
        public:
        Eigen::VectorXd LADnorm;
        Eigen::VectorXd LAIsun;
        Eigen::VectorXd LAIshade;

        EigenCanopyClass (int size)
        {
            LADnorm.resize(size);
            LAIsun.resize(size);
            LAIshade.resize(size);
        }
    };
#endif


#ifndef EIGENSOILCLASS_H
#define EIGENSOILCLASS_H
    class EigenSoilClass
    {
        public:
        Eigen::VectorXd Ts;

        EigenSoilClass (int size)
        {
            Ts.resize(size);
        }
    };
#endif