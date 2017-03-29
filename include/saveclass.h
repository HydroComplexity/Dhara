#ifndef OUTPUTCLASS_H
#define OUTPUTCLASS_H
    class OutputClass
    {
        public:
        double *An_can;
        double *LE_can;
        double *H_can;
        double *TR_can;
        double *Rnrad_can;

        // Soil total
        double *E_soil;
        double *qss;
        double *htb;
        double *diff;
        double *PH;
        double *G;
        double *T_surf;
        double *ppt_ground;

        // Mass balance
        double *mbw_can;

        // Canopy Profile
        double *An_sun;
        double *An_shade;
        double *LE_sun;
        double *LE_shade;
        double *H_sun;
        double *H_shade;
        double *Evap_prof;
        double *TR_sun;
        double *TR_shade;
        double *LAI_sun;
        double *LAI_shade;
        double *zhc;

        // Soil Profile
        double *volliq;
        double *smp;
        double *krad;
        double *rpp;
        double *Ts;
        double *zhs;
    };
#endif