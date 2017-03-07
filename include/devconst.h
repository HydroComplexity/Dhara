/////////////////////////////////////////////////////////////////////
// Flow model parameters                                           //
/////////////////////////////////////////////////////////////////////

// Universal parameters
__constant__ double dx = 1.2;   		   // space resolution x-direction [m]
__constant__ double dy = 1.2;   		   // space resolution y-direction [m]
__constant__ double dz = 0.2;   		   // space resolution z-direction [m]
__constant__ double dt = 0.5;   		   // time step [hour]
__constant__ double sec_p_mm2dt_p_m = 1.8; // sec/mm to dt/m, unit convertion

// Subsurface model
__constant__ double alpha   = 0.02;
__constant__ double poros   = 0.45;
__constant__ double theta_S = 0.45;
__constant__ double theta_R = 0.1;
__constant__ double Ss      = 5e-4;
__constant__ double n       = 1.8;
__constant__ double air_dry = -5.0;
__constant__ double psimin  = 0.005;
__constant__ double am      = 0.0;       // Jumping in the soil moisture dynamics for the next iteration. It shuold be "0" not to cause additional unstability


// Overland flow model
__constant__ double delta = 1e-7;
__constant__ double hmin  = 1e-5;
__constant__ double hcri  = 0.0;
__constant__ double K0    = 1e-5;
__constant__ double hn    = 0.0;
__constant__ double hs    = 0.0;
__constant__ double hw    = 0.0;
__constant__ double he    = 0.0;