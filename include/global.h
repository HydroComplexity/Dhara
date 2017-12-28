#ifndef __GLOBAL_H__
#define __GLOBAL_H__
#endif


__global__ void vanGenuchtenInverse(double *theta, double *psi, int size);

__global__ void vanGenuchten(double *C, double *theta, double *Ksat, double *K, double *psi, 
							 int3 globsize);

__global__ void WaterElevationOverland(double *waterelev, double *waterdepth, double *ztopo, 
									   int size);

__global__ void EstimateFluxes(double *ph, double *hpoten, double *qcapa, double *psinp1m, 
							   double *knp1m, double *ppt, double *et, double *ksat, int3 globsize);

__global__ void IdentifyTopBoundary(double *hpoten, double *qcapa, int *topbc, double *topqflux, 
									double *psinp1m, double *Psi_top, 
									double *thetan, double *ksat, int3 globsize);


