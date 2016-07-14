#ifndef __OVERLAND_H__
#define __OVERLAND_H__
#endif

__global__ void SetUpLinearSystemsOverland(double *a2d, double *rhs2d, double *waterelev, double *ztopo, double *kwest, double *keast, double *ksouth, double *knorth, double ppt, double et, double infil, int3 globsize);


void OverlandFlowModel(TimeForcingClass * &timeforcings, OverlandFlowClass * &overland_dev,
                       SubsurfaceFlowClass * &subsurface_dev, cuspdev_diamat &a2d_cusp,
                       cuspdev_1d &we_out_cusp, cuspdev_1d &rhs2d_cusp, cuspdev_idoper &id2d,
                       int rank, int procsize, int3 globsize, int t, int num_steps);