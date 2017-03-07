#ifndef __SUBSURFACE_H__
#define __SUBSURFACE_H__
#endif


__global__ void SetUpLinearSystemsSubsurface(
                double *a3d, double *rhs, double *psinp1m, double *psin, double *thetanp1m, 
                double *thetan, double *knp1m, double *cnp1m, int *wbc, int *ebc, int *sbc,
                int *nbc, int *tbc, int *bbc, double *psi_w, double *psi_e, double *psi_s,
                double *psi_n, double *psi_t, double *psi_b, double *qw, double *qe, double *qs,
                double *qn, double *qt, double *qb, int3 globsize);


__global__ void GetIterationDifference(double *psinp1m, double *psinp1mp1, double* deltam, 
                int3 globsize);


__global__ void ModifiedPicardUpdate(double *psinp1m, double *psinp1mp1, int3 globsize);


__global__ void SubsurfaceEstimateInfiltrationPonding(
                double *psinp1mp1, double *knp1m, double *qt, double *qss, double *psi_top, 
                double *ph, int *tbc, int3 globsize);

void GatherFluxesDomain(ProjectClass *project, VerticalCanopyClass *vertcanopies,
                        VerticalSoilClass *vertsoils, SubsurfaceFlowClass *subsurface_host, 
                        SubsurfaceFlowClass *subsurface_dev, int rank, int procsize, int3 globsize,
                        int3 domsize, int2 topolsize, int2 topolindex, MPI_Comm *cartComm);

void SubsurfaceFlowModel(TimeForcingClass * &timeforcings, OverlandFlowClass * &overland_host,
                         OverlandFlowClass * &overland_dev, SubsurfaceFlowClass * &subsurface_host,
                         SubsurfaceFlowClass * &subsurface_dev, cuspdev_diamat &a3d_cusp,
                         cuspdev_1d &psinp1mp1_cusp, cuspdev_1d &rhs3d_cusp, cuspdev_idoper &id3d,
                         cuspdev_1d &deltam_cusp, thrustdev_iter &maxError, 
                         thrustdev &quflux_thrust, thrustdev &qdflux_thrust, thrustdev &qwflux_thrust, thrustdev &qeflux_thrust, thrustdev &qsflux_thrust, thrustdev &qnflux_thrust,
                         thrustdev &dtheta_thrust, thrustdev &transp_thrust, thrustdev &evapo_thrust, thrustdev &ssflux_thrust,
                         int rank, int procsize,
                         int3 globsize, int t, int num_steps);