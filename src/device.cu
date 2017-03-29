/*
// Copyright (C) 2016, HydroComplexity Group
// All rights reserved.
//
// Distributed Hydrologicc and Regional Analysis (DHARA) Model
// DHARA model is made available as a restricted, non-exclusive, 
// non-transferable license for education and research purpose only, 
// and not for commercial use. See the LICENSE.txt for more details.
//
// Author: levuvietphong@gmail.com (Phong Le)
*/


#include "../include/main.h"
#include "../include/cusplib.h"
#include "../include/mlcan.h"
#include "../include/global.h"
#include "../include/overland.h"
#include "../include/subsurface.h"
#include "../include/devconst.h"


// ========================
// Cusp and Thrust wrappers
// ========================

void ResizeCuspMemory(int rank, int3 globsize, cuspdev_diamat &a2d_cusp, cuspdev_1d &we_out_cusp,
                      cuspdev_1d &rhs2d_cusp, cuspdev_idoper &id2d, cuspdev_diamat &a3d_cusp,
                      cuspdev_1d &psinp1mp1_cusp, cuspdev_1d &rhs3d_cusp, cuspdev_idoper &id3d,
                      cuspdev_1d &deltam_cusp, thrustdev &thetanp1mp1_thrust,
                      thrustdev &quflux_thrust, thrustdev &qdflux_thrust, thrustdev &qwflux_thrust, thrustdev &qeflux_thrust, thrustdev &qsflux_thrust, thrustdev &qnflux_thrust,
                      thrustdev &dtheta_thrust, thrustdev &transp_thrust, thrustdev &evapo_thrust, thrustdev &ssflux_thrust)
{
    int globsize_xy  = globsize.x * globsize.y;
    int globsize_xyz = globsize.x * globsize.y * globsize.z;

    if (rank == MPI_MASTER_RANK)
    {
        printf("\nRESIZING MEMORY . . . .");


        // 2d overland flow
        cusp::gallery::poisson5pt(a2d_cusp, globsize.x, globsize.y);
        rhs2d_cusp.resize(globsize_xy, 1.0);
        we_out_cusp.resize(globsize_xy, 0.0);
        id2d.resize(a2d_cusp.num_rows, a2d_cusp.num_rows, 0.0);

        // 3d subsurface flow
        cusp::gallery::poisson7pt(a3d_cusp, globsize.x, globsize.y, globsize.z);
        rhs3d_cusp.resize(globsize_xyz, 2.0);
        psinp1mp1_cusp.resize(globsize_xyz, 0.0);
        thetanp1mp1_thrust.resize(globsize_xyz, 0.0);
        deltam_cusp.resize(globsize_xyz, 0.0);
        id3d.resize(a3d_cusp.num_rows, a3d_cusp.num_rows, 0.0);
        
        quflux_thrust.resize(globsize_xyz, 0.0);
        qdflux_thrust.resize(globsize_xyz, 0.0);
        qwflux_thrust.resize(globsize_xyz, 0.0);
        qeflux_thrust.resize(globsize_xyz, 0.0);
        qsflux_thrust.resize(globsize_xyz, 0.0);
        qnflux_thrust.resize(globsize_xyz, 0.0);
        dtheta_thrust.resize(globsize_xyz, 0.0);
        transp_thrust.resize(globsize_xyz, 0.0);
        evapo_thrust.resize(globsize_xyz, 0.0);
        ssflux_thrust.resize(globsize_xyz, 0.0);        


        printf(" . . . . . . . . . . . . . . . . . . . . completed! \n");
        printf("--------------- \n");

    }
}


/**
 * @brief      Map cusp memory to plain memory for passing to CUDA kernels
 *
 * @param      forcings        The forcings
 * @param      overland_dev    The overland dev
 * @param      subsurface_dev  The subsurface dev
 * @param[in]  rank            The rank
 * @param[in]  globsize        The globsize
 * @param[in]  a2d_cusp        The 2D matrix in cusp
 * @param      we_out_cusp     Water elevation in cusp
 * @param      rhs2d_cusp      Right hand side of 2D cusp system
 * @param[in]  a3d_cusp        The 3D matrix in cusp
 * @param      psinp1mp1_cusp  Pressure head at n+1,m+1 in cusp
 * @param      rhs3d_cusp      Right hand side of 3D cusp system
 * @param      deltam_cusp     The deltam cusp
 */
void CastCuspToPlainMemory(TimeForcingClass * &timeforcings, OverlandFlowClass * &overland_dev,
                           SubsurfaceFlowClass * &subsurface_dev, int rank, int3 globsize,
                           cuspdev_diamat &a2d_cusp, cuspdev_1d &we_out_cusp,
                           cuspdev_1d &rhs2d_cusp, cuspdev_diamat &a3d_cusp,
                           cuspdev_1d &psinp1mp1_cusp, cuspdev_1d &rhs3d_cusp,
                           cuspdev_1d &deltam_cusp, thrustdev &thetanp1mp1_thrust,
                           thrustdev &quflux_thrust, thrustdev &qdflux_thrust, thrustdev &qwflux_thrust, thrustdev &qeflux_thrust, thrustdev &qsflux_thrust, thrustdev &qnflux_thrust,
                           thrustdev &dtheta_thrust, thrustdev &transp_thrust, thrustdev &evapo_thrust, thrustdev &ssflux_thrust)
{
    if (rank == MPI_MASTER_RANK)
    {
        printf("\nCASTING CUSP MEMORY . . . .");

        // 2D overland flow
        overland_dev->a2d           = thrust::raw_pointer_cast(&a2d_cusp.values(0,0));
        overland_dev->waterelev     = thrust::raw_pointer_cast(&we_out_cusp[0]);
        overland_dev->rhs2d         = thrust::raw_pointer_cast(&rhs2d_cusp[0]);

        // 3D subsurface flow
        subsurface_dev->a3d         = thrust::raw_pointer_cast(&a3d_cusp.values(0,0));
        subsurface_dev->psinp1mp1   = thrust::raw_pointer_cast(&psinp1mp1_cusp[0]);
        subsurface_dev->rhs3d       = thrust::raw_pointer_cast(&rhs3d_cusp[0]);
        subsurface_dev->deltam      = thrust::raw_pointer_cast(&deltam_cusp[0]);
        subsurface_dev->thetanp1mp1 = thrust::raw_pointer_cast(&thetanp1mp1_thrust[0]);

        subsurface_dev->quflux = thrust::raw_pointer_cast(&quflux_thrust[0]);
        subsurface_dev->qdflux = thrust::raw_pointer_cast(&qdflux_thrust[0]);
        subsurface_dev->qwflux = thrust::raw_pointer_cast(&qwflux_thrust[0]);
        subsurface_dev->qeflux = thrust::raw_pointer_cast(&qeflux_thrust[0]);
        subsurface_dev->qsflux = thrust::raw_pointer_cast(&qsflux_thrust[0]);
        subsurface_dev->qnflux = thrust::raw_pointer_cast(&qnflux_thrust[0]);
        subsurface_dev->dtheta = thrust::raw_pointer_cast(&dtheta_thrust[0]);
        subsurface_dev->transp = thrust::raw_pointer_cast(&transp_thrust[0]);
        subsurface_dev->evapo = thrust::raw_pointer_cast(&evapo_thrust[0]);
        subsurface_dev->ssflux = thrust::raw_pointer_cast(&ssflux_thrust[0]);
        

        printf(" . . . . . . . . . . . . . . . . . . completed! \n");
        printf("-------------------- \n");
    }
}



void CopyLayers(SubsurfaceFlowClass * &subsurface_host, cuspdev_1d psinp1mp1_cusp,
                thrustdev &thetanp1mp1_thrust, int t, int sizexy, int sizez)
{
    for (int k=0; k<sizez; k++)
    {
        subsurface_host->psi_col[t*sizez+k]   = thrust::reduce(
                                                psinp1mp1_cusp.begin() + k*sizexy,
                                                psinp1mp1_cusp.begin() + (k+1)*sizexy-1 );
        subsurface_host->theta_col[t*sizez+k] = thrust::reduce(
                                                thetanp1mp1_thrust.begin() + k*sizexy,
                                                thetanp1mp1_thrust.begin() + (k+1)*sizexy-1 );
    }
}


// =============
// Host wrappers
// =============

/**
 * @brief The host function for checking the result of a CUDA API call
 *
 * @param[in]   command         The result of the previously-issued CUDA API call
 * @param[in]   commandName     The name of the issued API call
 * @param[in]   fileName        The name of the file where the API call occurred
 * @param[in]   line            The line in the file where the API call occurred
 */
extern "C" \
void CheckCudaCall(cudaError_t command, const char * commandName, const char * fileName, int line)
{
    if (command != cudaSuccess)
    {
        fprintf(stderr, "Error: CUDA result \"%s\" for call \"%s\" in file \"%s\" at line %d. Terminating...\n",
            cudaGetErrorString(command), commandName, fileName, line);
        exit(STATUS_ERR);
    }
}


// This function is temporarily unused to due no support in cudaMemcpyToSymbol for linking
extern "C" \
void CopyConstantToDevice(ProjectClass *project, OverlandFlowClass * &overland,
                          SubsurfaceFlowClass *subsurface)
{
    ///////////////////////////////////////
    // Constants                         //
    ///////////////////////////////////////

    SafeCudaCall(cudaMemcpyToSymbol(dt, &project->dtimehr, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(dx, &project->dx_meter, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(dy, &project->dy_meter, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(dz, &project->dz_meter, sizeof(double)));

    SafeCudaCall(cudaMemcpyToSymbol(alpha, &subsurface->alpha, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(poros, &subsurface->porosity, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(theta_S, &subsurface->theta_s, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(theta_R, &subsurface->theta_r, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(Ss, &subsurface->Ss, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(n, &subsurface->poresize, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(air_dry, &subsurface->airdry, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(psimin, &subsurface->psimin, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(am, &subsurface->am, sizeof(double)));

    SafeCudaCall(cudaMemcpyToSymbol(delta, &overland->delta, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(hmin, &overland->hmin, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(hcri, &overland->hcri, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(K0, &overland->K0, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(hn, &overland->hn, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(hs, &overland->hs, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(hw, &overland->hw, sizeof(double)));
    SafeCudaCall(cudaMemcpyToSymbol(he, &overland->he, sizeof(double)));
}




/**
 * @brief      Copy data from host to device memory.
 *
 * @param      project          Class including project info
 * @param      timeforcings     forcings data
 * @param      overland_host    overland data in host memory
 * @param      overland_dev     overland data in dev memory
 * @param      subsurface_host  subsurface data in host memory
 * @param      subsurface_dev   subsurface data in dev memory
 * @param[in]  globsize         size of global domain
 * @param[in]  num_steps  number of time steps for simulation
 */
extern "C" \
void CopyDataToDevice(ProjectClass *project, TimeForcingClass *timeforcings,
                      OverlandFlowClass *overland_host, OverlandFlowClass *overland_dev,
                      SubsurfaceFlowClass *subsurface_host, SubsurfaceFlowClass *subsurface_dev,
                      int3 globsize)
{
    int sizexy  = globsize.x * globsize.y;
    int sizexz  = globsize.x * globsize.z;
    int sizeyz  = globsize.y * globsize.z;
    int sizexyz = globsize.x * globsize.y * globsize.z;

    printf("\nCOPYING DATA TO DEVICE");

    ///////////////////////////////////////
    // 2D overland flow                  //
    ///////////////////////////////////////

    SafeCudaCall( cudaMemcpy(overland_dev->mann, overland_host->mann,
            sizexy*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(overland_dev->waterelev, overland_host->waterelev,
            sizexy*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(overland_dev->waterdepth, overland_host->waterdepth,
            sizexy*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(overland_dev->ztopo, overland_host->ztopo,
            sizexy*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(overland_dev->ph, overland_host->ph,
            sizexy*sizeof(double), cudaMemcpyHostToDevice) );

    ///////////////////////////////////////
    // 3D subsurface flow                //
    ///////////////////////////////////////

    // Vegetation mapping
    SafeCudaCall( cudaMemcpy(subsurface_dev->procmap, subsurface_host->procmap,
            sizexy*sizeof(int), cudaMemcpyHostToDevice) );

    // Boundary types and conditions
    SafeCudaCall( cudaMemcpy(subsurface_dev->bcb, subsurface_host->bcb,
            sizexy*sizeof(int), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bct, subsurface_host->bct,
            sizexy*sizeof(int), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bce, subsurface_host->bce,
            sizeyz*sizeof(int), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcw, subsurface_host->bcw,
            sizeyz*sizeof(int), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcn, subsurface_host->bcn,
            sizexz*sizeof(int), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcs, subsurface_host->bcs,
            sizexz*sizeof(int), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcpsib, subsurface_host->bcpsib,
            sizexy*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcpsit, subsurface_host->bcpsit,
            sizexy*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcpsie, subsurface_host->bcpsie,
            sizeyz*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcpsiw, subsurface_host->bcpsiw,
            sizeyz*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcpsin, subsurface_host->bcpsin,
            sizexz*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcpsis, subsurface_host->bcpsis,
            sizexz*sizeof(double), cudaMemcpyHostToDevice) );


    SafeCudaCall( cudaMemcpy(subsurface_dev->bcqb, subsurface_host->bcqb,
            sizexy*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcqt, subsurface_host->bcqt,
            sizexy*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcqe, subsurface_host->bcqe,
            sizeyz*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcqw, subsurface_host->bcqw,
            sizeyz*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcqn, subsurface_host->bcqn,
            sizexz*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->bcqs, subsurface_host->bcqs,
            sizexz*sizeof(double), cudaMemcpyHostToDevice) );


    // Main 3d subsurface variables
    SafeCudaCall( cudaMemcpy(subsurface_dev->ksat, subsurface_host->ksat,
            sizexyz*sizeof(double), cudaMemcpyHostToDevice) );

    SafeCudaCall( cudaMemcpy(subsurface_dev->thetan, subsurface_host->thetan,
            sizexyz*sizeof(double), cudaMemcpyHostToDevice) );

    vanGenuchtenInverse<<<TSZ,BSZ>>>(subsurface_dev->thetan, subsurface_dev->psin, sizexyz);
    cudaCheckError("vanGenuchtenInverse");

    SafeCudaCall( cudaMemcpy(subsurface_dev->psinp1m, subsurface_dev->psin,
            sizexyz*sizeof(double), cudaMemcpyDeviceToDevice) );

    // If reach this point, print out info
    printf(" . . . . . . . . . . . . . . . . . . . . .completed! \n");
    printf("------------------------- \n");

}


/**
 * @brief      Run the model
 *
 * @param      timeforcings     The timeforcings
 * @param      overland_host    The overland host
 * @param      overland_dev     The overland dev
 * @param      subsurface_host  The subsurface host
 * @param      subsurface_dev   The subsurface dev
 * @param      files            The files
 * @param      project          The project
 * @param      forcings         The forcings
 * @param      switches         The switches
 * @param      constants        The constants
 * @param      canopies         The canopies
 * @param      soils            The soils
 * @param      radiation        The radiation
 * @param      photosynthesis   The photosynthesis
 * @param      respiration      The respiration
 * @param      stomaconduct     The stomaconduct
 * @param      microenviron     The microenviron
 * @param      vertcanopies     The vertcanopies
 * @param      vertsoils        The vertsoils
 * @param      eigencanopies    The eigencanopies
 * @param      eigensoils       The eigensoils
 * @param      outmlcan         The outmlcan
 * @param[in]  rank             The rank
 * @param[in]  procsize         The procsize
 * @param[in]  globsize         The globsize
 * @param[in]  domsize          The domsize
 * @param[in]  topolsize        The topolsize
 * @param[in]  topolindex       The topolindex
 * @param      cartComm         The cartesian communications
 * @param[in]  num_steps  The num steps
 */
extern "C" \
void RunCoupledFlowModel(TimeForcingClass *timeforcings, OverlandFlowClass *overland_host,
                         OverlandFlowClass *overland_dev, SubsurfaceFlowClass *subsurface_host,
                         SubsurfaceFlowClass *subsurface_dev, FileNameClass *files,
                         ProjectClass *project, ForcingClass * &forcings, SwitchClass *switches,
                         ConstantClass *constants, CanopyClass *canopies, SoilClass *soils,
                         RadiationClass *radiation, PhotosynthesisClass *photosynthesis,
                         RespirationClass *respiration, StomaConductClass *stomaconduct,
                         MicroEnvironmentClass *microenviron, VerticalCanopyClass *vertcanopies,
                         VerticalSoilClass *vertsoils, EigenCanopyClass *eigencanopies,
                         EigenSoilClass *eigensoils, OutputClass *outmlcan, int rank, int procsize,
                         int3 globsize, int3 domsize, int2 topolsize, int2 topolindex,
                         MPI_Comm *cartComm)
{
    int sizex = globsize.x;
    int sizey = globsize.y;
    int sizez = globsize.z;
    int sizexy = sizex * sizey;
    int time_print = project->printiterval;
    int num_steps = project->num_steps;
    int isroot = rank == MPI_MASTER_RANK;

    /* CUSP variables - 2D and 3D flow model.*/
    cuspdev_diamat a2d_cusp, a3d_cusp;                      // diagonal matrices
    cuspdev_idoper id2d, id3d;                              // identity operator
    thrustdev_iter maxError;                                // iterator error
    cuspdev_1d we_out_cusp, rhs2d_cusp;                     // 1d array - for 2d model
    cuspdev_1d psinp1mp1_cusp, rhs3d_cusp, deltam_cusp;     // 1d array - for 3d model
    thrustdev thetanp1mp1_thrust,
        quflux_thrust, qdflux_thrust, qwflux_thrust, qeflux_thrust, qsflux_thrust, qnflux_thrust,
        dtheta_thrust, transp_thrust, evapo_thrust, ssflux_thrust;  // 1d thrust array
                                                            // All: see definitions in cusplib.h
                                                            // cusp is based on thrust

    // Process CUSP memory.
    ResizeCuspMemory(rank, globsize, a2d_cusp, we_out_cusp, rhs2d_cusp, id2d, a3d_cusp,
                     psinp1mp1_cusp, rhs3d_cusp, id3d, deltam_cusp, thetanp1mp1_thrust,
                     quflux_thrust, qdflux_thrust, qwflux_thrust, qeflux_thrust, qsflux_thrust, qnflux_thrust,
                     dtheta_thrust, transp_thrust, evapo_thrust, ssflux_thrust);

    CastCuspToPlainMemory(timeforcings, overland_dev, subsurface_dev, rank, globsize, a2d_cusp,
                          we_out_cusp, rhs2d_cusp, a3d_cusp, psinp1mp1_cusp, rhs3d_cusp,
                          deltam_cusp, thetanp1mp1_thrust,
                          quflux_thrust, qdflux_thrust, qwflux_thrust, qeflux_thrust, qsflux_thrust, qnflux_thrust,
                          dtheta_thrust, transp_thrust, evapo_thrust, ssflux_thrust);

    // Initialization for flow model
    PreRunningFlowModel(project, subsurface_dev, rank, procsize, globsize, cartComm);

    // Time loop for simulation.
    double timerStart = MPI_Wtime();
    double elapsedTime;
    for (int t = 0; t < num_steps; t++)
    {
        elapsedTime = MPI_Wtime() - timerStart;
        OnePrintf(isroot && t%time_print == 0, "Simulations %5d of %5d completed. \
                  Elapsed time (s): %6.3f \n", t, num_steps, elapsedTime);

        // MLCan model on MPI
        CanopyModel(project, switches, constants, canopies, vertcanopies, soils, vertsoils,
                    timeforcings, forcings, radiation, photosynthesis, stomaconduct, respiration,
                    microenviron, outmlcan, t, rank, procsize);

        MPI_Barrier(*cartComm);              // Synchronize all MPIs

        

        // Gather transpiration to master
        GatherFluxesDomain(project, vertcanopies, vertsoils, subsurface_host, subsurface_dev, rank,
                           procsize, globsize, domsize, topolsize, topolindex, cartComm);

        MPI_Barrier(*cartComm);              // Synchronize all MPIs

        // Flow3D model in device. Run on root/master process only.
        if (isroot)
        {
            SubsurfaceFlowModel(timeforcings, overland_host, overland_dev, subsurface_host,
                                subsurface_dev, a3d_cusp, psinp1mp1_cusp, rhs3d_cusp, id3d,
                                deltam_cusp, maxError, 
                                quflux_thrust, qdflux_thrust, qwflux_thrust, qeflux_thrust, qsflux_thrust, qnflux_thrust,
                                dtheta_thrust, transp_thrust, evapo_thrust, ssflux_thrust,
                                rank, procsize, globsize, t, num_steps);

            OverlandFlowModel(timeforcings, overland_dev, subsurface_dev, a2d_cusp, we_out_cusp,
                              rhs2d_cusp, id2d, rank, procsize, globsize, t, num_steps);


            if (project->savestat)
            {            
                // Get mean value in each layers.
                CopyLayers(subsurface_host, psinp1mp1_cusp, thetanp1mp1_thrust, t, sizexy, sizez);
            }

            // Save 2D and 3D results.
            SaveModelResults(project, overland_host, overland_dev, subsurface_host, subsurface_dev,             globsize, t);
        }
    }

    SaveResultEntirePeriod(project, canopies, subsurface_host, outmlcan, rank ,procsize, 
                           globsize, cartComm);

}
