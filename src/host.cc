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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "../include/main.h"
#include "../include/parser.h"


/**
 * @brief      Set the initial conditions in host memory
 *
 * @param      forcings    Class including forcings info
 * @param      overland    Overland flow class
 * @param      subsurface  Subsurface flow class
 * @param[in]  num_steps   The number of timesteps
 * @param[in]  rank        Global rank of the current MPI process
 * @param[in]  procsize    Total number of MPI processes available
 * @param[in]  globsize    Size of the global domain
 * @param      cartComm    Carthesian MPI communicator
 */
void SetInitialConditionsHost(TimeForcingClass *timeforcings, FileNameClass *files,
                              OverlandFlowClass *overland, SubsurfaceFlowClass * &subsurface,
                              int num_steps, int rank, int procsize, int3 globsize)
{
    // Get sizes of domain
    int sizez   = globsize.z;
    int sizexy  = globsize.x * globsize.y;    

    printf("\nSETTING UP INITIAL CONDITIONS IN HOST \n");
    printf("--------------------------------------- \n");

    // Overland flow
    for (int i = 0; i < sizexy; i++)
    {
        overland->mann[i] = 0.025;
        overland->waterdepth[i] = 0.;
        overland->waterelev[i]  = overland->ztopo[i] + overland->waterdepth[i];
        overland->ph[i] = 0.;
    }

    // Soil moisture and hydraulic conductivity
    double ksat = 0.0025;//5e-3; // [m/dtime]
    for (int k = 0; k < sizez; k++)
    {
        for (int i = 0; i < sizexy; i++)
        {
            subsurface->ksat[k*sizexy+i] = ksat;
            //subsurface->ksat[k*sizexy+i] = ksat * std::exp(-(k*0.1));
            subsurface->thetan[k*sizexy+i] = 0.3 + 0.01*k;
        }
    }

    printf("\t Flow model - Initialization . . . . . . . . . . . . . .completed! \n");
}



/**
 * @brief      Set the boundary conditions host in host memory
 *
 * @param      forcings    Class including forcings info
 * @param      overland    Overland flow class
 * @param      subsurface  Subsurface flow class
 * @param[in]  num_steps   The number of timesteps
 * @param[in]  rank        Global rank of the current MPI process
 * @param[in]  procsize    Total number of MPI processes available
 * @param[in]  globsize    Size of the global domain
 * @param      cartComm    Carthesian MPI communicator
 */
void SetBoundaryConditionsHost(TimeForcingClass *timeforcings, FileNameClass *files,
                               OverlandFlowClass *overland, SubsurfaceFlowClass * &subsurface,
                               int num_steps, int rank, int procsize, int3 globsize)
{
    // Get sizes of domain
    int sizexy = globsize.x * globsize.y;
    int sizexz = globsize.x * globsize.z;
    int sizeyz = globsize.y * globsize.z;
    int bct, bcb, bce, bcw, bcn, bcs;
    double bcpsit, bcpsib, bcpsie, bcpsiw, bcpsin, bcpsis;
    double bcqt, bcqb, bcqe, bcqw, bcqn, bcqs;

    bct    = GetOptionToInt("BOUNDARY_TYPE_TOP");
    bcb    = GetOptionToInt("BOUNDARY_TYPE_BOTTOM");
    bce    = GetOptionToInt("BOUNDARY_TYPE_EAST");
    bcw    = GetOptionToInt("BOUNDARY_TYPE_WEST");
    bcn    = GetOptionToInt("BOUNDARY_TYPE_NORTH");
    bcs    = GetOptionToInt("BOUNDARY_TYPE_SOUTH");

    bcpsit = GetOptionToDouble("PRESSURE_HEAD_BOTTOM");
    bcpsib = GetOptionToDouble("PRESSURE_HEAD_TOP");
    bcpsie = GetOptionToDouble("PRESSURE_HEAD_EAST");
    bcpsiw = GetOptionToDouble("PRESSURE_HEAD_WEST");
    bcpsin = GetOptionToDouble("PRESSURE_HEAD_NORTH");
    bcpsis = GetOptionToDouble("PRESSURE_HEAD_SOUTH");

    bcqt   = GetOptionToDouble("FLUX_BOTTOM");
    bcqb   = GetOptionToDouble("FLUX_TOP");
    bcqe   = GetOptionToDouble("FLUX_EAST");
    bcqw   = GetOptionToDouble("FLUX_WEST");
    bcqn   = GetOptionToDouble("FLUX_NORTH");
    bcqs   = GetOptionToDouble("FLUX_SOUTH");

    ParserConfigFile(files->config);

    printf("\nSETTING UP BOUNDARY CONDITIONS IN HOST");

    ///////////////////////////////////////////
    // Overland flow model                   //
    ///////////////////////////////////////////

    // TODO: Add boundary conditions for overland flow model


    ///////////////////////////////////////////
    // Subsurface flow model                 //
    ///////////////////////////////////////////

    // Top and bottom boundary conditions
    for (int i = 0; i < sizexy; i++)
    {
        subsurface->bct[i]    = bct;
        subsurface->bcb[i]    = bcb;
        subsurface->bcpsit[i] = bcpsit;
        subsurface->bcpsib[i] = bcpsib;
        subsurface->bcqt[i]   = bcqt;
        subsurface->bcqb[i]   = bcqb;
    }

    // East and west boundary conditions
    for (int i = 0; i < sizeyz; i++)
    {
        subsurface->bce[i]    = bce;
        subsurface->bcw[i]    = bcw;
        subsurface->bcpsie[i] = bcpsie;
        subsurface->bcpsiw[i] = bcpsiw;
        subsurface->bcqe[i]   = bcqe;
        subsurface->bcqw[i]   = bcqw;
    }

    // North and south boundary conditions
    for (int i = 0; i < sizexz; i++)
    {
        subsurface->bcn[i]    = bcn;
        subsurface->bcs[i]    = bcs;
        subsurface->bcpsin[i] = bcpsin;
        subsurface->bcpsis[i] = bcpsis;
        subsurface->bcqn[i]   = bcqn;
        subsurface->bcqs[i]   = bcqs;
    }

    printf(" . . . . . . . . . . . . .completed! \n");
    printf("-------------------------------------- \n");
}



/**
 * @brief      Set conditions in host memory
 *
 * @param      forcings    Class including forcings info
 * @param      overland    Overland flow class
 * @param      subsurface  Subsurface flow class
 * @param[in]  num_steps   The number of timesteps
 * @param[in]  rank        Global rank of the current MPI process
 * @param[in]  procsize    Total number of MPI processes available
 * @param[in]  globsize    Size of the global domain
 * @param      cartComm    Carthesian MPI communicator
 */
void SetFlowModelConditions(TimeForcingClass *timeforcings, FileNameClass *files,
                            OverlandFlowClass *overland, SubsurfaceFlowClass * &subsurface,
                            int num_steps, int rank, int procsize, int3 globsize)
{
    SetInitialConditionsHost(timeforcings, files, overland, subsurface, num_steps, rank, procsize,
                             globsize);
    SetBoundaryConditionsHost(timeforcings, files, overland, subsurface, num_steps, rank, procsize,
                              globsize);
}



/**
 * @brief      Sets the MPI map used for flow model in GPU
 *
 * @param      project     Class including project info
 * @param      subsurface  Subsurface flow class
 * @param[in]  procsize    Total number of MPI processes available
 * @param[in]  globsize    Size of the global domain
 * @param[in]  topolsize   Size of topology nodes (mpi)
 */
void SetMPIGPUMapping(mpiClass *mpiobj, SubsurfaceFlowClass * &subsurface, int procsize,
                      int3 globsize, int2 topolsize)
{
    int imax, jmax, globind, proc;

    proc = 0;
    for (int jo = 0; jo < topolsize.y; jo++)
    {
        // use jmax to find upper bound of the loop
        jmax = jo < topolsize.y-1 ? mpiobj->ydispls[jo+1] : globsize.y;
        for (int io = 0; io < topolsize.x; io++)
        {
            // use imax to find upper bound of the loop as well
            imax = io < topolsize.x-1 ? mpiobj->xdispls[io+1] : globsize.x;
            for (int j = mpiobj->ydispls[jo]; j < jmax; j++)
            {
                for (int i = mpiobj->xdispls[io]; i < imax; ++i)
                {
                    globind = j * globsize.x + i;           // global index in the entire domain
                    subsurface->procmap[globind] = proc;      // assign process id to the map
                }
            }
            proc = proc + 1;
        }
    }
}


/**
 * @brief      Collect root density profile from all processes to master
 *
 * @param      project          Class including project info
 * @param      vertsoils        Class including vertical soil info
 * @param      soils            Class including soil variables
 * @param      subsurface_host  Subsurface flow class in host memory
 * @param      subsurface_dev   Subsurface flow class in device memory
 * @param[in]  rank             Global rank of the current MPI process
 * @param[in]  procsize         Size of the global domain
 * @param      cartComm         Carthesian MPI communicator
 */
void GatherRootDensityToMasterAndFilter(ProjectClass *project, VerticalSoilClass *vertsoils,
                                        SoilClass *soils, SubsurfaceFlowClass *subsurface_host,
                                        SubsurfaceFlowClass *subsurface_dev, int rank, 
                                        int procsize, MPI_Comm *cartComm)
{
    int nl_soil = soils->nl_soil;
    int isroot = rank == MPI_MASTER_RANK;

    // Gather root density to master process
    MPI_Gather(vertsoils->rootfr, nl_soil, MPI_DOUBLE, subsurface_host->rda, nl_soil, MPI_DOUBLE,
               MPI_MASTER_RANK, *cartComm);
    MPI_Gather(&project->planttype, 1, MPI_INT, subsurface_host->type, 1, MPI_INT, 0, *cartComm);

    if (isroot)
    {
        SafeCudaCall( cudaMemcpy(subsurface_dev->rda, subsurface_host->rda,
                                 procsize*nl_soil*sizeof(double), cudaMemcpyHostToDevice) );
    }
}
