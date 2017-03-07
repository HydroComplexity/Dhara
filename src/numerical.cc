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

#include <Eigen/Dense>
#include <iostream>
#include "../include/main.h"
#include "../include/vertclass.h"
#include "../include/parser.h"
//#include "../include/devconst.h"


///////////////////////////////////////////////////////
// Allocate and free memory used                     //
///////////////////////////////////////////////////////


/**
 * @brief      Allocate run-time memory for Flow3D model in host and device memory
 *
 * @param      timeforcings     Class including time forcings
 * @param      overland_host    Overland flow class in host memory
 * @param      overland_dev     Overland flow class in device memory
 * @param      subsurface_host  Subsurface flow class in host memory
 * @param      subsurface_dev   Subsurface flow class in device memory
 * @param[in]  procsize         Total number of MPI processes available
 * @param[in]  globsize         Size of the global domain
 * @param[in]  num_steps        Number of time steps for simulation
 */
void AllocateMemoryFlowModel(TimeForcingClass *timeforcings,
                             OverlandFlowClass * &overland_host,
                             OverlandFlowClass * &overland_dev,
                             SubsurfaceFlowClass * &subsurface_host,
                             SubsurfaceFlowClass * &subsurface_dev,
                             int procsize, int3 globsize, int num_steps)
{
    int sizez   = globsize.z;
    int sizexy  = globsize.x * globsize.y;
    int sizexz  = globsize.x * globsize.z;
    int sizeyz  = globsize.y * globsize.z;
    int sizexyz = globsize.x * globsize.y * globsize.z;

    printf("\nALLOCATING MEMORY \n");
    printf("------------------- \n");

    ///////////////////////////////////////////////
    // HOST MEMORY                               //
    ///////////////////////////////////////////////

    /* 2D spatial variables. */
    overland_host->a2d          = new double[5*sizexy];
    overland_host->hpoten       = new double[sizexy];
    overland_host->mann         = new double[sizexy];
    overland_host->ph           = new double[sizexy];
    overland_host->qcapa        = new double[sizexy];
    overland_host->rhs2d        = new double[sizexy];
    overland_host->ztopo        = new double[sizexy];
    overland_host->waterdepth   = new double[sizexy];
    overland_host->waterelev    = new double[sizexy];

    subsurface_host->bcpsib     = new double[sizexy];
    subsurface_host->bcpsit     = new double[sizexy];
    subsurface_host->bcqb       = new double[sizexy];
    subsurface_host->bcqt       = new double[sizexy];
    subsurface_host->bcpsie     = new double[sizeyz];
    subsurface_host->bcpsiw     = new double[sizeyz];
    subsurface_host->bcqe       = new double[sizeyz];
    subsurface_host->bcqw       = new double[sizeyz];
    subsurface_host->bcpsin     = new double[sizexz];
    subsurface_host->bcpsis     = new double[sizexz];
    subsurface_host->bcqn       = new double[sizexz];
    subsurface_host->bcqs       = new double[sizexz];
    subsurface_host->qss        = new double[sizexy];
    subsurface_host->ppt_ground = new double[sizexy];
    subsurface_host->TR         = new double[sizexy];

    subsurface_host->bcb        = new int[sizexy];
    subsurface_host->bct        = new int[sizexy];
    subsurface_host->bce        = new int[sizeyz];
    subsurface_host->bcw        = new int[sizeyz];
    subsurface_host->bcn        = new int[sizexz];
    subsurface_host->bcs        = new int[sizexz];
    subsurface_host->procmap    = new int[sizexy];
    subsurface_host->type       = new int[procsize];

    /* 2D temporal variables. */
    subsurface_host->psi_col    = new double[num_steps*sizez];
    subsurface_host->theta_col  = new double[num_steps*sizez];
    subsurface_host->rda        = new double[procsize*sizez];

    /* 3D spatial variables. */
    subsurface_host->a3d        = new double[7*sizexyz];
    subsurface_host->cnp1m      = new double[sizexyz];
    subsurface_host->knp1m      = new double[sizexyz];
    subsurface_host->ksat       = new double[sizexyz];
    subsurface_host->psin       = new double[sizexyz];
    subsurface_host->psinp1m    = new double[sizexyz];
    subsurface_host->psinp1mp1  = new double[sizexyz];
    subsurface_host->psin       = new double[sizexyz];
    subsurface_host->rhs3d      = new double[sizexyz];
    subsurface_host->thetan     = new double[sizexyz];
    subsurface_host->psiout     = new double[sizexyz];
    subsurface_host->thetaout   = new double[sizexyz];
    subsurface_host->TR_root    = new double[procsize];
    subsurface_host->ppt_root   = new double[procsize];
    subsurface_host->E_soil_root= new double[procsize];

    /* 1D temporal variable. */
    subsurface_host->mb_subsurfaceW = new double[num_steps];

    /* If reach this point, print out info. */
    printf("\t Flow model - host memory . . . . . . . . . . . . . . . completed! \n");


    ///////////////////////////////////////////////
    // DEVICE MEMORY                             //
    ///////////////////////////////////////////////

    /* 2D spatial variables. */
    SafeCudaCall(cudaMalloc((void**)&overland_dev->a2d        , 5*sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->hpoten       , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->ke           , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->kw           , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->kn           , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->ks           , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->mann         , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->ph           , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->qcapa        , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->rhs2d        , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->u            , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->v            , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->waterdepth   , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->waterelev    , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&overland_dev->ztopo        , sizexy*sizeof(double)));

    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcb        , sizexy*sizeof(int)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bct        , sizexy*sizeof(int)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bce        , sizeyz*sizeof(int)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcw        , sizeyz*sizeof(int)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcn        , sizexz*sizeof(int)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcs        , sizexz*sizeof(int)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->procmap    , sizexy*sizeof(int)));

    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->qss        , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcpsib     , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcpsit     , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcqb       , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcqt       , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcpsie     , sizeyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcpsiw     , sizeyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcqe       , sizeyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcqw       , sizeyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcpsin     , sizexz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcpsis     , sizexz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcqn       , sizexz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->bcqs       , sizexz*sizeof(double)));

    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->TR         , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->ppt_ground , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->E_soil     , sizexy*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->TR_root    , procsize*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->ppt_root   , procsize*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->E_soil_root, procsize*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->rda        , procsize*sizez*sizeof(double)));

    /* 3D spatial variables. */
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->a3d      , 7*sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->cnp1m      , sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->deltam     , sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->knp1m      , sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->ksat       , sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->psin       , sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->psinp1m    , sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->psinp1mp1  , sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->psiout     , sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->rhs3d      , sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->thetan     , sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->thetanp1m  , sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->thetanp1mp1, sizexyz*sizeof(double)));
    SafeCudaCall(cudaMalloc((void**)&subsurface_dev->thetaout   , sizexyz*sizeof(double)));

    /* If reach this point, print out info. */
    printf("\t Flow model - device memory . . . . . . . . . . . . . . completed! \n");
}


/**
 * @brief      Allocate memory for forcing variables.
 *
 * @param      timeforcings  Class including time forcings
 * @param[in]  rank          Global rank of the current MPI process
 * @param[in]  procsize      Total number of MPI processes available
 * @param[in]  num_steps     Number of time steps for simulation
 */
void AllocateMemoryForcing(TimeForcingClass *timeforcings, int rank, int procsize, int num_steps)
{
    // 1D temporal variables
    timeforcings->doy    = new int[num_steps];
    timeforcings->years  = new int[num_steps];
    timeforcings->rg     = new double[num_steps];
    timeforcings->pa     = new double[num_steps];
    timeforcings->lwdn   = new double[num_steps];
    timeforcings->zen    = new double[num_steps];
    timeforcings->u      = new double[num_steps];
    timeforcings->ppt    = new double[num_steps];
    timeforcings->ta     = new double[num_steps];
    timeforcings->ea     = new double[num_steps];
    timeforcings->lai    = new double[num_steps];
    timeforcings->vpd    = new double[num_steps];
    timeforcings->decdoy = new double[num_steps];
    timeforcings->hour   = new double[num_steps];
}



/**
 * @brief      Free all memory in host
 *
 * @param      timeforcings    Class including time forcings
 * @param      forcings        Class including forcings info
 * @param      overland        Overland flow class
 * @param      subsurface      Subsurface flow class
 * @param      switches        Class including switch info
 * @param      constants       Class including constant info
 * @param      canopies        Class including canopy variables
 * @param      soils           Class including soil variables
 * @param      radiation       Class including radiation variables
 * @param      photosynthesis  Class including photosynthesis variables
 * @param      respiration     Class including respiration variables
 * @param      stomaconduct    Class including stomaconductance variables
 * @param      microenviron    Class including microenrivonment variables
 * @param[in]  rank            Global rank of the current MPI process
 * @param[in]  procsize        Total number of MPI processes available
 */
void FreeHostMemory(TimeForcingClass *timeforcings, ForcingClass *forcings,
                    OverlandFlowClass *overland, SubsurfaceFlowClass *subsurface,
                    SwitchClass *switches, ConstantClass *constants,
                    CanopyClass *canopies, SoilClass *soils,
                    RadiationClass *radiation, PhotosynthesisClass *photosynthesis,
                    RespirationClass *respiration, StomaConductClass *stomaconduct,
                    MicroEnvironmentClass *microenviron, int rank, int procsize)
{
    if (rank == MPI_MASTER_RANK)
    {
        free(overland);
        free(subsurface);
    }
    free(timeforcings);
    free(forcings);
    free(switches);
    free(constants);
    free(canopies);
    free(soils);
    free(radiation);
    free(photosynthesis);
    free(respiration);
    free(stomaconduct);
    free(microenviron);
}



/**
 * @brief      To free all memory in device
 *
 * @param      overland    Overland flow class
 * @param      subsurface  Subsurface flow class
 * @param[in]  rank        Global rank of the current MPI process
 * @param[in]  procsize    Total number of MPI processes available
 */
void FreeDeviceMemory(OverlandFlowClass *overland, SubsurfaceFlowClass *subsurface,
                      int rank, int procsize)
{
    if (rank == MPI_MASTER_RANK)
    {
        ///////////////////////////
        // 2D overland flow      //
        ///////////////////////////

        // a2d, rhs2d, waterelev have been freed by CUSP
        SafeCudaCall(cudaFree(overland->hpoten));
        SafeCudaCall(cudaFree(overland->ke));
        SafeCudaCall(cudaFree(overland->kn));
        SafeCudaCall(cudaFree(overland->ks));
        SafeCudaCall(cudaFree(overland->kw));
        SafeCudaCall(cudaFree(overland->mann));
        SafeCudaCall(cudaFree(overland->ph));
        SafeCudaCall(cudaFree(overland->qcapa));
        SafeCudaCall(cudaFree(overland->u));
        SafeCudaCall(cudaFree(overland->v));
        SafeCudaCall(cudaFree(overland->waterdepth));
        SafeCudaCall(cudaFree(overland->ztopo));


        ///////////////////////////
        // 3D subsurface flow    //
        ///////////////////////////

        // a3d, rhs3d, psinp1mp1, deltam have been freed by CUSP
        SafeCudaCall(cudaFree(subsurface->bcb));
        SafeCudaCall(cudaFree(subsurface->bce));
        SafeCudaCall(cudaFree(subsurface->bcn));
        SafeCudaCall(cudaFree(subsurface->bcs));
        SafeCudaCall(cudaFree(subsurface->bct));
        SafeCudaCall(cudaFree(subsurface->bcw));
        SafeCudaCall(cudaFree(subsurface->bcpsib));
        SafeCudaCall(cudaFree(subsurface->bcpsie));
        SafeCudaCall(cudaFree(subsurface->bcpsin));
        SafeCudaCall(cudaFree(subsurface->bcpsis));
        SafeCudaCall(cudaFree(subsurface->bcpsit));
        SafeCudaCall(cudaFree(subsurface->bcpsiw));
        SafeCudaCall(cudaFree(subsurface->bcqb));
        SafeCudaCall(cudaFree(subsurface->bcqe));
        SafeCudaCall(cudaFree(subsurface->bcqn));
        SafeCudaCall(cudaFree(subsurface->bcqs));
        SafeCudaCall(cudaFree(subsurface->bcqt));
        SafeCudaCall(cudaFree(subsurface->bcqw));
        SafeCudaCall(cudaFree(subsurface->cnp1m));
        SafeCudaCall(cudaFree(subsurface->knp1m));
        SafeCudaCall(cudaFree(subsurface->ksat));
        SafeCudaCall(cudaFree(subsurface->psin));
        SafeCudaCall(cudaFree(subsurface->psinp1m));
        SafeCudaCall(cudaFree(subsurface->psiout));
        SafeCudaCall(cudaFree(subsurface->thetan));
        SafeCudaCall(cudaFree(subsurface->thetanp1m));
        SafeCudaCall(cudaFree(subsurface->thetaout));
    }
}


///////////////////////////////////////////////////////
// MLCan Models                                      //
///////////////////////////////////////////////////////


/**
 * @brief      Sets up ml can model.
 *
 * @param[in]  rank            Global rank of the current MPI process
 * @param[in]  procsize        Total number of MPI processes available
 * @param[in]  globsize        Size of global computational domain
 * @param[in]  domsize         Size of local domain in each MPI process
 * @param[in]  topolindex      The 2D index in the desired topology
 * @param      files           Class including files for loading data
 * @param      project         Class including project info
 * @param      switches        Class including switches
 * @param      constants       Class including constants
 * @param      canopies        Class including canopies
 * @param      soils           Class including soils
 * @param      radiation       Class including radiation
 * @param      photosynthesis  Class including photosynthesis
 * @param      respiration     Class including respiration
 * @param      stomaconduct    Class including stomaconduct
 * @param      microenviron    Class including microenviron
 * @param      cartComm        Carthesian MPI communicator
 */
void SetUpMLCanModel(int rank, int procsize, int3 globsize, int3 domsize, int2 topolindex,
                     FileNameClass *files, ProjectClass *project, SwitchClass *switches,
                     ConstantClass *constants, CanopyClass *canopies, SoilClass *soils,
                     RadiationClass *radiation, PhotosynthesisClass *photosynthesis,
                     RespirationClass *respiration, StomaConductClass *stomaconduct,
                     MicroEnvironmentClass *microenviron, MPI_Comm *cartComm)
{
    char plantname[64];

    // Assign plant to each process.
    if (topolindex.y < 2)
    {
        project->planttype = 0;
    }
    else
    {
        if (topolindex.x < 2)
        {
            project->planttype = 1;
        } else {
            project->planttype = 2;
        }
    }

    snprintf(plantname, sizeof(char) * 64, "PLANT_NAME_%d", project->planttype);
    files->plants = GetOptionToChar(plantname);

    if (project->verbose)
    {
        printf("\t rank: %3d \t Plant type: %3d \t Plant filename: %s \n", 
                rank, project->planttype, files->plants);
    }
    // Load MLCan configuration for each plant.
    LoadConfigMLCanModel(project, files, switches, constants, canopies, soils, radiation,
                         photosynthesis, respiration, stomaconduct, microenviron);

    MPI_Barrier(*cartComm); // sync all process.
}


/**
 * @brief      Two-way mapping eigen classes to regular classes in host memory.
 *             This allows changes in one are automatically linked to the other.
 *             See http://eigen.tuxfamily.org/
 *
 * @param      project        Class including project information
 * @param      timeforcings   Forcing variables over the entire simulation period.
 * @param      vertcanopies   Canopy classes including vertical variables
 * @param      eigencanopies  Canopy class defined in eigen structure
 * @param      vertsoils      Soil classes including vertical variables
 * @param      eigensoils     Soil class defined in eigen structure
 * @param[in]  nl_can         Number of canopy layers
 * @param[in]  nl_soil        Number of soil layers
 */
void MappingEigensToClasses(ProjectClass *project, TimeForcingClass *timeforcings,
                            VerticalCanopyClass *vertcanopies, EigenCanopyClass *eigencanopies,
                            VerticalSoilClass *vertsoils, EigenSoilClass *eigensoils, int nl_can,
                            int nl_soil)
{
    for (int i = 0; i < nl_can; i++)
    {
        vertcanopies->Tl_sun[i] = timeforcings->ta[0];
        vertcanopies->gsv_sun[i] = 0.01;
        vertcanopies->Ci_sun[i] = 0.7 * project->co2concentration;
        vertcanopies->Tl_shade[i] = timeforcings->ta[0];
        vertcanopies->gsv_shade[i] = 0.01;
        vertcanopies->Ci_shade[i] = 0.7 * project->co2concentration;
        vertcanopies->TR[i] = 0.0;
        vertcanopies->Sh2o_prof[i] = 0.0;
    }

    new (&eigencanopies->LAIsun) Map<VectorXd>(vertcanopies->LAIsun, nl_can);
    new (&eigencanopies->LAIshade) Map<VectorXd>(vertcanopies->LAIshade, nl_can);
    new (&eigensoils->Ts) Map<VectorXd>(vertsoils->Ts, nl_soil);
    eigensoils->Ts <<
        16.1925000000000, 18.4725000000000, 18.0112500000000, 17.2600000000000,
        16.3600000000000, 15.4600000000000, 14.6092187500000, 14.2014062500000,
        13.7935937500000, 13.3857812500000, 12.9779687500000, 12.5701562500000,
        12.1623437500000, 11.7545312500000, 11.3467187500000;
}


/**
 * @brief      Allocate memory to store outputs from MLCan model.
 *             Outputs may contain variables at different sizes.
 *
 * @param      canopies   Class including variables in canopy
 * @param      soils      Class including variables in soil
 * @param      outmlcan   Output class from MLCan model
 * @param[in]  num_steps  Number of simulation steps
 */
void AllocateMemoryOutput(CanopyClass *canopies, SoilClass *soils, OutputClass *outmlcan,
                          int num_steps)
{
    int nl_can = canopies->nl_can;
    int nl_soil = soils->nl_soil;

    outmlcan->An_can     = new double[num_steps];
    outmlcan->LE_can     = new double[num_steps];
    outmlcan->H_can      = new double[num_steps];
    outmlcan->TR_can     = new double[num_steps];
    outmlcan->Rnrad_can  = new double[num_steps];

    outmlcan->mbw_can    = new double[num_steps];

    outmlcan->E_soil     = new double[num_steps];
    outmlcan->qss        = new double[num_steps];
    outmlcan->htb        = new double[num_steps];
    outmlcan->diff       = new double[num_steps];
    outmlcan->PH         = new double[num_steps];
    outmlcan->G          = new double[num_steps];
    outmlcan->T_surf     = new double[num_steps];
    outmlcan->ppt_ground = new double[num_steps];

    outmlcan->An_sun     = new double[num_steps*nl_can];
    outmlcan->An_shade   = new double[num_steps*nl_can];
    outmlcan->LE_sun     = new double[num_steps*nl_can];
    outmlcan->LE_shade   = new double[num_steps*nl_can];
    outmlcan->H_sun      = new double[num_steps*nl_can];
    outmlcan->H_shade    = new double[num_steps*nl_can];
    outmlcan->Evap_prof  = new double[num_steps*nl_can];
    outmlcan->TR_sun     = new double[num_steps*nl_can];
    outmlcan->TR_shade   = new double[num_steps*nl_can];
    outmlcan->LAI_sun    = new double[num_steps*nl_can];
    outmlcan->LAI_shade  = new double[num_steps*nl_can];
    outmlcan->zhc        = new double[num_steps*nl_can];

    outmlcan->volliq     = new double[num_steps*nl_soil];
    outmlcan->smp        = new double[num_steps*nl_soil];
    outmlcan->krad       = new double[num_steps*nl_soil];
    outmlcan->rpp        = new double[num_steps*nl_soil];
    outmlcan->Ts         = new double[num_steps*nl_soil];
    outmlcan->zhs        = new double[num_steps*nl_soil];
}


///////////////////////////////////////////////////////
// Numerical Model gateway functions                 //
///////////////////////////////////////////////////////

/**
 * @brief      Main kernels calling numerical models including MLCan and Flow3D model. 
 *             While MLCan is implemented in all MPI processes in host, Flow3D model is run in GPU
 *             device controlled by the Master process. Communications between MPI and GPU for
 *             link MLCan and Flow3D is executed in Master process only. Simulations in other MPI
 *             processes are gatherred/scatterred to mMaster process every time step for 
 *             communicating with GPU device.
 *
 * @param[in]  rank      Global rank of the current MPI process
 * @param[in]  procsize  Total number of MPI processes available
 * @param      files     Number of time steps for simulation
 * @param      project   Filename for loading data
 * @param      mpiobj    The message passing interface objects
 */
void NumericalModelKernels(int rank, int procsize, FileNameClass *files, ProjectClass *project,
                           mpiClass *mpiobj)
{
    int num_steps     = project->num_steps;
    int3 globsize     = mpiobj->global_size;
    int3 domsize      = mpiobj->domain_size;
    int2 topolsize    = mpiobj->topology_size;
    int2 topolindex   = mpiobj->topology_index;
    MPI_Comm cartComm = mpiobj->cartComm;

    // Common time forcing variables.
    TimeForcingClass *timeforcings      = new TimeForcingClass;;

    // Classes in Flow model - host and device memory.
    OverlandFlowClass *olf_h            = new OverlandFlowClass;
    OverlandFlowClass *olf_d            = new OverlandFlowClass;
    SubsurfaceFlowClass *ssf_h          = new SubsurfaceFlowClass;
    SubsurfaceFlowClass *ssf_d          = new SubsurfaceFlowClass;

    // Classes in MLCan - host memory only.
    ForcingClass *mforcings             = new ForcingClass;;
    SwitchClass *switches               = new SwitchClass;
    ConstantClass *constants            = new ConstantClass;
    CanopyClass *canopies               = new CanopyClass;
    SoilClass *soils                    = new SoilClass;
    RadiationClass *radiation           = new RadiationClass;
    PhotosynthesisClass *photosynthesis = new PhotosynthesisClass;
    RespirationClass *respiration       = new RespirationClass;
    StomaConductClass *stomaconduct     = new StomaConductClass;
    MicroEnvironmentClass *microenviron = new MicroEnvironmentClass;
    OutputClass *outmlcan               = new OutputClass;


    // MLCan model set up.
    SetUpMLCanModel(rank, procsize, globsize, domsize, topolindex, files, project, switches,
                    constants, canopies, soils, radiation, photosynthesis, respiration,
                    stomaconduct, microenviron, &cartComm);

    // Forcing file for both models.
    AllocateMemoryForcing(timeforcings, rank, procsize, num_steps);
    LoadForcingData(files, canopies, timeforcings, rank, procsize, num_steps);

    // Vertical classes in canopy and soil.
    VerticalCanopyClass *vertcanopies = new VerticalCanopyClass(canopies->nl_can);
    VerticalSoilClass *vertsoils      = new VerticalSoilClass(canopies->nl_can);
    EigenCanopyClass *eigencanopies   = new EigenCanopyClass(canopies->nl_can);
    EigenSoilClass *eigensoils        = new EigenSoilClass(soils->nl_soil);

    // Load distributions and set up grids.
    LoadCanopyRootDistributions(project, canopies, vertcanopies, soils, vertsoils);

    MappingEigensToClasses(project, timeforcings, vertcanopies, eigencanopies, vertsoils,
                           eigensoils, canopies->nl_can, soils->nl_soil);

    AllocateMemoryOutput(canopies, soils, outmlcan, num_steps);     // for storing results
                                                                    // null vars for now

    // Set up flow model on device by master.
    if (rank == MPI_MASTER_RANK)
    {
        AllocateMemoryFlowModel(timeforcings, olf_h, olf_d, ssf_h, ssf_d, procsize, globsize,
                                num_steps);
        LoadTopography(files, olf_h, rank, procsize);
        LoadFlowModelConfig(project, files, olf_h, ssf_h);
        SetFlowModelConditions(timeforcings, files, olf_h, ssf_h, num_steps, rank, procsize,
                               globsize);
        SetMPIGPUMapping(mpiobj, ssf_h, procsize, globsize, topolsize);
        CopyDataToDevice(project, timeforcings, olf_h, olf_d, ssf_h, ssf_d, globsize);
    }
    MPI_Barrier(cartComm);     // sync all process

    GatherRootDensityToMasterAndFilter(project, vertsoils, soils, ssf_h, ssf_d, rank, procsize,
                                       &cartComm);

    /*
     * Run main numerical models.
     * All models are run within this function.
     * Time loop and export funcs are inside as well.
     */
    RunCoupledFlowModel(timeforcings, olf_h, olf_d, ssf_h, ssf_d, files, project, mforcings,
                        switches, constants, canopies, soils, radiation, photosynthesis,
                        respiration, stomaconduct,microenviron, vertcanopies, vertsoils,
                        eigencanopies, eigensoils, outmlcan, rank, procsize, globsize, domsize,
                        topolsize, topolindex, &cartComm);

    // Free all memory after completion.
    FreeHostMemory(timeforcings, mforcings, olf_h, ssf_h, switches, constants, canopies, soils,
                   radiation, photosynthesis, respiration, stomaconduct, microenviron, rank,
                   procsize);

    FreeDeviceMemory(olf_d, ssf_d, rank, procsize);
    MPI_Barrier(cartComm);     // sync all process
}
