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

template<typename type>
void SaveOutput1D(const char *file, const char *dataname, type *data, nc_type nctype, 
                  int Nx, int write);


// Print the usage information for this application
void PrintUsage(const char * appName)
{
    printf("\n");
    printf("Dhara Modeling System \n\n");

    printf("Usage: mpirun -np [numprocs] %s -[options] \n\n", appName);

    printf("Options:\n");
    printf("    Startup:\n");
    printf("    -V,  --version          Print version information of MLCan-Flow3D \n");
    printf("    -h,  --help             Print help information. \n");

    printf("\n");
    printf("    Mandatory:\n");
    printf("    -c [config file]        Choose configuration file \n");
    printf("    -t [topo.x] [topo.y]    Set MPI topology size \n");
    printf("                            if \"topo.y\" is missing, it will be set to 1 \n");
    printf("                            topo.x and topo.y must be positive integers\n\n");

    printf("    Optional:\n");
    printf("    -q,  --quite            quite mode (no display) \n");
    printf("    -v,  --verbose          Enable and use verbose mode \n\n");

    printf("Documentation:\n");
    printf("    http://hydrocomplexity.github.io/Dhara \n");
    printf("    https://github.com/HydroComplexity/Dhara \n\n\n");
}


// Print the usage information for this application
void PrintVersion()
{
    printf("Dhara Modeling System \n");
    printf("version 1.0 (Build 0316) \n");
    printf("Copyright (c) 2016 HydroComplexity Group \n");
    printf("\n");
}


// Check argument in the command line
int CheckArgument(const char * argName, int argc, char ** argv)
{
    for(int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], argName) == 0)
        {
            strcpy(argv[i], "");
            return i;
        }
    }
    return STATUS_ERR;
}


/**
 * Extract a number given as a command-line argument
 *
 * @param[in] argIdx        Argument index in command line
 * @param[in] argc          The number of input arguments
 * @param[in] argv          The input arguments
 * @return                  Value of the the argument
 */
int ExtractNumber(int argIdx, int argc, char ** argv)
{
    int result = 0;

    if (argIdx < argc){
        result = atoi(argv[argIdx]);
        if (result > 0) {
            strcpy(argv[argIdx], "");
        }
    }
    return result;
}


/**
 * @brief      Loads user configurations to MLCanmodel.
 *
 * @param      project         Class including project info
 * @param      files           Class including filename info
 * @param      switches        Class including switch info
 * @param      constants       Class including constant info
 * @param      canopies        Class including canopy variables
 * @param      soils           Class including soil variables
 * @param      radiation       Class including radiation variables
 * @param      photosynthesis  Class including photosynthesis variables
 * @param      respiration     Class including respiration variables
 * @param      stomaconduct    Class including stomaconductance variables
 * @param      microenviron    Class including microenrivonment variables
 */
void LoadConfigMLCanModel(ProjectClass *project, FileNameClass *files, SwitchClass *switches,
                          ConstantClass *constants, CanopyClass *canopies, SoilClass *soils,
                          RadiationClass *radiation, PhotosynthesisClass *photosynthesis,
                          RespirationClass *respiration, StomaConductClass *stomaconduct,
                          MicroEnvironmentClass *microenviron)
{
    // Parsing configuration file
    ParserConfigFile(files->plants);

  // SWITCHES . . . . . . . . .
    switches->PH_type                 = GetOptionToInt("PHOTOSYNTHESIS_TYPE");
    switches->Turbulence              = GetOptionToInt("TURBULENCE");
    switches->HydraulicRedistribution = GetOptionToInt("HYDRAULIC_REDISTRIBUTION");
    switches->RootConductivity        = GetOptionToInt("ROOT_CONDUCTIVITY");
    switches->RootHydrauConduct       = GetOptionToInt("ROOT_HYDRAUCONDUCT");
    switches->SoilHeat                = GetOptionToInt("SOIL_HEAT_MODEL");
    switches->LWequation              = GetOptionToInt("LONGWAVE_EQUATION");

    // CONSTANTS . . . . . . . . .
    constants->umoltoWm2   = GetOptionToDouble("UMOL_TO_WM2");
    constants->Wm2toumol   = 1.0 / constants->umoltoWm2;
    constants->mmH2OtoMPa  = GetOptionToDouble("MMH2O_TO_MPA");
    constants->R           = GetOptionToDouble("GAS_CONSTANT");
    constants->R_kJ        = constants->R / 1000.0;
    constants->Lv          = GetOptionToDouble("LATENT_HEAT_VAPORIZATION");
    constants->Lv_g        = GetOptionToDouble("LATENT_HEAT_VAPORIZATION_GRAM");
    constants->cp_mol      = GetOptionToDouble("SPECIFIC_HEAT_MOL");
    constants->cp_JkgK     = GetOptionToDouble("SPECIFIC_HEAT_JKGK");
    constants->boltz       = GetOptionToDouble("STEFAN_BOLTZMANN_CONSTANT");
    constants->vonk        = GetOptionToDouble("VONKARMAN_CONSTANT");
    constants->rho_dry_air = GetOptionToDouble("DENSITY_DRY_AIR");
    constants->grav        = GetOptionToDouble("GRAVITY");
    constants->timestep    = project->dtimehr * 60;         // minute
    constants->dtime       = constants->timestep * 60;      // second

    // CANOPY STRUCTURE . . . . . .
    canopies->name         = GetOptionToChar("PLANT_NAME");
    canopies->dists        = GetOptionToChar("DISTRIBUTIONS_FILE");
    canopies->nl_can       = GetOptionToDouble("CANOPY_LAYERS");
    canopies->LEfact       = GetOptionToInt("LATENT_FACT");
    canopies->Hfact        = GetOptionToInt("SENSIBLE_FACT");
    canopies->LWfact       = GetOptionToInt("LONGWAVE_FACT");
    canopies->hcan         = GetOptionToDouble("CANOPY_HEIGHT");
    canopies->z0           = 0.13 * canopies->hcan;
    canopies->d0           = 2.0 / 3.0 * canopies->hcan;
    canopies->leaftype     = GetOptionToInt("LEAF_TYPE");
    canopies->ld           = GetOptionToDouble("LEAF_DIMENSION");
    canopies->lw           = GetOptionToDouble("LEAF_WIDTH");
    canopies->Smax         = GetOptionToDouble("STORAGE_MAX");
    canopies->Ffact        = GetOptionToDouble("WET_FRACTION");
    canopies->pptintfact   = GetOptionToDouble("PPT_EXTINCTION_COEFF");

    // SOIL . . . . . . . . . .
    soils->nl_soil         = NUM_SOIL_LAYERS;
    soils->z0              = GetOptionToDouble("SURFACE_ROUGHNESS_LENGTH");
    soils->sand            = GetOptionToDouble("SAND_PERCENTAGE");
    soils->clay            = GetOptionToDouble("CLAY_PERCENTAGE");

    // RADIATION . . . . . . . . . .
    radiation->transmiss   = GetOptionToDouble("ATMOSPHERIC_TRANSMISSIVITY");
    radiation->epsv        = GetOptionToDouble("VEGETATION_EMISSIVITY");
    radiation->epss        = GetOptionToDouble("SOIL_EMISSIVITY");
    radiation->epsa        = GetOptionToDouble("ATMOSPHERIC_EMISSIVITY");
    radiation->xx          = GetOptionToDouble("LEAF_ANGLE_DISTRIBUTION");
    radiation->clump       = GetOptionToDouble("LEAF_CLUMP");
    radiation->Kdf         = GetOptionToDouble("EXTINCTION_COEFF_DIFFUSION");
    radiation->absorp_PAR  = GetOptionToDouble("PAR_ABSORPTIVITY");
    radiation->absorp_NIR  = GetOptionToDouble("NIR_ABSORPTIVITY");
    radiation->refl_PAR    = (1.0-sqrt(radiation->absorp_PAR)) / (1.0+sqrt(radiation->absorp_PAR));
    radiation->refl_NIR    = (1.0-sqrt(radiation->absorp_NIR)) / (1.0+sqrt(radiation->absorp_NIR));
    radiation->refl_soil   = GetOptionToDouble("SOIL_REFLECTION");
    radiation->trans_PAR   = 1.0 - radiation->absorp_PAR - radiation->refl_PAR;
    radiation->trans_NIR   = 1.0 - radiation->absorp_NIR - radiation->refl_NIR;;

    // PHOTOSYNTHESIS . . . . . . . .
    photosynthesis->ph_type      = switches->PH_type;
    photosynthesis->Vmax_C4      = GetOptionToDouble("SATURATED_RUBISCO_CAPACITY");
    photosynthesis->Rd_C4        = GetOptionToDouble("LEAF_RESPIRATION_C4");
    photosynthesis->Q10_C4       = GetOptionToDouble("Q10_C4");
    photosynthesis->kk_C4        = GetOptionToDouble("KK_C4");
    photosynthesis->theta_C4     = GetOptionToDouble("THETA_C4");
    photosynthesis->beta_C4      = GetOptionToDouble("BETA_C4");
    photosynthesis->al_C4        = GetOptionToDouble("AL_C4");

    photosynthesis->beta_ph_C3   = GetOptionToDouble("FRACTION_ABSORP_Q");
    photosynthesis->Vcmax25_C3   = GetOptionToDouble("MAXIMUM_RUBISCO_AT_25");
    photosynthesis->Vcmax25_fact = GetOptionToDouble("MAXIMUM_RUBISCO25_FACTOR");
    photosynthesis->Jmax25_C3    = 2.0 * photosynthesis->Vcmax25_C3;
    photosynthesis->Rd25         = 0.015 * photosynthesis->Vcmax25_C3;
    
    photosynthesis->kn_canopy    = GetOptionToDouble("KN_CANOPY");
    photosynthesis->ap           = GetOptionToDouble("AP");
    photosynthesis->bp           = GetOptionToDouble("BP");
    photosynthesis->Oi           = GetOptionToDouble("INTERNAL_OXYGEN_CONCENTRATION");

    // RESPIRATION . . . . . . . . . . .
    respiration->Ro           = GetOptionToDouble("RESPIRATION_RATE_10C");
    respiration->Q10          = GetOptionToDouble("RESPIRATION_Q10");

    // STOMATA CONDUCTANCE . . . . . . .
    stomaconduct->mslope      = GetOptionToDouble("MSLOPE");
    stomaconduct->bint        = GetOptionToDouble("BINT");
    stomaconduct->sf          = GetOptionToDouble("SF");
    stomaconduct->psif        = GetOptionToDouble("LEAF_POTENTIAL_HALF");
    stomaconduct->Rp          = GetOptionToDouble("RESISTANCE_FLOW");

    // MICRO-ENVIRONMENT . . . . . . . .
    microenviron->Cd          = GetOptionToDouble("DRAG_COEFFICIENT");
    microenviron->alph        = constants->vonk/3.0;
}


/**
 * @brief      Loads user configurations to flow model.
 *
 * @param      project     Class including project info
 * @param      files       Class including filename info
 * @param      overland    Overland flow class
 * @param      subsurface  Subsurface flow class
 */
void LoadFlowModelConfig(ProjectClass *project, FileNameClass *files, OverlandFlowClass *overland,
                         SubsurfaceFlowClass *subsurface)
{
    // Parsing configuration file
    ParserConfigFile(files->config);

    // SUBSURFACE FLOW . . . . . . . . . .
    subsurface->maxiter         = GetOptionToInt("MAXIMUM_ITERATIONS");
    subsurface->picardmethod    = GetOptionToInt("PICARD_METHOD");
    subsurface->airdry          = GetOptionToDouble("DRY_AIR_PRESSURE");
    subsurface->tolerance_psi   = GetOptionToDouble("TOLERANCE_PRESSURE");
    subsurface->tolerance_theta = GetOptionToDouble("TOLERANCE_THETA");
    subsurface->psimin          = GetOptionToDouble("MINIMUM_PRESSURE_HEAD");
    subsurface->am              = GetOptionToDouble("WEIGHTED_SCHEME");
    subsurface->alpha           = GetOptionToDouble("ALPHA");
    subsurface->poresize        = GetOptionToDouble("PORE_SIZE_DISTRIBUTION");
    subsurface->porosity        = GetOptionToDouble("POROSITY");
    subsurface->Ss              = GetOptionToDouble("SPECIFIC_STORAGE");
    subsurface->theta_r         = GetOptionToDouble("RESIDUAL_WATER_CONTENT");
    subsurface->theta_s         = GetOptionToDouble("SATURATED_WATER_CONTENT");

    subsurface->dt              = GetOptionToDouble("TIME_STEP_DATA");
    subsurface->dx              = GetOptionToDouble("GRID_SIZE.DX");
    subsurface->dy              = GetOptionToDouble("GRID_SIZE.DY");
    subsurface->dz              = GetOptionToDouble("GRID_SIZE.DZ");

    // OVERLAND FLOW . . . . . . . . . .
    overland->delta = GetOptionToDouble("DELTA");
    overland->hmin  = GetOptionToDouble("WATERDEPTH_MINIMUM");
    overland->hcri  = GetOptionToDouble("WATERDEPTH_CRITICAL");
    overland->K0    = GetOptionToDouble("KZERO");
    overland->hn    = GetOptionToDouble("WATERDEPTH_NORTH");
    overland->hs    = GetOptionToDouble("WATERDEPTH_SOUTH");
    overland->hw    = GetOptionToDouble("WATERDEPTH_WEST");
    overland->he    = GetOptionToDouble("WATERDEPTH_EAST");
    overland->dt    = GetOptionToDouble("TIME_STEP_DATA");
    overland->dx    = GetOptionToDouble("GRID_SIZE.DX");
    overland->dy    = GetOptionToDouble("GRID_SIZE.DY");
}


/**
 * @brief      Load NetCDF file format
 *
 * @param[in]  file_name  NetCDF file name
 * @param[in]  var_name   Variable name
 * @param      data       Data assigned
 */
template<typename Type>
void LoadFileNetCDF(const char *file_name, const char *var_name, Type *data)
{
    int ncid, varid, status;

    // Open Netcdf file with NC_NOWRITE options
    status = nc_open(file_name, NC_NOWRITE, &ncid);
    if (status != NC_NOERR)
    {
        printf("Load NetCDF File Error: No file \"%s\" found. Exit program now! \n ", file_name);
        exit(1);
    }

    // Get variable ID
    status = nc_inq_varid(ncid, var_name, &varid);
    if (status != NC_NOERR)
    {
        printf("Error: No variable \"%s\" found in file \"%s\" found. Exit program now! \n ",
                var_name, file_name);
        exit(1);
    }
    nc_get_var(ncid, varid, &data[0]);          // Get data based on data_type

    nc_close(ncid);                             // Close the NetCDF file
}


/**
 * @brief      To extract information and dimension of NetCDF files.
 *
 * @param[in]  file_name  NetCDF file
 * @param[in]  var_name   Variable inspecting
 * @param[in]  ndims      Number of dimension of the vars
 * @param      dim        Dimension of the vars
 */
void GetFileInfo(const char *file_name, const char *var_name, int ndims, int *dim)
{
    int ncid, status, var_id, dimids[4];
    size_t length;

    // Open Netcdf file with NC_NOWRITE options
    status = nc_open(file_name, NC_NOWRITE, &ncid);
    if (status != NC_NOERR) {
        printf("GetFileInfo Error: No file \"%s\" found. Quitting program! \n ", file_name);
        exit(1);
    }

    // Get variable id, dimension, size
    nc_inq_varid(ncid, var_name, &var_id);
    nc_inq_vardimid(ncid, var_id, dimids);

    for (int i = 0; i < ndims; i++) {
        nc_inq_dimlen(ncid, dimids[i], &length);
        dim[i] = length;
    }

    // Close Netcdf file
    nc_close(ncid);
}


/**
 * @brief      Parse command line arguments and read configuration file
 *
 * @param[in]  rank       Global rank of the current MPI process
 * @param[in]  size       The total number of MPI processes available
 * @param      globsize   Size of the entire domain
 * @param      domsize    Size of local domain in each MPI
 * @param      topolsize  The desired topology size of MPI
 */
void ParsingCommandsAndConfigurations(int argc, char **argv, const char * &file_config, int rank,
                                      int procsize, ProjectClass * &project, 
                                      FileNameClass * &files, mpiClass * &mpiobj)
{
    int isroot = (rank == MPI_MASTER_RANK);
    int dim_topography[2];
    int argIdx;

    if (argc < 2) {
        if (isroot)
            PrintUsage(argv[0]);

        exit(STATUS_ERR);
    }

    // If help is requested, all other arguments will be ignored
    if ((CheckArgument("-h", argc, argv) != -1) || (CheckArgument("--help", argc, argv) != -1))
    {
        if (isroot)
            PrintUsage(argv[0]);

        // This simply prevents the application from continuing
        exit(STATUS_ERR);
    }


    // If version is requested, all other arguments will be ignored
    if ((CheckArgument("-V", argc, argv) != -1) || (CheckArgument("--version", argc, argv) != -1))
    {
        if (isroot)
            PrintVersion();

        // Prevents the application from continuing
        exit(STATUS_ERR);
    }


    // Verbose
    if ((CheckArgument("-v", argc, argv) != -1) || (CheckArgument("--verbose", argc, argv) != -1))
    {
        project->verbose = 1;
        if (isroot)
            printf("\t Verbose mode is on . . . . . .\n");        
    } else {
        project->verbose = 0;
    }

    // Configuration file must always be present
    argIdx = CheckArgument("-c", argc, argv);
    if (argIdx == STATUS_ERR)
    {
        OneErrPrintf(isroot, "Error: Could not find the configuration file information.\n");
        exit(STATUS_ERR);
    }
    else
    {
        file_config = argv[argIdx+1];
        files->config = argv[argIdx+1];
    }


    // Topology information must always be present
    argIdx = CheckArgument("-t", argc, argv);
    if (argIdx == STATUS_ERR)
    {
        OneErrPrintf(isroot, "Error: Could not find the topology information.\n");
        exit(STATUS_ERR);
    }
    else
    {
        mpiobj->topology_size.x = ExtractNumber(argIdx + 1, argc, argv);
        mpiobj->topology_size.y = ExtractNumber(argIdx + 2, argc, argv);

        // At least the first topology dimension must be specified
        if (mpiobj->topology_size.x <= 0)
        {
            OneErrPrintf(isroot, "Error: The topology size is invalid (first value: %d)\n", mpiobj->topology_size.x);
            exit(STATUS_ERR);
        }

        // If the second topology dimension is missing, set default to 1
        if (mpiobj->topology_size.y <= 0)
        {
            mpiobj->topology_size.y = 1;
        }
    }

    // The number of MPI processes must fill the topology
    int topologySize = mpiobj->topology_size.x * mpiobj->topology_size.y;
    if (procsize != topologySize)
    {
        OneErrPrintf(rank == MPI_MASTER_RANK,
            "Error: The number of MPI processes (%d) does not match topology size (%d x %d).\n",
            procsize, mpiobj->topology_size.x, mpiobj->topology_size.y);
        exit(STATUS_ERR);
    }

    // Parse config file to process
    ParserConfigFile(file_config);

    // General project
    project->dtimehr  = GetOptionToDouble("TIME_STEP_DATA");
    project->dx_meter = GetOptionToDouble("GRID_SIZE.DX");
    project->dy_meter = GetOptionToDouble("GRID_SIZE.DY");
    project->dz_meter = GetOptionToDouble("GRID_SIZE.DZ");

    // Number of plants modeled
    project->numberofplants   = GetOptionToInt("PLANT_NUMBERS");
    project->co2elevation     = GetOptionToInt("CO2_ELEVATION");
    project->co2concentration = GetOptionToDouble("CO2_CONCENTRATION");
    project->mlcanoutput      = GetOptionToChar("OUTPUT_NAME_MLCAN");
    project->savemlcan        = GetOptionToInt("SAVE_MLCAN");
    project->saveinterval     = GetOptionToInt("SAVE_INTERVAL");
    project->printiterval     = GetOptionToInt("PRINT_INTERVAL");

    // Load forcing files
    files->forcings = GetOptionToChar("FILE_FORCING");
    LoadFileNetCDF(files->forcings, "NumTimeSteps", &project->num_steps);

    // Specify folder to save results for all processes
    project->folderoutput = GetOptionToChar("OUTPUT_FOLDER_NAME");

    // Get topography from configuration file and specify outputs
    if (rank == MPI_MASTER_RANK)
    {
        files->topography = GetOptionToChar("FILE_TOPOGRAHY");
        GetFileInfo(files->topography, "Topography", 2, dim_topography);

        // Get domainSize then estimate size of global domain.
        mpiobj->global_size.x = dim_topography[1];
        mpiobj->global_size.y = dim_topography[0];
        mpiobj->global_size.z = NUM_SOIL_LAYERS;

        project->saveolf     = GetOptionToInt("SAVE_OVERLAND");
        project->savessf     = GetOptionToInt("SAVE_SUBSURFACE");
        project->savestat    = GetOptionToInt("SAVE_STATISTICS");
        project->olfoutput   = GetOptionToChar("OUTPUT_NAME_OVERLAND");
        project->ssfoutput   = GetOptionToChar("OUTPUT_NAME_SUBSURFACE");
        project->ssf1doutput = GetOptionToChar("OUTPUT_NAME_SUBSURFACE_1D");
        project->statoutput  = GetOptionToChar("OUTPUT_NAME_STATISTICS");
    }

    if (isroot)
    {
        printf("\n");
        printf("INITIALIZE MPI PARALLEL ENVIRONMENT . . . . . . . . . . . . . . completed! \n");
        printf("----------------------------------- \n");

        printf("\nMODEL INFORMATION AND PARAMETERS \n");
        printf("-------------------------------- \n");
        printf("\t Topology size: %d x %d\n", mpiobj->topology_size.x, 
                                              mpiobj->topology_size.y);
        printf("\t Global domain size (all nodes): %d x %d\n", mpiobj->global_size.x, 
                                                               mpiobj->global_size.y);
        printf("\t Number of soil layers: %d              \n", mpiobj->global_size.z);
        printf("\t Number of time steps simulation: %d \n \n", project->num_steps);
    }
}


/**
 * @brief      load forcing data from NetCDF file.
 *
 * @param      files         NetCDF file
 * @param      canopies      The canopies
 * @param      timeforcings  Forcings time class
 * @param[in]  rank          Global rank of the current MPI process
 * @param[in]  procsize      Total number of MPI processes available
 * @param[in]  num_steps     The number of timesteps
 */
void LoadForcingData(FileNameClass *files, CanopyClass *canopies, TimeForcingClass *timeforcings,
                     int rank, int procsize, int num_steps)
{
    int isroot = (rank == MPI_MASTER_RANK);
    if (isroot)
    {
        printf("\nLOADING FORCINGS \n");
        printf("------------------ \n");
    }

    LoadFileNetCDF(files->forcings , "AirTemperature"      , timeforcings->ta);
    LoadFileNetCDF(files->forcings , "AtmosphericPressure" , timeforcings->pa);
    LoadFileNetCDF(files->forcings , "DayOfYear"           , timeforcings->doy);
    LoadFileNetCDF(files->forcings , "DecimalDayOfYear"    , timeforcings->decdoy);
    LoadFileNetCDF(files->forcings , "GlobalRadiation"     , timeforcings->rg);
    LoadFileNetCDF(files->forcings , "Hour"                , timeforcings->hour);
    LoadFileNetCDF(files->forcings , "LongwaveDownward"    , timeforcings->lwdn);
    LoadFileNetCDF(files->forcings , "Precipitation"       , timeforcings->ppt);
    LoadFileNetCDF(files->forcings , "VaporPressureAir"    , timeforcings->ea);
    LoadFileNetCDF(files->forcings , "WindVelocity"        , timeforcings->u);
    LoadFileNetCDF(files->forcings , "Year"                , timeforcings->years);
    LoadFileNetCDF(files->forcings , "ZenithAngle"         , timeforcings->zen);
    LoadFileNetCDF(canopies->dists , "LeafAreaIndex"       , timeforcings->lai);

    double pptsum  = 0.0;
    double rgsum   = 0.0;
    double pasum   = 0.0;
    double lwdnsum = 0.0;
    double zensum  = 0.0;
    double usum    = 0.0;
    double tasum   = 0.0;
    double easum   = 0.0;
    double laisum  = 0.0;

    for (int i = 0; i < num_steps; ++i)
    {
        timeforcings->ppt[i]; // ppt in [mm] /= 1000.;
        pptsum += timeforcings->ppt[i];
        rgsum += timeforcings->rg[i];
        pasum += timeforcings->pa[i];
        lwdnsum += timeforcings->lwdn[i];
        zensum += timeforcings->zen[i];
        usum += timeforcings->u[i];
        tasum += timeforcings->ta[i];
        easum += timeforcings->ea[i];
        laisum += timeforcings->lai[i];
    }

    if (isroot)
    {
        //printf("\t Total precipitation [mm]:           % 5.2f \n", pptsum * 1000.);
        printf("\t Total precipitation [mm]:           % 5.2f \n", pptsum);
        printf("\t Mean global radiation [W/m^2]:      % 5.2f \n", rgsum/num_steps);
        printf("\t Mean air pressure [MPa]:            % 5.2f \n", pasum/num_steps);
        printf("\t Mean longwave downward [W/m^2]:     % 5.2f \n", lwdnsum/num_steps);
        printf("\t Mean zenith angle [degree]:         % 5.2f \n", zensum/num_steps);
        printf("\t Mean wind velocity [m/s]:           % 5.2f \n", usum/num_steps);
        printf("\t Mean air temperature [^oC]:         % 5.2f \n", tasum/num_steps);
        printf("\t Mean pressure of water vapor [kPa]: % 5.2f \n", easum/num_steps);
    }

}


/**
 * @brief      To load topography info from NetCDF file.
 *             Topography is only used for flow modling in GPU device through master process.
 *
 * @param      files     NetCDF file
 * @param      overland  Overland flow class
 * @param[in]  rank      Global rank of the current MPI process
 * @param[in]  procsize  Total number of MPI processes available
 */
void LoadTopography(FileNameClass *files, OverlandFlowClass *overland, int rank, int procsize)
{
    printf("\nLOADING TOPOGRAPHY \n");
    printf("-------------------- \n");
    LoadFileNetCDF(files->topography, "Topography", overland->ztopo);
    printf("\t Flow model . . . . . completed! \n");
}



/**
 * @brief      Set up canopy grid for simulation. Grid is constructed based
 *             on canopy height and the number of canopy layers specified by users.
 *             Uniform canopy grid is used.
 *
 * @param      canopies      Class including variables in canopy
 * @param      vertcanopies  Canopy classes including vertical variables
 */
void SetCanopyGrid(CanopyClass *canopies, VerticalCanopyClass *vertcanopies)
{
    int nl_can = canopies->nl_can;
    double hcan = canopies->hcan;

    for (int i=0; i< nl_can; i++)
        vertcanopies->zhc[i] = (i+1) * (hcan / nl_can);

    vertcanopies->dzc = 0.5 * vertcanopies->zhc[0];

    for (int i=0; i< nl_can; i++)
        vertcanopies->znc[i] = vertcanopies->zhc[i] - vertcanopies->dzc;
}


void SetSoilGrid( ProjectClass *project, SoilClass *soils, VerticalSoilClass *vertsoils)
{
    int nl_soil = soils->nl_soil;
    double rhod;
    double sand = 5.0;
    double clay = 25.0;

    vertsoils->zhs[0] = project->dz_meter;
    for (int i=1; i<nl_soil; i++)
    {
        vertsoils->zhs[i] = vertsoils->zhs[i] + project->dz_meter;
    }

    for (int i=0; i<nl_soil; i++)
    {
        vertsoils->zns[i] = vertsoils->zhs[i] - 0.5 * project->dz_meter;
        vertsoils->dzs[i] = project->dz_meter;
        vertsoils->sand[i] = sand;
        vertsoils->clay[i] = clay;
        vertsoils->HC_sol[i] = (2.128*sand + 2.385*clay) / (sand+clay) * 1e6;
        vertsoils->porosity[i] = 0.45;
        vertsoils->psi0[i] = -10 * pow( 10, 1.88-0.0131*sand );
        vertsoils->TK_sol[i] = (8.80*sand + 2.92*clay) / (sand+clay);
        rhod = 2700*(1 - vertsoils->porosity[i]);
        vertsoils->TK_dry[i]  = (0.135*rhod + 64.7) / (2700 - 0.947*rhod);
    }
}

/**
 * @brief      Load canopy and root distribution from text files specified in
 *             configuration file. If more than one vegetation, distributions
 *             for each species must be provided.
 *
 * @param      canopies      Class including variables in canopy
 * @param      vertcanopies  Canopy classes including vertical variables
 * @param      soils         Class including variables in soil
 * @param      vertsoils     Soil classes including vertical variables
 */
void LoadCanopyRootDistributions(ProjectClass *project, CanopyClass *canopies, 
                                 VerticalCanopyClass *vertcanopies, SoilClass *soils,
                                 VerticalSoilClass *vertsoils)
{
    LoadFileNetCDF(canopies->dists, "LeafAreaDensity", vertcanopies->LADnorm);
    LoadFileNetCDF(canopies->dists, "RootFraction", vertsoils->rootfr);

    SetCanopyGrid(canopies, vertcanopies);
    SetSoilGrid(project, soils, vertsoils);
    // TODO: check if nl_can and nl_soil are matched with LAD and rootfr
}



///////////////////////////////////
// SAVE DATA                     //
///////////////////////////////////

/**
 * @brief      Save 3D outputs to NetCDF files in double format.
 *
 * @param[in]  file      File will be saved in NetCDF format
 * @param[in]  dataname  Variable name saved
 * @param      data      The dataset
 * @param[in]  globsize  Size of global computational domain
 * @param[in]  write     Option to decide write new or append
 */
template<typename type>
void SaveOutput3D(const char *file, const char *dataname, type *data, nc_type nctype, 
                  int3 globsize, int write)
{
    int ncid, x_dimid, y_dimid, z_dimid, varid, dimids[3], status;
    int sizex  = globsize.x;
    int sizey  = globsize.y;
    int sizez  = globsize.z;

    // Set up NetCDF file for writing
    if (write == 0)
    {
        // Create a new NetCDF file
        status = nc_create(file, NC_CLOBBER, &ncid);
        if (status != NC_NOERR)
        {
            printf("Error in %s: No file \"%s\" found. No data is save. Quitting program! \n ",
                    __func__, file);
        }

    } else {
        // Open and re-define an existing NetCDF file
        status = nc_open(file, NC_WRITE, &ncid);
        if (status != NC_NOERR)
        {
            printf("Error in %s: No file \"%s\" found. No data is save. Quitting program! \n ",
                    __func__, file);
        }
        nc_redef(ncid);
    }

    nc_def_dim(ncid, "x", sizex, &x_dimid);
    nc_def_dim(ncid, "y", sizey, &y_dimid);
    nc_def_dim(ncid, "z", sizez, &z_dimid);

    dimids[0] = z_dimid;
    dimids[1] = y_dimid;
    dimids[2] = x_dimid;

    nc_def_var(ncid, dataname, nctype, 3, dimids, &varid);
    nc_enddef(ncid);
    nc_put_var(ncid, varid, &data[0]);
    nc_close(ncid);
}



/**
 * @brief      Save 2D outputs to NetCDF files in double format.
 *
 * @param[in]  file      File will be saved in NetCDF format
 * @param[in]  dataname  Variable name saved
 * @param      data      The dataset
 * @param[in]  nctype    NetCDF data type used
 * @param[in]  My        Data dimension in y
 * @param[in]  Nx        Data dimension in x
 * @param[in]  write     Option to decide write new or append
 */
template<typename type>
void SaveOutput2D(const char *file, const char *dataname, type *data, nc_type nctype, 
                  int My, int Nx, int write)
{
    int ncid, x_dimid, y_dimid, varid, dimids[2], status;

    // Set up NetCDF file for writing
    if (write == 0)
    {
        // Create a new NetCDF file
        status = nc_create(file, NC_CLOBBER, &ncid);
        if (status != NC_NOERR)
        {
            printf("Error in %s: No file \"%s\" found. No data is save. Quitting program! \n ",
                    __func__, file);
        }

    } else {
        // Open and re-define an existing NetCDF file
        status = nc_open(file, NC_WRITE, &ncid);
        if (status != NC_NOERR)
        {
            printf("Error in %s: No file \"%s\" found. No data is save. Quitting program! \n ",
                    __func__, file);
        }
        nc_redef(ncid);
    }

    nc_def_dim(ncid, "x", Nx, &x_dimid);
    nc_def_dim(ncid, "y", My, &y_dimid);
    dimids[1] = y_dimid;
    dimids[0] = x_dimid;
    nc_def_var(ncid, dataname, nctype, 2, dimids, &varid);
    nc_enddef(ncid);
    nc_put_var(ncid, varid, &data[0]);              // Write data to NetCDF file
    nc_close(ncid);                                 // Close the NetCDF file
}



/**
 * @brief      Saves one-dimensional variables into a NetCDF file.
 *
 * @param[in]  file      File will be saved in NetCDF format
 * @param[in]  dataname  Variable name saved
 * @param      data      The dataset
 * @param[in]  nctype    NetCDF data type used
 * @param[in]  dimname   Name of dimension in netcdf
 * @param[in]  Nx        The length data
 * @param[in]  write     Option to decide write new or append
 *
 * @tparam     type      { description }
 */
template<typename type>
void SaveOutput1D(const char *file, const char *dataname, type *data, nc_type nctype, 
                  const char *dimname, int Nx, int write)
{
    int ncid, varid, t_dimid, dimids[1], status;

    // Set up NetCDF file for writing
    if (write == 0)
    {
        // Create a new NetCDF file
        status = nc_create(file, NC_CLOBBER, &ncid);
        if (status != NC_NOERR)
        {
            printf("Error in %s: No file \"%s\" found. No data is save. Quitting program! \n ",
                    __func__, file);
        }

    } else {
        // Open and re-define an existing NetCDF file
        status = nc_open(file, NC_WRITE, &ncid);
        if (status != NC_NOERR)
        {
            printf("Error in %s: No file \"%s\" found. No data is save. Quitting program! \n ",
                    __func__, file);
        }
        nc_redef(ncid);
    }

    nc_def_dim(ncid, dimname, Nx, &t_dimid);
    dimids[0] = t_dimid;
    nc_def_var(ncid, dataname, nctype, 1, dimids, &varid);
    nc_enddef(ncid);
    nc_put_var(ncid, varid, &data[0]);              // Write data to NetCDF file
    nc_close(ncid);                                 // Close the NetCDF file
}


/**
 * @brief      Copy to host and save 2D and 3D results to netcdf
 *
 * @param      project          Class including project info
 * @param      overland_host    Overland flow class in host memory
 * @param      overland_dev     Overland flow class in device memory
 * @param      subsurface_host  Subsurface flow class in host memory
 * @param      subsurface_dev   Subsurface flow class in device memory
 * @param[in]  globsize         Size of the global domain
 * @param[in]  timestep         Current time step
 */
void SaveModelResults(ProjectClass *project, OverlandFlowClass *overland_host,
                      OverlandFlowClass *overland_dev, SubsurfaceFlowClass *subsurface_host,
                      SubsurfaceFlowClass *subsurface_dev, int3 globsize, int timestep)
{
    char file2D[64], file3D[64];
    int time_save = project->saveinterval;
    int sizex = globsize.x;
    int sizey = globsize.y;
    int sizez = globsize.z;
    int sizexy = sizex * sizey;
    int sizexyz = sizex * sizey * sizez;

    if ( timestep % time_save == 0 )
    {
        // 2D overland flow
        if (project->saveolf)
        {
            SafeCudaCall( cudaMemcpy(overland_host->waterelev, overland_dev->waterelev,
                                     sizexy*sizeof(double), cudaMemcpyDeviceToHost) );

            SafeCudaCall( cudaMemcpy(overland_host->waterdepth, overland_dev->waterdepth,
                                     sizexy*sizeof(double), cudaMemcpyDeviceToHost) );

            snprintf(file2D, sizeof(char) * 64, "%s/%s_%d.nc", project->folderoutput,
                                                               project->olfoutput,
                                                               timestep);

            SaveOutput2D(file2D, "water_elevation", overland_host->waterelev, NC_DOUBLE, 
                         sizex, sizey, 0);
            SaveOutput2D(file2D, "water_depth", overland_host->waterdepth, NC_DOUBLE, 
                         sizex, sizey, 1);
        }

        // 3D subsurface flow
        if (project->savessf)
        {
            SafeCudaCall( cudaMemcpy(subsurface_host->psiout, subsurface_dev->psin,
                                     sizexyz*sizeof(double), cudaMemcpyDeviceToHost) );

            SafeCudaCall( cudaMemcpy(subsurface_host->thetaout, subsurface_dev->thetan,
                                     sizexyz*sizeof(double), cudaMemcpyDeviceToHost) );

            snprintf(file3D, sizeof(char) * 64, "%s/%s_%d.nc", project->folderoutput,
                                                               project->ssfoutput, timestep);

            SaveOutput3D(file3D, "pressure_head", subsurface_host->psiout, NC_DOUBLE, globsize, 0);
            SaveOutput3D(file3D, "moisture", subsurface_host->thetaout, NC_DOUBLE, globsize, 1);
        }
    }
}


/**
 * @brief      Save MLCan results to netCDF for all processes
 *
 * @param      project    Class including project info
 * @param      outmlcan   Class including output for MLCan model
 * @param      canopies   Class including canopy variables
 * @param[in]  num_steps  The number of timesteps for simulation
 * @param[in]  rank       Global rank of the current MPI process
 * @param[in]  procsize   Total number of MPI processes available
 */
void SaveMLCanResults(ProjectClass *project, OutputClass *outmlcan, CanopyClass *canopies,
                      int num_steps, int rank, int procsize)
{
    char file_mlcan[64], dim1Dname[64];
    int nl_can = canopies->nl_can;
    snprintf(file_mlcan, sizeof(char) * 64, "%s/proc%d_%s_2d.nc", project->folderoutput, rank,
                                                                  project->mlcanoutput);

    SaveOutput2D(file_mlcan, "An_sun", outmlcan->An_sun, NC_DOUBLE, nl_can, num_steps, 0);
    SaveOutput2D(file_mlcan, "LE_sun", outmlcan->LE_sun, NC_DOUBLE, nl_can, num_steps, 1);
    SaveOutput2D(file_mlcan, "H_sun", outmlcan->H_sun, NC_DOUBLE, nl_can, num_steps, 1);
    SaveOutput2D(file_mlcan, "TR_sun", outmlcan->TR_sun, NC_DOUBLE, nl_can, num_steps, 1);
    SaveOutput2D(file_mlcan, "LAI_sun", outmlcan->LAI_sun, NC_DOUBLE, nl_can, num_steps, 1);
    SaveOutput2D(file_mlcan, "An_shade", outmlcan->An_shade, NC_DOUBLE, nl_can, num_steps, 1);
    SaveOutput2D(file_mlcan, "LE_shade", outmlcan->LE_shade, NC_DOUBLE, nl_can, num_steps, 1);
    SaveOutput2D(file_mlcan, "H_shade", outmlcan->H_shade, NC_DOUBLE, nl_can, num_steps, 1);
    SaveOutput2D(file_mlcan, "TR_shade", outmlcan->TR_shade, NC_DOUBLE, nl_can, num_steps, 1);
    SaveOutput2D(file_mlcan, "LAI_shade", outmlcan->LAI_shade, NC_DOUBLE, nl_can, num_steps, 1);

    snprintf(file_mlcan, sizeof(char) * 64, "%s/proc%d_%s_1d.nc", project->folderoutput,
                                                                  rank,
                                                                  project->mlcanoutput);

    snprintf(dim1Dname, sizeof(char) * 64, "time_series");
    SaveOutput1D(file_mlcan, "An_can", outmlcan->An_can, NC_DOUBLE, dim1Dname, num_steps, 0);
    SaveOutput1D(file_mlcan, "LE_can", outmlcan->LE_can, NC_DOUBLE, dim1Dname, num_steps, 1);
    SaveOutput1D(file_mlcan, "H_can", outmlcan->H_can, NC_DOUBLE, dim1Dname, num_steps, 1);
    SaveOutput1D(file_mlcan, "TR_can", outmlcan->TR_can, NC_DOUBLE, dim1Dname, num_steps, 1);
    SaveOutput1D(file_mlcan, "Rnrad_can", outmlcan->Rnrad_can, NC_DOUBLE, dim1Dname, num_steps, 1);
    SaveOutput1D(file_mlcan, "mbw_can", outmlcan->mbw_can, NC_DOUBLE, dim1Dname, num_steps, 1);
    
}


/**
 * @brief      Save mean results for entire simulation period
 *
 * @param      project          Class including project info
 * @param      canopies         Class including canopy variables
 * @param      subsurface_host  Subsurface flow class in host memory
 * @param      outmlcan         Class including output for MLCan model
 * @param[in]  rank             Global rank of the current MPI process
 * @param[in]  procsize         Total number of MPI processes available
 * @param[in]  globsize         Size of the global domain
 * @param      cartComm         The cartesian communications
 */
void  SaveResultEntirePeriod(ProjectClass *project, CanopyClass *canopies,
                             SubsurfaceFlowClass *subsurface_host, OutputClass *outmlcan,
                             int rank, int procsize, int3 globsize, MPI_Comm *cartComm)
{
    int isroot = rank == MPI_MASTER_RANK;
    int num_steps = project->num_steps;
    int sizez = globsize.z;
    char file1D[64], file2D[64], dim1Dname[64];

    // Save MLCan model on all processes
    if (project->savemlcan)
        SaveMLCanResults(project, outmlcan, canopies, num_steps, rank, procsize);

    // Save Flow model on MASTER only
    if (isroot)
    {
        if (project->savestat)
        {
            snprintf(file1D, sizeof(char) * 64, "%s/%s.nc", project->folderoutput,
                                                            project->ssf1doutput);
            snprintf(dim1Dname, sizeof(char) * 64, "time_series");
            SaveOutput1D(file1D, "mb_subsurfaceW", subsurface_host->mb_subsurfaceW, NC_DOUBLE, dim1Dname, num_steps, 0);


            snprintf(file2D, sizeof(char) * 64, "%s/%s.nc", project->folderoutput,
                                                            project->statoutput);
            SaveOutput2D(file2D, "pressure_mean", subsurface_host->psi_col, NC_DOUBLE, sizez,
                         num_steps, 0);
            SaveOutput2D(file2D, "moisture_mean", subsurface_host->theta_col, NC_DOUBLE, sizez,
                         num_steps, 1);
        }
    }

    MPI_Barrier(*cartComm); // sync all process
                            //
}