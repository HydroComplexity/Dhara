#ifndef __MAIN_H__
#define __MAIN_H__
#endif

#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <netcdf.h>
#include <Eigen/Dense>
#include "class.h"
#include "vertclass.h"
#include "saveclass.h"

using namespace std;
using namespace Eigen;


///////////////////////////////////////////
// User-define values                    //
///////////////////////////////////////////

// MPI extra definitions used              
#define MPI_MASTER_RANK             0

// Soil layers
#define NUM_SOIL_LAYERS             15

// Options for printing out                
#define PRINT_NEIGHBORS             0
#define PRINT_TOPOLOGY              0
#define PRINT_INITIALDATA           0
#define PRINT_RESULTS               0

// Indicator for surrounding neighbors     
#define DIR_TOP                     0
#define DIR_RIGHT                   1
#define DIR_BOTTOM                  2
#define DIR_LEFT                    3

// Status value that indicates a successful operation or an error
#define STATUS_OK                   0
#define STATUS_ERR                  -1

// Extra MLCan constant
#define pinum 3.1415926535897932384626433832795028841971693993751058209749445923
#define Jmax_Bisec 100

// Support functions for printing & checking neighbors
#define cudaCheckError(msg) \
  do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                            msg, cudaGetErrorString(__err), \
                            __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
      } while (0)

#define SafeCudaCall(call)          CheckCudaCall(call, #call, __FILE__, __LINE__)
#define SafeHostFree(block)         { if (block) free(block); }
#define SafeDevFree(block)          { if (block) SafeCudaCall(cudaFree(block)); }
#define OnePrintf(allow, ...)       { if (allow) printf(__VA_ARGS__); }
#define OneErrPrintf(allow, ...)    { if (allow) fprintf(stderr, __VA_ARGS__); }
#define HasNeighbor(neighbors, dir) (neighbors[dir] != MPI_PROC_NULL)

// GPU block and grid sizes
#define BSZ 256
#define TSZ 256


///////////////////////////////////////////////////
// MPI support functions                         //
///////////////////////////////////////////////////

void InitializeMPI(int *argc, char ***argv, int *rank, int *procsize);

void FinalizeMPI(int rank, int procsize, mpiClass *mpiobj);

int SetTopologyNetwork(int rank, int procsize, ProjectClass * &project, mpiClass * &mpiobj);



///////////////////////////////////////////////////
// Input and output funtions                     //
///////////////////////////////////////////////////

void ParsingCommandsAndConfigurations(int argc, char **argv, const char * &file_config, int rank,
                                      int procsize, ProjectClass * &project, 
                                      FileNameClass * &files, mpiClass * &mpiobj);

void LoadConfigMLCanModel(ProjectClass *project, FileNameClass *files, SwitchClass *switches,
                          ConstantClass *constants, CanopyClass *canopies, SoilClass *soils,
                          RadiationClass *radiation, PhotosynthesisClass *photosynthesis,
                          RespirationClass *respiration, StomaConductClass *stomaconduct,
                          MicroEnvironmentClass *microenviron);

void LoadCanopyRootDistributions(ProjectClass *project, CanopyClass *canopies, 
                                 VerticalCanopyClass *vertcanopies, SoilClass *soils,
                                 VerticalSoilClass *vertsoils);

void LoadForcingData(FileNameClass *files, CanopyClass *canopies, TimeForcingClass *timeforcings,
                     int rank, int procsize, int num_steps);

void LoadFlowModelConfig(ProjectClass *project, FileNameClass *files, OverlandFlowClass *overland,
                         SubsurfaceFlowClass *subsurface );

void LoadTopography(FileNameClass *files, OverlandFlowClass *overland, int rank, int procsize);


void GetFileInfo(const char *file_name, const char *var_name, int ndims, int *dim);


void SetFlowModelConditions(TimeForcingClass *timeforcings, FileNameClass *files, 
                            OverlandFlowClass *overland, SubsurfaceFlowClass * &subsurface, 
                            int num_steps, int rank, int procsize, int3 globsize);

void SetMPIGPUMapping(mpiClass *mpiobj, SubsurfaceFlowClass * &subsurface, int procsize,
                      int3 globsize, int2 topolsize);

void GatherRootDensityToMasterAndFilter(ProjectClass *project, VerticalSoilClass *vertsoils,
                                        SoilClass *soils, SubsurfaceFlowClass *subsurface_host,
                                        SubsurfaceFlowClass *subsurface_dev, int rank, 
                                        int procsize, MPI_Comm *cartComm);

///////////////////////////////////////////////////
// Numerical functions                           //
///////////////////////////////////////////////////

void NumericalModelKernels(int rank, int procsize, FileNameClass *files, ProjectClass *project,
                           mpiClass *mpiobj);


///////////////////////////////////////////////////
// Device function wrappers                      //
///////////////////////////////////////////////////

#ifdef __cplusplus
extern "C" 
{
#endif

void CheckCudaCall(cudaError_t command, const char * commandName, const char * fileName, int line);

void PreRunningFlowModel(ProjectClass *project, SubsurfaceFlowClass * &subsurface_dev, int rank,
                         int procsize, int3 globsize, MPI_Comm *cartComm);

void SaveModelResults(ProjectClass *project, OverlandFlowClass *overland_host, 
                      OverlandFlowClass *overland_dev, SubsurfaceFlowClass *subsurface_host,
                      SubsurfaceFlowClass *subsurface_dev, int3 globsize, int t);

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
                         MPI_Comm *cartComm);

/*
void CopyConstantToDevice(ProjectClass *project, OverlandFlowClass * &overland,
                          SubsurfaceFlowClass *subsurface);
*/
                          
void CopyDataToDevice(ProjectClass *project, TimeForcingClass *timeforcings,
                      OverlandFlowClass *overland_host, OverlandFlowClass *overland_dev,
                      SubsurfaceFlowClass *subsurface_host, SubsurfaceFlowClass *subsurface_dev,
                      int3 globsize);

void  SaveResultEntirePeriod(ProjectClass *project, CanopyClass *canopies,
                             SubsurfaceFlowClass *subsurface_host, OutputClass *outmlcan,
                             int rank, int procsize, int3 globsize, MPI_Comm *cartComm);

#ifdef __cplusplus
}   // extern "C"
#endif