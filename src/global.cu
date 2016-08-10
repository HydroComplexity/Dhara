/*
// Copyright (C) 2016, HydroComplexity Group
// All rights reserved.
//
// Distributed Hydrologicc and Regional Analysis (DHARA) Model
// DHARA model is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
//
// DHARA model is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with the software; if not, see <http://www.gnu.org/licenses/>.
//
// Author: levuvietphong@gmail.com (Phong Le)
*/

#include "../include/main.h"
#include "../include/cusplib.h"
#include "../include/devconst.h"


/**
 * @brief      Inverse Genuchten conversion of soil moisture - pressure
 *
 * @param      theta  Soil moisture [-]
 * @param      psi    Pressure head [L]
 * @param[in]  size   Size of the domain
 */
__global__ void vanGenuchtenInverse(double *theta, double *psi, int size)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    double lambda, m;

    while ( i < size )
    {
        lambda = n - 1.0;
        m = lambda/n;

        if (theta[i] < theta_S)
            psi[i] = -(1/alpha) * pow(pow((theta_S-theta_R)/(theta[i]-theta_R), 1/m) - 1.0, 1/n) 
                                * 0.01; // [m]
        else
            psi[i] = 0;

        // Update threads if vector is long
        i += blockDim.x * gridDim.x;

    }
}



/**
 * @brief       Genuchten conversion of soil moisture - pressure
 *
 * @param      C      Specific soil moisture capacity [1/L]
 * @param      theta  Soil moisture [-]
 * @param      Ksat   Saturated hydraulic conductivity in soil [L/T]
 * @param      K      Hydraulic conductivity in soil at moisture theta [L/T]
 * @param      psi    Pressure head [L]
 * @param[in]  size   Size of the domain
 */
__global__ 
void vanGenuchten(double *C, double *theta, double *Ksat, double *K, double *psi, int3 globsize)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int sizexyz = globsize.x * globsize.y * globsize.z;
    double Se, _theta, _psi, lambda, m;

    while ( i < sizexyz )
    {
        lambda = n - 1.0;
        m = lambda/n;

        // Compute the volumetric moisture content [eqn 21]
        _psi = psi[i] * 100;
        if ( _psi < 0 )
            _theta = (theta_S - theta_R) / pow(1.0 + pow((alpha*(-_psi)),n), m) + theta_R;
        else
            _theta = theta_S;

        theta[i] = _theta;

        // Compute the effective saturation [eqn 2]
        Se = (_theta - theta_R)/(theta_S - theta_R);

        /* . . .Compute the hydraulic conductivity [eqn 8] . . .*/
        K[i] = Ksat[i] * sqrt(Se) * (1.0 - pow( 1.0-pow(Se,1.0/m), m) ) * (1.0 - pow( 1.0-pow( Se, 1.0/m), m ));

        // Compute the specific moisture storage derivative of eqn (21).
        // So we have to calculate C = d(theta)/dh. Then the unit is converted into [1/m].
        if (_psi < 0)
            C[i] = 100 * alpha * n * (1.0/n-1.0)*pow(alpha*abs(_psi), n-1)
                * (theta_R-theta_S) * pow(pow(alpha*abs(_psi), n)+1, 1.0/n-2.0);
        else
            C[i] = 0.0;

        // Update threads if vector is long
        i += blockDim.x * gridDim.x;

    }
}


/**
 * @brief      Set boundary conditions for 6 faces of the subsurface domain
 *
 * @param      west_bc    The west boundary condition
 * @param      east_bc    The east boundary condition
 * @param      south_bc   The south boundary condition
 * @param      north_bc   The north boundary condition
 * @param      top_bc     The top boundary condition
 * @param      bottom_bc  The bottom boundary condition
 * @param[in]  sizex      The sizeof domain in x-direction
 * @param[in]  sizey      The sizeof domain in y-direction
 * @param[in]  sizez      The sizeof domain in z-direction
 */
__global__ 
void SubsurfaceSetBoundaryConditionType(int *west_bc, int *east_bc, int *south_bc, int *north_bc,
                                        int *top_bc, int *bottom_bc, int3 globsize)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int sizex = globsize.x;
    int sizey = globsize.y;
    int sizez = globsize.z;

    while (tid < sizey * sizez)
    {
        west_bc[tid] = 1;
        east_bc[tid] = 1;
        __syncthreads();               // All thread must sync at this point
        tid += blockDim.x * gridDim.x; // Update threads
    }

    tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < sizex * sizez)
    {
        south_bc[tid] = 1;
        north_bc[tid] = 1;
        __syncthreads();               // All thread must sync at this point
        tid += blockDim.x * gridDim.x; // Update threads
    }

    tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < sizex * sizey)
    {
        top_bc[tid] = 1;
        bottom_bc[tid] = 1;
        __syncthreads();               // All thread must sync at this point

        tid += blockDim.x * gridDim.x; // Update threads
    }
}


/**
 * @brief      Estimate the flux entering 1st soil layer.
 *
 * @param      ph       Ponding height [L]
 * @param      hpoten   The potential pressure head [L]
 * @param      qcapa    The flux capacity of soil [L/T]
 * @param      psinp1m  Pressure head at n+1, m [L]
 * @param      knp1m    Hydraulic conductivity at n+1,m [L]
 * @param[in]  ppt      Precipitation [L]
 * @param[in]  et       Evapotranspiration [L]
 * @param[in]  sizex    Domain size in x-direction
 * @param[in]  sizey    Domain size in y-direction
 */
__global__ 
void EstimateFluxes(double *ph, double *hpoten, double *qcapa, double *psinp1m, double *knp1m,
                    double *ppt, double et, int3 globsize)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int sizexy = globsize.x * globsize.y;

    while (tid < sizexy) {
        hpoten[tid] = ph[tid] + ppt[tid] + et;
        qcapa[tid] = -knp1m[tid]*(psinp1m[tid]-hpoten[tid]-0.5*dz) / (0.5*dz);

        // Update threads if vector is long
        tid += blockDim.x * gridDim.x;
    }
}



/**
 * @brief      Identify the top boundary type of subsurface flow3D model for the entire surface
 *
 * @param      hpoten    The potential pressure head [L]
 * @param      qcapa     The flux capacity of soil [L/T]
 * @param      topbc     The top boundary condition
 * @param      topqflux  The flux go through top soil layer
 * @param      psinp1m   Pressure head at n+1, m [L]
 * @param      Psi_top   Pressure head at top layer [L]
 * @param[in]  globsize  Size of the global domain
 */
__global__ 
void IdentifyTopBoundary(double *hpoten, double *qcapa, int *topbc, double *topqflux, 
                         double *psinp1m, double *Psi_top, int3 globsize)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int sizex = globsize.x;
    int sizey = globsize.y;

    while (tid < sizex * sizey)
    {
        if (hpoten[tid] > 0.0)
        {
            if (hpoten[tid]/dt > qcapa[tid])
            {
                topbc[tid] = 0;
                Psi_top[tid] = hpoten[tid];
            } else {
                topbc[tid] = 1;
                topqflux[tid] = hpoten[tid]/dt;
            }
        } else {
            if (psinp1m[tid] > air_dry) {
                topbc[tid] = 1;
                topqflux[tid] = hpoten[tid]/dt;
            } else {
                topbc[tid] = 0;
                Psi_top[tid] = air_dry + 0.1;
            }
        }

        // Update threads if vector is long
        tid += blockDim.x * gridDim.x;
    }
}


/**
 * @brief      Set up before running flow model
 *
 * @param      project         Class including project info
 * @param      subsurface_dev  The subsurface class in device memory
 * @param[in]  rank            Global rank of the current MPI process
 * @param[in]  procsize        Total number of MPI processes available
 * @param[in]  globsize        Size of the global domain
 * @param      cartComm        Carthesian MPI communicator
 */
void PreRunningFlowModel(ProjectClass *project, SubsurfaceFlowClass * &subsurface_dev, int rank,
                         int procsize, int3 globsize, MPI_Comm *cartComm)
{

    if (rank == MPI_MASTER_RANK)
    {
        // If any saving any data is switched on, check if output folder exist
        struct stat st = {0};
        if (stat(project->folderoutput, &st) == -1) 
        {
            mkdir(project->folderoutput, 0700);
        }

        // Print out information
        printf("\n");
        printf("\nSIMULATION STARTS \n");
        printf("------------------ \n");

        // Convert pressure head (psi) to moisture (theta)
        vanGenuchten<<<TSZ,BSZ>>>(subsurface_dev->cnp1m, subsurface_dev->thetan, 
                                  subsurface_dev->ksat, subsurface_dev->knp1m, 
                                  subsurface_dev->psin, globsize );
        cudaCheckError("vanGenuchten");
    }

    MPI_Barrier(*cartComm);
}
