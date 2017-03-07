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
#include "../include/devconst.h"
#include "../include/global.h"


__device__ double maxcompolf (double a, double b)
{
    return (a < b) ? b : a;
}


/**
 * @brief      { function_description }
 *
 * @param      waterelev   The waterelev
 * @param      waterdepth  The waterdepth
 * @param      ztopo       The ztopo
 * @param[in]  size        The size
 */
__global__ void WaterElevationOverland(double *waterelev, double *waterdepth, double *ztopo,
                int size)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    while (tid < size)
    {
        // Transfer data and memory . . .
        waterelev[tid] = waterdepth[tid] + ztopo[tid];

        // Update threads if vector is long . . .
        tid += blockDim.x * gridDim.x;
    }
}



/**
 * @brief      Set the up linear systems overland.
 *
 * @param      a2d        The a2 d
 * @param      rhs2d      The rhs2 d
 * @param      waterelev  The waterelev
 * @param      ztopo      The ztopo
 * @param      kwest      The kwest
 * @param      keast      The keast
 * @param      ksouth     The ksouth
 * @param      knorth     The knorth
 * @param[in]  ppt        The ppt
 * @param[in]  et         { parameter_description }
 * @param[in]  infil      The infil
 * @param[in]  globsize   The globsize
 */
__global__ void SetUpLinearSystemsOverland(double *a2d, double *rhs2d, double *waterelev,
                double *ztopo, double *kwest, double *keast, double *ksouth, double *knorth,
                double ppt, double et, double infil, int3 globsize)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int sizex = globsize.x;
    int sizey = globsize.y;
    int sizexy = sizex * sizey;
    int i, j;
    double dtodx2, dtody2, waterdepth;

    while (tid < sizexy)
    {
        i = tid % sizex;
        j = tid / sizex;
        dtodx2 = -dt / (dx*dx);
        dtody2 = -dt / (dy*dy);
        waterdepth = 0.0;

        // Setup right hand side and boundary conditions
        rhs2d[tid] = waterelev[tid] + ppt + et + infil;

        // Set the boundary condition
        // At East boundary
        if (i == sizex-1){
            keast[tid] = 0;
        }

        // At Wets boundary
        if (i == 0){
            kwest[tid] = 0;
        }

        // At North boundary
        if (j == sizey-1){
            knorth[tid] = 0;
        }

        // At South boundary
        if (j == 0){
            ksouth[tid] = 0;
        }

        // Setup Block matrices
        a2d[0*sizexy + tid] = dtody2 * ksouth[tid];
        a2d[1*sizexy + tid] = dtodx2 * kwest[tid];
        a2d[2*sizexy + tid] = 1.0 - dtodx2 * (keast[tid] + kwest[tid]) - dtody2 * (knorth[tid] + ksouth[tid]);
        a2d[3*sizexy + tid] = dtodx2 * keast[tid];
        a2d[4*sizexy + tid] = dtody2 * knorth[tid];

        ///////////////////////////////////////////////////////
        // We need to modify matrix A & rhs2d for BCs        //
        ///////////////////////////////////////////////////////

        // Inner bound (i, j > 0 and i or j equal 1)
        // if (j == 1)
        // {
        //     a2d[0*sizexy + tid] = 0.0;
        //     rhs2d[tid] += -dtody2 * ksouth[tid] * (ztopo[tid] + waterdepth);
        // }

        // if (j == sizey-2)
        // {
        //     a2d[4*sizexy + tid] = 0.0;
        //     rhs2d[tid] += -dtody2 * knorth[tid] * (ztopo[tid] + waterdepth);
        // }

        // if (i == 1)
        // {
        //     a2d[1*sizexy + tid] = 0.0;
        //     rhs2d[tid] += -dtodx2 * kwest[tid] * (ztopo[tid] + waterdepth);
        // }

        // if (i == sizex-2)
        // {
        //     a2d[3*sizexy + tid] = 0.0;
        //     rhs2d[tid] += -dtodx2 * keast[tid] * (ztopo[tid] + waterdepth);
        // }

        // Outer bound (i, j > 0 and i or j equal 1)
        if (j == 0)
        {
            for (int k=0; k<5; k++)
                a2d[k*sizexy + tid] = 0.0;

            a2d[2*sizexy + tid] = 1.0;
            //rhs2d[tid] = ztopo[tid] + waterdepth;
            rhs2d[tid] = waterelev[tid] + waterdepth;  
        }

        if (j == sizey-1)
        {
            for (int k=0; k<5; k++)
                a2d[k*sizexy + tid] = 0.0;

            a2d[2*sizexy + tid] = 1.0;
            //rhs2d[tid] = ztopo[tid] + waterdepth;
            rhs2d[tid] = waterelev[tid] + waterdepth;
        }

        if (i == 0)
        {
            for (int k=0; k<5; k++)
                a2d[k*sizexy + tid] = 0.0;

            a2d[2*sizexy + tid] = 1.0;
            //rhs2d[tid] = ztopo[tid] + waterdepth;
            rhs2d[tid] = waterelev[tid] + waterdepth;
        }

        if (i == sizex-1)
        {
            for (int k=0; k<5; k++)
                a2d[k*sizexy + tid] = 0.0;

            a2d[2*sizexy + tid] = 1.0;
            //rhs2d[tid] = ztopo[tid] + waterdepth;
            rhs2d[tid] = waterelev[tid] + waterdepth;
        }

        __syncthreads();    // All thread must sync at this point

        tid += blockDim.x * gridDim.x;
    }
}



__global__ void UpdateKOverland(double *Hs, double *h, double *mann, double *kw, double *ke,
                double *ks, double *kn, int3 globsize )
{
    int i0 = blockIdx.x * blockDim.x + threadIdx.x;
    int j0 = blockIdx.y * blockDim.y + threadIdx.y;
    int i = i0;
    int j = j0;
    int N = globsize.x;
    int M = globsize.y;
    int I;

    double Ssw, Sse, Sss, Ssn, Kwest, Keast, Ksouth, Knorth;
    double h_, h_p1, h_m1, h_pN, h_mN;
    double mann_, mann_p1, mann_m1, mann_pN, mann_mN;
    double Hs_, Hs_p1, Hs_pN, Hs_pNp1, Hs_pNm1, Hs_m1, Hs_mN, Hs_mNp1, Hs_mNm1;

    while (i < N && j < M)
    {
        I = j * N + i;
        h_ = h[I];
        Hs_ = Hs[I];
        mann_ = mann[I];

        // At East boundary
        if (i == N-1)
        {
            h_p1 = h_;
            Hs_p1 = Hs_;
            mann_p1 = mann_;
        } else {
            h_p1 = h[I+1];
            Hs_p1 = Hs[I+1];
            mann_p1 = mann[I+1];
        }

        // At Wets boundary
        if (i == 0)
        {
            h_m1 = h_;
            Hs_m1 = Hs_;
            mann_m1 = mann_;
        } else {
            h_m1 = h[I-1];
            Hs_m1 = Hs[I-1];
            mann_m1 = mann[I-1];
        }

        // At North boundary
        if (j == M-1)
        {
            h_pN = h_;
            Hs_pN = Hs_;
            mann_pN = mann_;
        } else {
            h_pN = h[I+N];
            Hs_pN = Hs[I+N];
            mann_pN = mann[I+N];
        }

        // At South boundary
        if (j == 0)
        {
            h_mN = h_;
            Hs_mN = Hs_;
            mann_mN = mann_;
        } else {
            h_mN = h[I-N];
            Hs_mN = Hs[I-N];
            mann_mN = mann[I-N];
        }

        Hs_pNp1 = Hs[I];
        Hs_pNm1 = Hs[I];

        Hs_mNp1 = Hs[I];
        Hs_mNm1 = Hs[I];

        if (i > 0 && i < N-1 && j > 0 && j < M-1)
        {
            Hs_pNp1 = Hs[I+N+1];
            Hs_pNm1 = Hs[I+N-1];

            Hs_mNp1 = Hs[I-N+1];
            Hs_mNm1 = Hs[I-N-1];
        }

        // Calculate Keast
        Sse = sqrt( pow((Hs_p1 - Hs_)/dx, 2.) + pow((Hs_pNp1+Hs_pN-Hs_mNp1-Hs_mN)/(4*dy), 2.) );

        if (Hs_p1 > Hs_)
        {
            if (h_p1 > hmin && (h_p1+h_)/2. > hmin && abs(Sse) > delta)
            {
                Keast = pow((h_p1+h_)/2., 5./3.)/(0.5*(mann_p1+mann_)*sqrt(Sse));
            } else {
                Keast = K0;
            }
        } else {
            if (h_ > hmin && (h_p1+h_)/2. > hmin && abs(Sse) > delta)
            {
                Keast = pow((h_p1 + h_)/2., 5./3.)/(0.5*(mann_p1+mann_)*sqrt(Sse));
            } else {
                Keast = K0;
            }
        }
        ke[I] = Keast;

        // Calculate Kwest
        Ssw = sqrt( pow((Hs_ - Hs_m1)/dx, 2.) + pow((Hs_pN+Hs_pNm1-Hs_mN-Hs_mNm1)/(4*dy), 2.) );
        if (Hs_m1 > Hs_)
        {
            if (h_m1>hmin && (h_m1+h_)/2. > hmin && abs(Ssw) > delta)
            {
                Kwest = pow((h_m1+h_)/2., 5./3.)/(0.5*(mann_m1+mann_)*sqrt(Ssw));
            } else {
                Kwest = K0;
            }
        } else {
            if (h_>hmin && (h_m1+h_)/2. > hmin && abs(Ssw) > delta)
            {
                Kwest = pow((h_m1+h_)/2., 5./3.)/(0.5*(mann_m1+mann_)*sqrt(Ssw));
            } else {
                Kwest = K0;
            }
        }
        kw[I] = Kwest;

        // Calculate Knorth
        Ssn = sqrt( pow((Hs_pN - Hs_)/dy, 2.) + pow((Hs_pNp1+Hs_p1-Hs_pNm1-Hs_m1)/(4*dx), 2.) );

        if (Hs_pN > Hs_)
        {
            if (h_pN>hmin && (h_pN+h_)/2. > hmin && abs(Ssn) > delta)
            {
                Knorth = pow((h_pN+h_)/2., 5./3.)/(0.5*(mann_pN+mann_)*sqrt(Ssn));
            } else {
                Knorth = K0;
            }
        } else {
            if (h_>hmin && (h_pN+h_)/2. > hmin && abs(Ssn) > delta)
            {
                Knorth = pow((h_pN+h_)/2., 5./3.)/(0.5*(mann_pN+mann_)*sqrt(Ssn));
            } else {
                Knorth = K0;
            }
        }
        kn[I] = Knorth;

        // Calculate Ksouth
        Sss = sqrt( pow((Hs_ - Hs_mN)/dy, 2.) + pow((Hs_p1+Hs_mNp1-Hs_m1-Hs_mNm1)/(4*dx), 2.) );
        if (Hs_mN > Hs_)
        {
            if (h_mN>hmin && (h_mN+h_)/2. > hmin && abs(Sss) > delta)
            {
                Ksouth = pow((h_mN+h_)/2., 5./3.)/(0.5*(mann_mN+mann_)*sqrt(Sss));
            } else {
                Ksouth = K0;
            }
        } else {
            if (h_>hmin && (h_mN+h_)/2. > hmin && abs(Sss) > delta)
            {
                Ksouth = pow((h_mN+h_)/2., 5./3.)/(0.5*(mann_mN+mann_)*sqrt(Sss));
            } else {
                Ksouth = K0;
            }
        }
        ks[I] = Ksouth;

        __syncthreads();    // All thread must sync at this point
        i += gridDim.x * blockDim.x;
        if (i >= N) {
            j += gridDim.y * blockDim.y;
            i = i0;
        }
    }
}


__global__ void WaterDepthOverland(double *waterelev, double *waterdepth, double *ztopo, int size)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    
    while (tid < size)
    {
        // Transfer data and memory . . .
        waterdepth[tid] = maxcompolf(waterelev[tid] - ztopo[tid], 0.0);

        // Update threads if vector is long
        tid += blockDim.x * gridDim.x;
    }
}



void OverlandFlowModel(TimeForcingClass * &timeforcings, OverlandFlowClass * &overland_dev,
                       SubsurfaceFlowClass * &subsurface_dev, cuspdev_diamat &a2d_cusp,
                       cuspdev_1d &we_out_cusp, cuspdev_1d &rhs2d_cusp, cuspdev_idoper &id2d,
                       int rank, int procsize, int3 globsize, int t, int num_steps)
{
    int sizexy  = globsize.x * globsize.y;

    WaterElevationOverland<<<TSZ,BSZ>>>(overland_dev->waterelev, overland_dev->waterdepth,
                          overland_dev->ztopo, sizexy);
    cudaCheckError("WaterElevationOverland");

    UpdateKOverland<<<TSZ, BSZ>>>(overland_dev->waterelev, overland_dev->waterdepth, 
                   overland_dev->mann, overland_dev->kw, overland_dev->ke, overland_dev->ks,
                   overland_dev->kn, globsize);
    cudaCheckError("UpdateKOverland");

    SetUpLinearSystemsOverland<<<TSZ, BSZ>>>(overland_dev->a2d, overland_dev->rhs2d, 
                              overland_dev->waterelev, overland_dev->ztopo, overland_dev->kw,
                              overland_dev->ke, overland_dev->ks, overland_dev->kn, 
                              0., 0., 0., globsize);
    cudaCheckError("SetUpLinearSystemsOverland");

    // Set A, b, and boundary conditions
    cusp::monitor <double> monitor(rhs2d_cusp, 100, 1e-9);
    cusp::krylov::cg(a2d_cusp, we_out_cusp, rhs2d_cusp, monitor, id2d);

    WaterDepthOverland<<<TSZ,BSZ>>>(overland_dev->waterelev, overland_dev->waterdepth,
                                    overland_dev->ztopo, sizexy);
    cudaCheckError("WaterDepthOverland");

    SafeCudaCall( cudaMemcpy(overland_dev->ph, overland_dev->waterdepth, sizexy*sizeof(double),
            cudaMemcpyDeviceToDevice) );
}
