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

 

///////////////////////////////////////////////
// MPI initialization and finalization       //
///////////////////////////////////////////////


/**
 * @brief      Initialize MPI processes.
 *
 * @param      argc  The number of command-line arguments
 * @param      argv  The list of command-line arguments
 * @param      rank  The global rank of the current MPI process
 * @param      size  The total number of MPI processes available
 */
void InitializeMPI( int *argc, char ***argv, int *rank, int *procsize )
{
    // Get the rank and size in the original communicator
    MPI_Init(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, rank);
    MPI_Comm_size(MPI_COMM_WORLD, procsize);
    MPI_Barrier(MPI_COMM_WORLD);
}


/**
 * @brief      : Close (finalize) the MPI environment
 *
 * @param[in]  rank      The global rank of the current MPI process
 * @param[in]  procsize  The total number of MPI processes available
 * @param      mpiobj    message passing interface object class
 */
void FinalizeMPI(int rank, int procsize, mpiClass *mpiobj)
{
    if (rank == MPI_MASTER_RANK)
    {
        printf("\n");
        printf("FINALIZING MPI PARALLEL ENVIRONMENT . . . . . . . . . . . . . . completed! \n");
        printf("----------------------------------- \n");
        printf("\n");
    }

    MPI_Barrier(mpiobj->cartComm);              // Synchronize all MPIs
    MPI_Comm_free(&mpiobj->cartComm);

    MPI_Finalize();                     // Finalize all MPIs

}



///////////////////////////////////////////////
// Topology and domain decomposition         //
///////////////////////////////////////////////


/**
 * @brief      Partitioning the entire domain to processes
 *
 * @param[in]  rank       The global rank of the current MPI process
 * @param[in]  procsize   The total number of MPI processes available
 * @param      colx       Topology size in x (horizontal) direction
 * @param      xoffset    Offset in x-direction
 * @param      xdispls    Displacements in x-direction
 * @param      xcounts    Counts in x-direction
 * @param      rowy       Topology size in y (vertical) direction
 * @param      yoffset    Offset in y-direction
 * @param      ydispls    Displacements in y-direction
 * @param      ycounts    Counts in y-direction
 * @param[in]  topolsize  Size of topology nodes (mpi)
 * @param      globsize   Size of the global domain
 * @param      cartComm   Carthesian MPI communicator
 */
void PartitionDomainToProcess(int rank, int procsize, int *colx, int *xoffset, int *xdispls,
                              int *xcounts, int *rowy, int *yoffset, int *ydispls, int *ycounts,
                              int2 topolsize, int3 *globsize, MPI_Comm *cartComm)
{
    int rowy0 = 0, colx0 = 0;
    int averowy, extrarowy, avecolx, extracolx;
    int dest, source, msgtype, tag = 1;
    int sizex = globsize->x;
    int sizey = globsize->y;
    int topolx = topolsize.x;
    int xcomm_rank, xcomm_size, ycomm_rank, ycomm_size;
    MPI_Status status;
    MPI_Comm xcomm, ycomm;

    ////////////////////////////////////////////////
    // Broadcast globsize from MASTER rank to all //
    ////////////////////////////////////////////////

    MPI_Bcast(&sizex, 1, MPI_INT, 0, *cartComm);
    MPI_Bcast(&sizey, 1, MPI_INT, 0, *cartComm);
    globsize->x = sizex;
    globsize->y = sizey;

    ////////////////////////////////////////////////
    // Split communicator for smart decomposition //
    ////////////////////////////////////////////////

    int xcolor = rank / topolx;
    int ycolor = rank % topolx;

    MPI_Comm_split(*cartComm, xcolor, rank, &xcomm);
    MPI_Comm_split(*cartComm, ycolor, rank, &ycomm);

    MPI_Comm_rank(xcomm, &xcomm_rank);
    MPI_Comm_size(xcomm, &xcomm_size);

    MPI_Comm_rank(ycomm, &ycomm_rank);
    MPI_Comm_size(ycomm, &ycomm_size);

    ///////////////////////////////////////
    // Split by column (x-direction)     //
    ///////////////////////////////////////

    if (xcomm_rank == 0)
    {   // Master process
        avecolx   = sizex / xcomm_size;
        extracolx = sizex % xcomm_size;
        *xoffset  = 0;

        for (int i = 0; i < xcomm_size; i++)
        {
            *colx = (i < xcomm_size - extracolx) ? avecolx : avecolx+1;
            xcounts[i] = *colx;

            if (i == 0)
            {
                colx0 = *colx;

                // count for offset in rank 0
                *xoffset += *colx;
                xdispls[i] = 0;
            } else {
                dest = i;
                MPI_Send(xoffset, 1, MPI_INT, dest, tag, xcomm);
                MPI_Send(colx, 1, MPI_INT, dest, tag, xcomm);
                *xoffset += *colx;
                xdispls[i] = xdispls[i-1] + xcounts[i-1];
                
            }
        }
        *xoffset = 0;
        *colx = colx0;
    }
    MPI_Barrier(*cartComm);

    if (xcomm_rank != 0)
    {   // workers
        source = 0;
        msgtype = tag;
        MPI_Recv(xoffset, 1 , MPI_INT, source, msgtype, xcomm, &status);
        MPI_Recv(colx, 1, MPI_INT, source, msgtype, xcomm, &status);
    }
    MPI_Barrier(*cartComm);

    ///////////////////////////////////////
    // Split by row (y-direction)        //
    ///////////////////////////////////////

    if (ycomm_rank == 0)
    {   // Master process
        averowy   = sizey / ycomm_size;
        extrarowy = sizey % ycomm_size;
        *yoffset  = 0;

        for (int i = 0; i < ycomm_size; i++)
        {
            *rowy = (i < ycomm_size - extrarowy) ? averowy : averowy+1;
            ycounts[i] = *rowy;


            if (i == 0)
            {
                rowy0 = *rowy;

                // count for offset in rank 0
                *yoffset += *rowy;
                ydispls[i] = 0;

            } else {
                dest = i;
                MPI_Send(yoffset, 1, MPI_INT, dest, tag, ycomm);
                MPI_Send(rowy, 1, MPI_INT, dest, tag, ycomm);
                *yoffset += *rowy;
                ydispls[i] = ydispls[i-1] + ycounts[i-1];
            }
        }
        *yoffset = 0;
        *rowy = rowy0;
    }
    MPI_Barrier(*cartComm);

    if (ycomm_rank != 0)
    {   // workers
        source = 0;
        msgtype = tag;
        MPI_Recv(yoffset, 1 , MPI_INT, source, msgtype, ycomm, &status);
        MPI_Recv(rowy, 1, MPI_INT, source, msgtype, ycomm, &status);
    }
    MPI_Barrier(*cartComm);

    MPI_Comm_free(&xcomm);
    MPI_Comm_free(&ycomm);
}



/**
 * @brief      Set up topolgy network (domain decomposition) and obtained neighbors.
 *             The size of topology must be matched with the size of MPI processes used.
 *
 * @param[in]  rank        The global rank of the current MPI process
 * @param[in]  procsize    The total number of MPI processes available
 * @param      project     Class including project info
 * @param      mpiobj      message passing interface object class
 *
 * @return     function status
 */
int SetTopologyNetwork(int rank, int procsize, ProjectClass * &project, mpiClass * &mpiobj)
{
    int isroot = (rank == MPI_MASTER_RANK);
    int topox = mpiobj->topology_size.x;
    int topoy = mpiobj->topology_size.y;
    int dimsize[2] = {topoy, topox};        // dimension of topology network used
    int usePeriods[2] = {0, 0};             // non-periodic connections
    int newCoords[2];
    int oldRank = rank;
    int sizex, sizey, xoffset=0, yoffset=0;

    mpiobj->neighbors = new int[4];
    mpiobj->xdispls   = new int[topox];
    mpiobj->ydispls   = new int[topoy];
    mpiobj->xcounts   = new int[topox];
    mpiobj->ycounts   = new int[topoy];
    mpiobj->cartComm  = MPI_COMM_WORLD;

    PartitionDomainToProcess(rank, procsize, 
                             &sizex, &xoffset, 
                             mpiobj->xdispls, 
                             mpiobj->xcounts,
                             &sizey, &yoffset, 
                             mpiobj->ydispls, 
                             mpiobj->ycounts, 
                             mpiobj->topology_size, 
                             &mpiobj->global_size, 
                             &mpiobj->cartComm);

    mpiobj->domain_size.x = sizex;
    mpiobj->domain_size.y = sizey;
    mpiobj->domain_size.z = NUM_SOIL_LAYERS;
    mpiobj->offset.x = xoffset;
    mpiobj->offset.y = yoffset;
    
    // Create a carthesian communicator
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimsize, usePeriods, 1, &mpiobj->cartComm);

    // Update the rank to be relevant to the new communicator
    MPI_Comm_rank(mpiobj->cartComm, &rank);

    // Show changes in ranks if happen
    if (rank != oldRank)
        printf("Rank change: from %d to %d\n", oldRank, rank);

    // Obtain the 2D coordinates & topology in the new cartesian communicator
    MPI_Cart_coords(mpiobj->cartComm, rank, 2, newCoords);
    mpiobj->topology_index = make_int2(newCoords[1], newCoords[0]);

    // Obtain the direct neighbor ranks
    MPI_Cart_shift(mpiobj->cartComm, 1, 1, mpiobj->neighbors + DIR_LEFT, 
                                           mpiobj->neighbors + DIR_RIGHT);
    MPI_Cart_shift(mpiobj->cartComm, 0, 1, mpiobj->neighbors + DIR_TOP, 
                                           mpiobj->neighbors + DIR_BOTTOM);

    // If reaching this point, all parameters were set up successfully
    if (isroot)
    {
        printf("\n");
        printf("DOMAIN DECOMPOSITION AND FINDING NEIGHBORS . . . . . . . . . . . completed! \n");
        printf("------------------------------------------ \n");
    }
    MPI_Barrier(mpiobj->cartComm);


    // Print the index of topology
    if (PRINT_TOPOLOGY || project->verbose)
    {
        printf("\t rank: %2d\t topoindex.X: %2d\t topoindex.Y: %2d\n", rank, 
                                                                       mpiobj->topology_index.x, 
                                                                       mpiobj->topology_index.y);
    }
    MPI_Barrier(mpiobj->cartComm);

    // Print 4 neighbors around each MPI called
    if (PRINT_NEIGHBORS || project->verbose)
    {
        if (isroot)
            printf("\n");
        MPI_Barrier(mpiobj->cartComm);

        printf("\t rank: %2d\t neighbors: %2d\t %2d\t %2d\t %2d\n", rank, 
                                                                    mpiobj->neighbors[0],
                                                                    mpiobj->neighbors[1],
                                                                    mpiobj->neighbors[2],
                                                                    mpiobj->neighbors[3]);
    }
    MPI_Barrier(mpiobj->cartComm);

    return STATUS_OK;
}
