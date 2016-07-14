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
#include "../include/class.h"


/*
 * Main application entry point
 * ----------------------------
 *      argc   : The number of command-line arguments
 *      argv   : The list of command-line arguments
 *      return : Returns zero
 */
int main(int argc, char ** argv)
{
    int rank, procsize;                         // Rank and size of all processes (WORLD)
    const char *file_config;                    // configuration file
    mpiClass *mpiobj = new mpiClass;            // user-defined class for message passing interface
    ProjectClass *project = new ProjectClass;   // user-defined class for project info
    FileNameClass *files = new FileNameClass;   // user-defined class for file info


    // Initialize MPI parallel environment.
    InitializeMPI(&argc, &argv, &rank, &procsize);


    // Always check and parse configuration from command line for running model.
    // This provide flexibility for users to modify parameters and configs.
    ParsingCommandsAndConfigurations(argc, argv, file_config, rank, procsize, project, 
                                     files, mpiobj);

    
    // Partition MPI process to the entire domain. 2D cartesian communicator is used for domain
    // decomposition. Each process will run simulation in an averaged sub-domain.
    SetTopologyNetwork(rank, procsize, project, mpiobj);
    
   
    // Model kernel functions
    //     - MLCan initialization is in MPI host memory 
    //     - Flow model initialization is in device memory (CUDA)
    //       Note: The device is control by Master process only
    NumericalModelKernels(rank, procsize, files, project, mpiobj);


    // Close MPI parallel environment and free communicators.
    FinalizeMPI(rank, procsize, mpiobj);

    return 0;
}

