/*
 *  File: petscSolverSerBlockCSR.cpp
 *  Author: Mhamad Mahdi Alloush, PhD
 *  Affiliation: Computational Mechanics Lab, American University of Beirut
 *  Description: 
 *      This example demonstrates solving a sparse block CSR matrix using PETSc.
 *      The sparse 3x3 block matrix consists of blocks of size 2, represented in 
 *      the following format:
 *      
 *       ----            -----
 *      |3   0|  0  0   |1   0|
 *      |0  -2|  0  0   |0  -1|
 *       ----            ------
 *              ------               
 *       0  0   |5  0|   0   0
 *       0  0   |0  4|   0   0
 *              ------
 *       ----    ----    -----
 *      |1  0|  |1  0|  |3   0|      
 *      |0 -3|  |0 -1|  |0   2|
 *       ----    ----    -----
 *      
 *      The right-hand side (RHS) vector is:
 *      
 *       -1
 *        0.5
 *       -0.5
 *        1
 *        1  
 *        0.5
 *      
 *      The expected solution vector is:
 *      
 *      -0.5125
 *      -0.2500
 *      -0.1000
 *       0.2500
 *       0.5375
 *      -0.0000
 */

// Include built-in libraries
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <mpi.h>

// PETSc
#include <petscsf.h>

// Include all required non built-in header files 
#include "Types.h" // to print vector
#include "PetscSolver.h"


int main(int argc, char **argv)
{       
    // Initialize PETSc
    PetscInitialize(&argc, &argv, NULL, NULL);

    PetscInt i[] = {0, 2, 3, 6}; // Indices into j for the start of each block row
    PetscInt j[] = {0, 2, 1, 0, 1, 2}; // Column indices for each block row
    PetscScalar v[] = {3,0,0,-2,1,0,0,-1,5,0,0,4,1,0,0,-3,1,0,0,-1,3,0,0,2}; // values in the matrix (block-row-wise)

    std::vector<PetscScalar> rhs = {-1,0.5,-0.5,1,1,0.5};
    std::vector<PetscScalar> sol = {0,0,0,0,0,0}; // initially 0

    int n_rows = 3; // number of block rows
    int n_cols = 3; // number of block columns

    // solve
    {
        std::cout << "Solving with PETSc ..." << std::endl;
        
	    PetscSolver<2> solver(n_rows,n_cols,i,j,v);
	    solver.solve(sol.data(), rhs.data());

        // perform an update
        {             
            // PetscInt row_index = 0;
            // PetscInt n_cols = 2;
            // PetscInt col_indices[] = {0,2};
            // PetscScalar values[] = {8,3,2,7,10,11,12,13};

            // solver.update_row(row_index, n_cols, col_indices, values);   
        }

        // write matrix to file
        solver.write_matrix("matrix.txt", true);
    }

    // print solution
    std::cout << "Done solving the system with PETSC" << std::endl;
    std::cout << "The solution is: " << std::endl;
    std::cout << sol << std::endl;

    // Finalize PETSc
    PetscFinalize();

    return 0;
}

