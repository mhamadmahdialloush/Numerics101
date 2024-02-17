/*
 *  File: petscSolverParBlockCSR.cpp
 *  Author: Mhamad Mahdi Alloush, PhD
 *  Affiliation: Computational Mechanics Lab, American University of Beirut
 *  Description: 
 *      This example demonstrates solving a sparse parallel block CSR matrix using PETSc.
 *      The sparse parallel matrix consists of blocks of size 2 with two partitions, represented as follows:
 *      
 *       ----            -----
 *      |3   0|  0  0   |1   0|
 *      |0  -2|  0  0   |0  -1|
 *       ----            ------
 *              ------               partition 0
 *       0  0   |5  0|   0   0
 *       0  0   |0  4|   0   0
 *              ------
 * -----------------------------------------
 *       ----    ----    -----
 *      |1  0|  |1  0|  |3   0|      partition 1
 *      |0 -3|  |0 -1|  |0   2|
 *       ----    ----    -----
 *      
 *      The right-hand side (RHS) vector is:
 *      
 *       -1
 *        0.5     partition 0
 *       -0.5
 *        1
 *       ------
 *        1       partition 1
 *       0.5
 *      
 *      The expected solution vector is:
 *      
 *       -1.19445
 *       -1.11683      partition 0
 *       -0.181542
 *        0.453856
 * -----------------------------------------
 *        1.37124       partition 1
 *       -0.872093
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

    PetscMPIInt    rank;
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    int m; // number of local block rows of paritioned matrix
    int n; // number of local block columns of paritioned matrix
    int nnz; // number of all nonzero blocks in the partitioned matrix
    const int bs = 2; // block size

    if(rank==0)
    {
    	m = 2;
    	n = 2;
    	nnz = 3;
    } 
    else if(rank==1)
    {
    	m = 1;
    	n = 1;
    	nnz = 3;
    }

    PetscInt i[m+1];
    PetscInt j[nnz];
    PetscScalar v[nnz*bs*bs];

	std::vector<PetscScalar> sol(m*bs);
    std::vector<PetscScalar> rhs(m*bs);

    if(rank==0)
    {
	    double i_values[] = {0, 2, 3};
	    std::copy(std::begin(i_values), std::end(i_values), i);

	    double j_values[] = {0, 2, 1};
	    std::copy(std::begin(j_values), std::end(j_values), j);

	    double v_values[] = {3,0,0,-2,1,0,0,-1,5,0,0,4};
	    std::copy(std::begin(v_values), std::end(v_values), v);

	    double sol_values[] = {0,0,0,0};
	    std::copy(std::begin(sol_values), std::end(sol_values), sol.data());	    

	    double rhs_values[] = {-1,0.5,-0.5,1};
	    std::copy(std::begin(rhs_values), std::end(rhs_values), rhs.data());	    	    
	    
    } 
    else if(rank==1)
    {
	    double i_values[] = {0, 3};
	    std::copy(std::begin(i_values), std::end(i_values), i);

	    double j_values[] = {0, 1, 2};
	    std::copy(std::begin(j_values), std::end(j_values), j);

	    double v_values[] = {1,0,0,-3,1,0,0,-1,3,0,0,2};
	    std::copy(std::begin(v_values), std::end(v_values), v);

	    double sol_values[] = {0,0};
	    std::copy(std::begin(sol_values), std::end(sol_values), sol.data());	    

	    double rhs_values[] = {1,0.5};
	    std::copy(std::begin(rhs_values), std::end(rhs_values), rhs.data());
    }

    // solve
    {
	    if(rank==0)
	    {
	    	std::cout << "Solving with PETSc ..." << std::endl;
	    } 

	    PetscSolver<bs> solver(m,n,i,j,v);
	    solver.solve(sol.data(), rhs.data());

	    // Write matrix and source in ascii format
	    solver.write_matrix("matrix.txt", true);
	    solver.write_rhs("source.txt", true);   
    }

    if(rank==0)
    {
    	std::cout << "Done solving the system with PETSC" << std::endl;
    	std::cout << "The solution is: " << std::endl;
    } 

    MPI_Barrier(MPI_COMM_WORLD);  

    // print solution
	std::cout << sol << std::endl;   

    // Finalize PETSc
    PetscFinalize();

    return 0;
}

