// Include built-in libraries
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <mpi.h>

// Include all non built-in header files 
#include "Matrix.h"
#include "Vector.h"
#include "DirectSolver.h"
#include "GaussSeidelSolver.h"
#include "JacobiSolver.h"
#include "SORSolver.h"
#include "CSRMatrix.h"
#include "PetscSolver.h"

using namespace std;


void TestParBlockCRSPetscSolver(int argc, char **argv);
void TestSerBlockCRSPetscSolver(int argc, char **argv);

void TestMatrixOperations(); // Matrix operations
void TestSparseMatrixOperations();
void TestSolvers(); // Solvers for linear systems 
void TestSSHeatEquation();
void TestMatrixMultiplication();

int main(int argc, char **argv)
{
	TestParBlockCRSPetscSolver(argc, argv);
	//TestMatrixOperations();
	//TestSparseMatrixOperations();
	//TestSolvers();
	//TestSSHeatEquation();
	//TestMatrixMultiplication();
	return 0;
}

void TestMatrixMultiplication()
{
	Matrix A1(4, 3);
	Matrix A2(3, 4);

	A1[0][0] = 11.0;
	A1[0][1] = 22.0;
	A1[0][2] = 3.0;

	A1[1][0] = 1.0;
	A1[1][1] = 7.0;
	A1[1][2] = 3.0;

	A1[2][0] = 1.0;
	A1[2][1] = 2.0;
	A1[2][2] = 20.0;

	A1[3][0] = 4.0;
	A1[3][1] = 4.0;
	A1[3][2] = 30.0;

	A2[0][0] = 1.0;
	A2[0][1] = 2.0;
	A2[0][2] = 3.0;
	A2[0][3] = 3.0;

	A2[1][0] = 8.0;
	A2[1][1] = 2.0;
	A2[1][2] = 3.0;
	A2[1][3] = 2.0;

	A2[2][0] = 1.0;
	A2[2][1] = 2.0;
	A2[2][2] = 3.0;
	A2[2][3] = 1.0;

	cout << "A1xA2:" << endl << A1 * A2 << endl;
}


void TestMatrixOperations()
{
	Matrix A1(4, 4);
	Matrix A2(4, 4);
	Vector v(4);

	A1[0][0] = 1.0;
	A1[0][1] = 2.0;
	A1[0][2] = 3.0;
	A1[0][3] = 1.0;

	A1[1][0] = 1.0;
	A1[1][1] = 7.0;
	A1[1][2] = 3.0;
	A1[1][3] = 2.0;

	A1[2][0] = 1.0;
	A1[2][1] = 2.0;
	A1[2][2] = 20.0;
	A1[2][3] = 3.0;

	A1[3][0] = 4.0;
	A1[3][1] = 4.0;
	A1[3][2] = 30.0;
	A1[3][3] = 4.0;

	A2[0][0] = 1.0;
	A2[0][1] = 2.0;
	A2[0][2] = 3.0;
	A2[0][3] = 3.0;

	A2[1][0] = 8.0;
	A2[1][1] = 2.0;
	A2[1][2] = 3.0;
	A2[1][3] = 2.0;

	A2[2][0] = 1.0;
	A2[2][1] = 2.0;
	A2[2][2] = 3.0;
	A2[2][3] = 1.0;

	A2[3][0] = 9.0;
	A2[3][1] = 4.0;
	A2[3][2] = 3.0;
	A2[3][3] = 4.0;

	v[0] = 5;
	v[1] = 33;
	v[2] = 10;
	v[3] = 6;

	int s1 = 4;

	cout << "A1xA2:" << endl << A1 * A2 << endl;
	cout << "A1xv:" << endl << A1 * v << endl;
	cout << "Inverse:" << endl << A1.inverse() << endl;
	cout << "Determinant: " << A1.determinant() << endl;
	cout << "A1 x " << s1 << " = " << endl << A1 * s1 << endl;
}

void TestSparseMatrixOperations()
{
	Matrix A
	(
		{
			{1,0,0,2},
			{5,0,3,1},
			{6,8,0,0},
			{0,0,0,4}
		}
	);

	cout << "A:" << endl << A << endl;

	CSRMatrix Asparse({ 1, 2, 5, 3, 1, 6, 8, 4 }, { 0, 2, 5, 7, 8 }, { 0, 3, 0, 2, 3, 0, 1, 3 });

	cout << "Asparse:" << endl << Asparse << endl;

	Vector v({ 1,2,3,4 });

	cout << "A x v:" << endl << A * v << endl;
	cout << "Asparse x v:" << endl << Asparse * v << endl;
}

void TestSolvers()
{
	Matrix A(3, 3);
	Vector x(3);
	Vector b(3);

	A[0][0] = 4;
	A[0][1] = 1;
	A[0][2] = 2;

	A[1][0] = 3;
	A[1][1] = 5;
	A[1][2] = 1;

	A[2][0] = 1;
	A[2][1] = 1;
	A[2][2] = 3;

	b[0] = 4;
	b[1] = 7;
	b[2] = 3;

	// Solve with LU-decomposition method
	{
		DirectSolver dsolver(A, x, b, direct_solver::lu_decomposition);
		dsolver.solve();
		cout << endl << "x:" << endl << x << endl;
	}

	// Solve with inverse matrix method
	{
		x = 0; // reset vector to 0
		DirectSolver dsolver(A, x, b, direct_solver::inverse_matrix);
		dsolver.solve();
		cout << endl << "x:" << endl << x << endl;
	}

	// Solve with Gauss-Seidel iterative method
	{
		x = 0; // reset to 0
		GaussSeidelSolver gssolver(A, x, b);
		gssolver.solve();
		cout << endl << "x:" << endl << x << endl;
	}

	// Solve with Jacobi iteration method
	{
		x = 0; // reset to 0
		JacobiSolver jsolver(A, x, b, 5, 30, 1e-12);
		jsolver.solve();
		cout << endl << "x:" << endl << x << endl;
	}
}


void TestSSHeatEquation()
{
	// Domain settings

	double L = 1; // Length of domain (x-dir)
	double W = 0.5; // Width of domain (y-dir)

	const int m = 20; // number of mesh spacings in y
	const int n = 20; // number of mesh spacings in x

	// Mesh settings

	const int N = (m + 1) * (n + 1); // total number of grid points in the mesh

	double dx = L / static_cast<double>(n);
	double dy = W / static_cast<double>(m);

	// System

	Matrix A(N, N);
	Vector x(N, 300);
	Vector b(N);

	// BC
	double T_u = 600;
	double T_d = 900;
	double T_l = 400;

	// Boundary conditions

	// Up
	for (int j = 0; j < n + 1; j++)
	{
		int i = m;
		int index = i * (n + 1) + j;
		x[index] = T_u;
	}

	// Down
	for (int j = 0; j < n + 1; j++)
	{
		int i = 0;
		int index = i * (n + 1) + j;
		x[index] = T_d;
	}

	// Left
	for (int i = 0; i < m + 1; i++)
	{
		int j = 0;
		int index = i * (n + 1) + j;
		x[index] = T_l;
	}

	for (int i = 0; i < m + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			if (i == 0 || i == m || j == 0)
			{
				// Dirichlet BC up, down and left
				int index = i * (n + 1) + j;
				A[index][index] = 1;
				b[index] = x[index];
			}
			else if (j == n)
			{
				// Zero-gradient BC at the right side. Assume first-order interpolation
				int index = i * (n + 1) + j;
				int index_l = i * (n + 1) + j - 1;

				A[index][index] = 1;
				A[index][index_l] = -1;
				b[index] = 0;
			}
			else
			{
				int index = i * (n + 1) + j;

				int index_l = i * (n + 1) + j - 1;
				int index_r = i * (n + 1) + j + 1;
				int index_d = (i - 1) * (n + 1) + j;
				int index_u = (i + 1) * (n + 1) + j;

				A[index][index] += 2.0 * (1.0 / pow(dx, 2.0) + 1.0 / pow(dy, 2.0));
				A[index][index_l] -= 1.0 / pow(dx, 2.0);
				A[index][index_r] -= 1.0 / pow(dx, 2.0);
				A[index][index_d] -= 1.0 / pow(dy, 2.0);
				A[index][index_u] -= 1.0 / pow(dy, 2.0);
			}
		}
	}

	GaussSeidelSolver solver(A, x, b, 10, 1000, 1e-5);
	//DirectSolver solver(A, x, b, direct_solver::lu_decomposition);
	solver.solve();


	int I = m / 2;
    Vector sample(n + 1);
    for (int J = 0; J < n + 1; J++)
    {
        int index = I * (n + 1) + J;
        sample[J] = x[index];
    }

    ofstream out("out.txt");
    for (int i = 0; i < sample.size(); i++)
    {
        out << sample[i] << endl;
    }

    out.close();

}


void TestSerBlockCRSPetscSolver(int argc, char **argv)
{
    // This is the sparse parallel matrix:
    
    //  ----            -----
    // |3   0|  0  0   |1   0|
    // |0  -2|  0  0   |0  -1|
    //  ----            ------
    //         ------               
    //  0  0   |5  0|   0   0
    //  0  0   |0  4|   0   0
    //         ------
    //  ----    ----    -----
    // |1  0|  |1  0|  |3   0|      
    // |0 -3|  |0 -1|  |0   2|
    //  ----    ----    -----


    // This is the rhs vector

	//  -1
	//   0.5
	//  -0.5
	//   1
	//   1  
	// 0.5


    // the expected solution:

	// -1.19445
	// -1.11683
	// -0.181542
	// 0.453856
	// 1.37124
	// -0.872093

    
    // Initialize PETSc
    PetscInitialize(&argc, &argv, NULL, NULL);

    PetscInt i[] = {0, 2, 3, 6}; // Indices into j for the start of each local row
    PetscInt j[] = {0, 2, 1, 0, 1, 2}; // Column indices for each local row (block positions)
    PetscScalar v[] = {3,0,0,-2,1,0,0,-1,5,0,0,4,1,0,0,-3,1,0,0,-1,3,0,0,2}; // Optional values in the matrix

    std::vector<PetscScalar> rhs = {-1,0.5,-0.5,1,1,0.5};
    std::vector<PetscScalar> sol = {0,0,0,0,0,0}; // initially 0

    // solve
    {
	    PetscSolver<2> solver(3,3,i,j,v);
	    solver.solve(sol.data(), rhs.data());
    }

    // print solution
    std::cout << sol << std::endl;

    // Finalize PETSc
    PetscFinalize();
}





void TestParBlockCRSPetscSolver(int argc, char **argv)
{
    // This is the sparse parallel matrix:
    
    //  ----            -----
    // |3   0|  0  0   |1   0|
    // |0  -2|  0  0   |0  -1|
    //  ----            ------
    //         ------               partition 0
    //  0  0   |5  0|   0   0
    //  0  0   |0  4|   0   0
    //         ------
	// -----------------------------------------
    //  ----    ----    -----
    // |1  0|  |1  0|  |3   0|      partition 1
    // |0 -3|  |0 -1|  |0   2|
    //  ----    ----    -----


    // This is the rhs vector

	//  -1
	//   0.5     partition 0
	//  -0.5
	//   1
	// ------
	//   1       partition 1
	// 0.5


    // the expected solution:

	// -1.19445
	// -1.11683      partition 0
	// -0.181542
	// 0.453856
	// -----------------------------------------
	// 1.37124       partition 1
	// -0.872093

    
    // Initialize PETSc
    PetscInitialize(&argc, &argv, NULL, NULL);

    PetscMPIInt    rank;
    
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    int m; // number of local block rows of paritioned matrix
    int nnz; // number of nonzero blocks in the partitioned matrix
    const int bs = 2; // block size

    if(rank==0)
    {
    	m = 2;
    	nnz = 3;
    } 
    else if(rank==1)
    {
    	m = 1;
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
	    PetscSolver<bs> solver(m,m,i,j,v);
	    solver.solve(sol.data(), rhs.data());

	    solver.write_matrix("matrix.txt", true);
	    solver.write_matrix("matrix_bin.txt", false);
	    solver.write_rhs("source.txt", true);

	    // print solution
		std::cout << sol << std::endl;      
    }

    // Finalize PETSc
    PetscFinalize();
}

