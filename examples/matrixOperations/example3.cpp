/*
 *  File: example1.cpp
 *  Author: Mhamad Mahdi Alloush, PhD
 *  Affiliation: Computational Mechanics Lab, American University of Beirut
 *  Description: 
 *      This example demonstrates solving a linear system Ax=b with multiple
 *      linear solvers
 */

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


int main(int argc, char **argv)
{
	// Declare and allocate
	Matrix A(3, 3);
	Vector x(3);
	Vector b(3);

	// Fill
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

	return 0;
}