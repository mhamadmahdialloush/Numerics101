/*
 *  File: solveHeatEquation.cpp
 *  Author: Mhamad Mahdi Alloush, PhD
 *  Affiliation: Computational Mechanics Lab, American University of Beirut
 *  Description: 
 *      This example demonstrates solving the 2D heat equation using finite
 *		difference method. The domain is 1x0.5 in size, and the mesh involves
 *		20x20 spacings. The temperature is fixes on the top, bottom and left 
 *		sides, while is set to zero-gradient on the right side. The solver
 *		will output a file containing the solution over a line which cuts
 *		the domain in the middle from left to right
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

    return 0;
}
