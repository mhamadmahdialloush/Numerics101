/*
 *  File: example1.cpp
 *  Author: Mhamad Mahdi Alloush, PhD
 *  Affiliation: Computational Mechanics Lab, American University of Beirut
 *  Description: 
 *      This example demonstrates applying a matrix multiplication process
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
	// Decalre matrices and allocate size
	Matrix A1(4, 3); // 4x3 matrix
	Matrix A2(3, 4); // 3x4 matrix

	// Fill the entries of the matrices
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

	// Print matrices
	cout << "A1:" << endl << A1 << endl << endl;
	cout << "A2:" << endl << A2 << endl << endl;

	// Print the product of their multiplication
	cout << "A1xA2:" << endl << A1 * A2 << endl;

	return 0;
}