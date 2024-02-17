/*
 *  File: example1.cpp
 *  Author: Mhamad Mahdi Alloush, PhD
 *  Affiliation: Computational Mechanics Lab, American University of Beirut
 *  Description: 
 *      This example demonstrates applying matrix operations like
 *      matrix-matrix multiplication
 *      matrix-vector multiplication
 *      matrix inverse
 *      matrix determinant
 *      
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

	// Print matrices
	cout << "A1:" << endl << A1 << endl << endl;
	cout << "A2:" << endl << A2 << endl << endl;
	cout << "v:" << endl << v << endl << endl;

	cout << "A1xA2:" << endl << A1 * A2 << endl;
	cout << "A1xv:" << endl << A1 * v << endl;
	cout << "Inverse:" << endl << A1.inverse() << endl;
	cout << "Determinant: " << A1.determinant() << endl;

	return 0;
}