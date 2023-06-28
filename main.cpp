// Include built-in libraries
#include <iostream>
#include <fstream>

// Include all non built-in header files 
#include "Matrix.h"
#include "Vector.h"
#include "DirectSolver.h"
#include "GaussSeidelSolver.h"
#include "JacobiSolver.h"
#include "SORSolver.h"

using namespace std;


void TestMatrixOperations(); // Matrix operations
void TestSolvers(); // Solvers for linear systems 

int main()
{
	//TestMatrixOperations();
	TestSolvers();

	return 0;
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

	cout << "A1xA2:" << endl << A1 * A2 << endl;
	cout << "A1xv:" << endl << A1 * v << endl;
	cout << "Inverse:" << endl << A1.inverse() << endl;
	cout << "Determinant: " << A1.determinant() << endl;
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