#include "DirectSolver.h"


DirectSolver::DirectSolver(Matrix& A, Vector& x, Vector& b)
	: Solver(A, x, b, linear_solver::direct)
{

}

DirectSolver::DirectSolver(Matrix& A, Vector& x, Vector& b, direct_solver ds)
	: Solver(A, x, b, linear_solver::direct, ds)
{

}


double DirectSolver::solve()
{
	cout << "Direct solver" << endl;

	if (direct_solver_ == direct_solver::inverse_matrix)
	{
		return solve_inverse_matrix();
	}
	else
	{
		return solve_lu_decomposition();
	}
}


double DirectSolver::solve_lu_decomposition()
{
	int n = x_.size();

	init_residual_ = calc_res();

	cout << "\tLU decomposition method" << endl;

	Matrix L(n, n);
	Matrix U(n, n);
	A_.lu_decomposition(L, U);

	// Forward substitution
	Vector y(n);

	y[0] = b_[0] / L[0][0];
	for (int i = 1; i < n; i++)
	{
		y[i] = b_[i] / L[i][i];
		for (int j = 0; j < i; j++)
		{
			y[i] -= L[i][j] / L[i][i] * y[j];
		}
	}


	// Backward substitution
	x_[n - 1] = y[n - 1] / U[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--)
	{
		x_[i] = y[i] / U[i][i];
		for (int j = i + 1; j < n; j++)
		{
			x_[i] -= U[i][j] / U[i][i] * x_[j];
		}
	}

	final_residual_ = calc_res();

	return final_residual_;
}

double DirectSolver::solve_inverse_matrix()
{
	int n = x_.size();

	init_residual_ = calc_res();

	cout << "\tinverse matrix method" << endl;
	x_ = (A_.inverse() * b_);

	final_residual_ = calc_res();

	return final_residual_;
}
