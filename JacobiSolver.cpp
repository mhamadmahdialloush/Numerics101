#include "JacobiSolver.h"


JacobiSolver::JacobiSolver(Matrix& A, Vector& x, Vector& b)
	: Solver(A, x, b, linear_solver::jacobi)
{

}

JacobiSolver::JacobiSolver(Matrix& A, Vector& x, Vector& b, int min_iters, int max_iters, double tolerance)
	: Solver(A, x, b, linear_solver::jacobi, min_iters, max_iters, tolerance, 1)
{

}


double JacobiSolver::solve()
{
	cout << "Jacobi solver" << endl;

	int n = x_.size();

	init_residual_ = calc_res();

	int iter;
	for (iter = 0; iter < max_iters_; iter++)
	{
		Vector x_new(n);

		for (int i = 0; i < n; i++)
		{
			x_new[i] = b_[i];
			for (int j = 0; j < n; j++)
			{
				if (j != i)
				{
					x_new[i] = x_new[i] - A_[i][j] * x_[j];
				}
			}
			x_new[i] = x_new[i] / A_[i][i];
		}

		x_ = x_new;

		final_residual_ = calc_res();

		cout << "  " << "iter: " << iter << ", error: " << final_residual_ << endl;

		if (final_residual_ < tolerance_ && iter>min_iters_)
		{
			break;
		}
	}

	return final_residual_;
}
