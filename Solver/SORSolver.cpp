#include "SORSolver.h"


SORSolver::SORSolver(Matrix& A, Vector& x, Vector& b)
	: Solver(A, x, b, linear_solver::sor)
{

}

SORSolver::SORSolver(Matrix& A, Vector& x, Vector& b, int min_iters, int max_iters, double tolerance, double relaxation_factor)
	: Solver(A, x, b, linear_solver::sor, min_iters, max_iters, tolerance, relaxation_factor)
{

}


double SORSolver::solve()
{
	cout << "SOR solver" << endl;

	int n = x_.size();

	init_residual_ = calc_res();

	int iter;
	for (iter = 0; iter < max_iters_; iter++)
	{
		Vector y(n);
		Vector x_new = x_;

		for (int i = 0; i < n; i++)
		{
			y[i] = (b_[i] / A_[i][i]);

			for (int j = 0; j < n; j++)
			{
				if (j == i)
					continue;

				y[i] = y[i] - ((A_[i][j] / A_[i][i]) * x_new[j]);
				x_new[i] = y[i];
			}

			x_[i] = (1.0 - relaxation_factor_) * x_[i] + relaxation_factor_ * x_new[i];
		}

		final_residual_ = calc_res();

		cout << "  " << "iter: " << iter << ", error: " << final_residual_ << endl;

		if (final_residual_ < tolerance_ && iter>min_iters_)
		{
			break;
		}
	}

	return final_residual_;
}
