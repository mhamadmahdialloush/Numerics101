#include "GaussSeidelSolver.h"


GaussSeidelSolver::GaussSeidelSolver(Matrix& A, Vector& x, Vector& b) 
	: Solver(A, x, b, linear_solver::gauss_seidel)
{

}

GaussSeidelSolver::GaussSeidelSolver(Matrix& A, Vector& x, Vector& b, int min_iters, int max_iters, double tolerance)
	: Solver(A, x, b, linear_solver::gauss_seidel, min_iters, max_iters, tolerance, 1)
{

}


double GaussSeidelSolver::solve()
{
	cout << "Gauss-seidel solver" << endl;

	int n = x_.size();

	init_residual_ = calc_res();

	for (int iter = 0; iter < max_iters_; iter++)
	{
		for (int j = 0; j < n; j++)
		{
			double d = b_[j];

			for (int i = 0; i < n; i++)
			{
				if (j != i)
				{
					d -= A_[j][i] * x_[i];
				}

				x_[j] = d / A_[j][j];
			}
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
