#include "Solver.h"


Solver::Solver(Matrix& A, Vector& x, Vector& b, linear_solver solver) : A_(A), x_(x), b_(b), solver_(solver)
{

}


Solver::Solver(Matrix& A, Vector& x, Vector& b, linear_solver solver, direct_solver direct_solver)
	: A_(A), x_(x), b_(b), solver_(solver), direct_solver_(direct_solver)
{

}


Solver::Solver(Matrix& A, Vector& x, Vector& b, linear_solver solver, int min_iters, int max_iters, double tolerance, double relaxation_factor)
	: A_(A), x_(x), b_(b), solver_(solver), min_iters_(min_iters), max_iters_(max_iters), tolerance_(tolerance), relaxation_factor_(relaxation_factor)
{

}


double Solver::solve()
{
	cout << "Should not reach here\n." << endl;
	exit(1);
	return 0.0;
}


double Solver::calc_res()
{
	int n = x_.size();

	Vector res = b_ - (A_ * x_);

	double error_rms = 0;

	for (int i = 0; i < n; i++)
	{
		error_rms += pow(res[i], 2.0);
	}

	error_rms = sqrt(error_rms / static_cast<double>(n));

	return error_rms;
}