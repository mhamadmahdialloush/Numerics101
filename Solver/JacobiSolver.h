#ifndef JacobiSolver_h
#define JacobiSolver_h


#include <iostream>
#include <vector>

#include "Solver.h"


using namespace std;



class JacobiSolver : public Solver
{

public:
	JacobiSolver(Matrix& A, Vector& x, Vector& b);
	JacobiSolver(Matrix& A, Vector& x, Vector& b, int min_iters, int max_iters, double tolerance);

	double solve() override;
};



#endif // JacobiSolver_h