#ifndef SORSolver_h
#define SORSolver_h


#include <iostream>
#include <vector>

#include "Solver.h"


using namespace std;



class SORSolver : public Solver
{

public:
	SORSolver(Matrix& A, Vector& x, Vector& b);
	SORSolver(Matrix& A, Vector& x, Vector& b, int min_iters, int max_iters, double tolerance, double relaxation_factor);



	double solve() override;


};



#endif // SORSolver_h