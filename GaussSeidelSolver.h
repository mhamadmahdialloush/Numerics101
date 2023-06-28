#ifndef GaussSeidelSolver_h
#define GaussSeidelSolver_h


#include <iostream>
#include <vector>

#include "Solver.h"


using namespace std;



class GaussSeidelSolver: public Solver
{

public:
	GaussSeidelSolver(Matrix& A, Vector& x, Vector& b);
	GaussSeidelSolver(Matrix& A, Vector& x, Vector& b, int min_iters, int max_iters, double tolerance, double relaxation_factor);



	double solve() override;


};



#endif // GaussSeidelSolver_h