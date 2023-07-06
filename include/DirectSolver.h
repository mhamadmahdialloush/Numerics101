#ifndef DirectSolver_h
#define DirectSolver_h


#include <iostream>
#include <vector>

#include "Solver.h"


using namespace std;



class DirectSolver : public Solver
{
private:
	double solve_lu_decomposition();
	double solve_inverse_matrix();

public:
	DirectSolver(Matrix& A, Vector& x, Vector& b);
	DirectSolver(Matrix& A, Vector& x, Vector& b, direct_solver direct_solver);

	double solve() override;

};



#endif // DirectSolver_h