#ifndef Solver_h
#define Solver_h


#include <iostream>
#include <vector>

#include "Matrix.h"
#include "Vector.h"
#include "Types.h"


using namespace std;



class Solver
{
protected:

    Matrix& A_;
    Vector& x_;
    Vector& b_;

    double relaxation_factor_ = 0.75;

    int max_iters_ = 20;
    int min_iters_ = 3;

    double tolerance_ = 1e-9;

    linear_solver solver_ = linear_solver::jacobi;
    direct_solver direct_solver_ = direct_solver::inverse_matrix;

    double init_residual_ = -1;
    double final_residual_ = -1;

public:

    Solver(Matrix& A, Vector& x, Vector& b, linear_solver solver);
    Solver(Matrix& A, Vector& x, Vector& b, linear_solver solver, direct_solver direct_solver);
    Solver(Matrix& A, Vector& x, Vector& b, linear_solver solver, int min_iters, int max_iters, double tolerance, double relaxation_factor);

    double calc_res();

    virtual double solve();

    const double& init_residual() { return init_residual_; };
    const double& final_residual() { return final_residual_; };

};



#endif // Solver_h