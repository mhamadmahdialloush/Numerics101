#ifdef USE_PETSC

#ifndef PetscSolver_h
#define PetscSolver_h


#include <petscmat.h>
#include <petscblaslapack.h>
#include <petscsf.h>
#include "petscmat.h"
#include "petscksp.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>


template<size_t bs>
class PetscSolver
{
protected:

	Mat A_;	

	Vec rhs_;

	Vec sol_;

    KSP ksp_;

public:

    // Constructors

    PetscSolver(PetscInt m, PetscInt n, PetscInt i[], PetscInt j[], PetscScalar v[]);


    ~PetscSolver();


    // Methods

    void solve(PetscScalar sol[], PetscScalar rhs[]);


    // Update

    void update_row(PetscInt row_index, PetscInt col_indices[], PetscScalar vals[]);


    // IO

    void write_matrix(std::string file_name, bool ascii = false);

    void write_sol(std::string file_name, bool ascii = false);

    void write_rhs(std::string file_name, bool ascii = false);

};

#include "PetscSolver.hpp"

#endif

#endif // USE_PETSC
