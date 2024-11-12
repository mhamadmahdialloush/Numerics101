#ifdef USE_PETSC

#include "PetscSolver.h"

// Constructors

template<size_t bs>
PetscSolver<bs>::PetscSolver(PetscInt m, PetscInt n, PetscInt i[], PetscInt j[], PetscScalar v[])
{
    MatCreate(PETSC_COMM_WORLD, &A_);
    MatSetSizes(A_, m*bs, n*bs, PETSC_DETERMINE, PETSC_DETERMINE);
    MatSetType(A_, MATMPIBAIJ);
    MatMPIBAIJSetPreallocationCSR(A_, bs, i, j, v);     
}


template<size_t bs>
PetscSolver<bs>::~PetscSolver()
{
    VecDestroy(&rhs_); 
    VecDestroy(&sol_);        
    KSPDestroy(&ksp_); 
    MatDestroy(&A_); 
}


// Methods

template<size_t bs>
void PetscSolver<bs>::solve(PetscScalar sol[], PetscScalar rhs[])
{
    MatCreateVecs(A_, &sol_, &rhs_);

    PetscObjectReference((PetscObject)rhs_);
    VecPlaceArray(rhs_, reinterpret_cast<double*>(const_cast<double*>(rhs))); 

    PetscObjectReference((PetscObject)sol_);
    VecPlaceArray(sol_, reinterpret_cast<double*>(const_cast<double*>(sol)));
    

    KSPCreate(PetscObjectComm((PetscObject)A_), &ksp_);
    KSPSetOperators(ksp_, A_, A_);

    KSPSetType(ksp_, KSPFGMRES); // Try a different solver type (flexible GMRES)
    KSPSetTolerances(ksp_, 0.01, 1e-9, PETSC_DEFAULT, 10);
    KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE); 

    // Solve A x = b
    KSPSolve(ksp_, rhs_, sol_);

    PetscInt nIters;
    KSPGetIterationNumber(ksp_, &nIters); 

    // Compute the residual norm
    PetscScalar res;
    KSPGetResidualNorm(ksp_, &res);

    int rank;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank==0)    
    {
        std::cout << "nIters: " << nIters << std::endl;    
        std::cout << "residual: " << res << std::endl;   
    } 
}


template<size_t bs>
void PetscSolver<bs>::update_row(PetscInt row_index, PetscInt n_cols, PetscInt col_indices[], PetscScalar vals[])
{
    for(int j=0; j<n_cols; j++)
    {
        MatSetValuesBlocked(A_, 1, &row_index, 1, &col_indices[j], &vals[bs*bs*j], INSERT_VALUES); 
    }

    // Essential after setting values: assemble the matrix
    MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);      
}


template<size_t bs>
void PetscSolver<bs>::update_values(PetscScalar vals[], PetscInt nnz)
{
    PetscScalar *array;
    MatSeqBAIJGetArray(A_, &array);   

    for(int i=0; i<nnz*bs*bs; i++)
    {
        array[i] = vals[i];
    } 
}


template<size_t bs>
void PetscSolver<bs>::write_matrix(std::string file_name, bool ascii)
{
    if(ascii)
    {
        // Print the matrix to a file in ASCII format
        PetscViewer viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
        MatView(A_, viewer);
        PetscViewerDestroy(&viewer);
    }
    else
    {
        // Print the matrix to a file in binary format
        PetscViewer viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, file_name.c_str(), FILE_MODE_WRITE, &viewer);
        MatView(A_, viewer);
        PetscViewerDestroy(&viewer);  
    }
}


template<size_t bs>
void PetscSolver<bs>::write_sol(std::string file_name, bool ascii)
{
    if(ascii)
    {
        // Print the vector to a file in ASCII format
        PetscViewer viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
        VecView(sol_, viewer);
        PetscViewerDestroy(&viewer);
    }
    else
    {
        // Print the vector to a file in binary format
        PetscViewer viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, file_name.c_str(), FILE_MODE_WRITE, &viewer);
        VecView(sol_, viewer);
        PetscViewerDestroy(&viewer);  
    }
}


template<size_t bs>
void PetscSolver<bs>::write_rhs(std::string file_name, bool ascii)
{
    if(ascii)
    {
        // Print the vector to a file in ASCII format
        PetscViewer viewer;
        PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
        VecView(rhs_, viewer);
        PetscViewerDestroy(&viewer);
    }
    else
    {
        // Print the vector to a file in binary format
        PetscViewer viewer;
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, file_name.c_str(), FILE_MODE_WRITE, &viewer);
        VecView(rhs_, viewer);
        PetscViewerDestroy(&viewer);  
    }
}

#endif // USE_PETSC