#ifndef Types_h
#define Types_h


enum linear_solver
{
	direct,
	jacobi,
	gauss_seidel,
	sor,
	ilu
};


enum direct_solver
{
	lu_decomposition,
	inverse_matrix
};



#endif
