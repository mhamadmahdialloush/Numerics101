#ifndef Types_h
#define Types_h


#include <mpi.h>
#include <vector>


// print out std::vector
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& data)
{
    int rank, size;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(size>1)
    {
		for(int irank=0; irank<size; irank++)
		{
			if(irank==rank)
			{
				os << "Rank " << irank << std::endl;
			    os << data.size() << std::endl << "(\n";
			    for(int i=0; i<data.size(); i++)
			    {
			        os << data[i] << std::endl;
			    }
			    os << ")\n";
			    os << "-----------------" << std::endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
    }
    else
    {
	    os << data.size() << std::endl << "(\n";
	    for(int i=0; i<data.size(); i++)
	    {
	        os << data[i] << std::endl;
	    }
	    os << ")\n";
    }

    return os;
}


enum linear_solver
{
	direct,
	jacobi,
	gauss_seidel,
	sor,
	ilu,
	petsc
};


enum direct_solver
{
	lu_decomposition,
	inverse_matrix
};



#endif
