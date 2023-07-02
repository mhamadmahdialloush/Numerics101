#include "Matrix.h"
#include "Vector.h"


Matrix::Matrix(int m, int n)
{
    data_.resize(m);
    for (auto& elem : data_)
    {
        elem.resize(n, 0);
    }
}

double Matrix::determinant()
{
	return determinant((*this), this->cols());
}

Matrix Matrix::inverse()
{
	int n = this->cols();
	Matrix inv_matrix(n, n);
	inverse((*this), inv_matrix);
	return inv_matrix;
}

void Matrix::get_cofactor(Matrix& A, Matrix& temp, int p, int q, int n)
{
	int N = A.rows();

	int i = 0, j = 0;

	// Looping for each element of the matrix
	for (int row = 0; row < n; row++) 
	{
		for (int col = 0; col < n; col++) 
		{
			// Copying into temporary matrix only those
			// element which are not in given row and
			// column
			if (row != p && col != q) 
			{
				temp[i][j++] = A[row][col];

				// Row is filled, so increase row index and
				// reset col index
				if (j == n - 1) 
				{
					j = 0;
					i++;
				}
			}
		}
	}
}

/* Recursive function for finding determinant of matrix.
n is current dimension of A[][]. */
double Matrix::determinant(Matrix& A, int n)
{
	int N = A.rows();

	double D = 0; // Initialize result

	// Base case : if matrix contains single element
	if (n == 1)
		return A[0][0];

	Matrix temp(N, N); // To store cofactors

	int sign = 1; // To store sign multiplier

	// Iterate for each element of first row
	for (int f = 0; f < n; f++) 
	{
		// Getting Cofactor of A[0][f]
		get_cofactor(A, temp, 0, f, n);
		D += sign * A[0][f] * determinant(temp, n - 1);

		// terms are to be added with alternate sign
		sign = -sign;
	}

	return D;
}

void Matrix::adjoint(Matrix& A, Matrix& adj)
{
	int N = A.rows();

	if (N == 1) 
	{
		adj[0][0] = 1;
		return;
	}

	// temp is used to store cofactors of A[][]
	int sign = 1;
	Matrix temp(N, N);

	for (int i = 0; i < N; i++) 
	{
		for (int j = 0; j < N; j++) 
		{
			// Get cofactor of A[i][j]
			get_cofactor(A, temp, i, j, N);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j][i] = (sign) * (determinant(temp, N - 1));
		}
	}
}

// Function to calculate and store inverse, returns false if
// matrix is singular
bool Matrix::inverse(Matrix& A, Matrix& inverse)
{
	int N = A.rows();

	// Find determinant of A[][]
	double det = determinant(A, N);
	if (det == 0) 
	{
		cout << "Singular matrix, can't find its inverse";
		return false;
	}

	// Find adjoint
	Matrix adj(N,N);
	adjoint(A, adj);

	// Find Inverse using formula "inverse(A) =
	// adj(A)/det(A)"
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			inverse[i][j] = adj[i][j] / float(det);

	return true;
}

Vector operator*(const Matrix& A, const Vector& v)
{
	int m = A.rows();
	int n = A.cols();

    Vector p(m, 0);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            p[i] += A[i][j] * v[j];
        }
    }
    return p;
}

Matrix operator*(const Matrix& A1, const Matrix& A2)
{

	if (A1.cols() == A2.rows())
	{
		int m = A1.rows();
		int n = A1.cols();
		int p = A2.cols();
		Matrix A(m, p);

		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < p; j++)
			{
				for (int k = 0; k < n; k++)
				{
					A[i][j] += A1[i][k] * A2[k][j];
				}
			}
		}

		return A;
	}
	else
	{
		cout << "The multiplication is not accepted as the number of coloumns of the first matrix does not equal the number of rows of the 2nd matrix";
	}
}

Matrix operator*(const Matrix& A_, const double s)
{
	int m = A_.rows();
	int n = A_.cols();
	Matrix A(m, n);
	for (int i =0 ; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j]=A_[i][j]*s;
		}
	}
	return A;
}

ostream& operator<<(ostream& os, const Matrix& A)
{
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < A.cols() - 1; j++)
        {
            os << A[i][j] << "\t";
        }
        os << A[i][A.cols() - 1] << endl << endl;
    }
    return os;
}

void Matrix::lu_decomposition(Matrix& L, Matrix& U)
{
	Matrix& A = *this;
	int n = A.rows();


	int i = 0, j = 0, k = 0;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (j < i)
			{
				L[j][i] = 0;
			}
			else
			{
				L[j][i] = A[j][i];
				for (k = 0; k < i; k++)
				{
					L[j][i] = L[j][i] - L[j][k] * U[k][i];
				}
			}
		}
		for (j = 0; j < n; j++)
		{
			if (j < i)
			{
				U[i][j] = 0;
			}
			else if (j == i)
			{
				U[i][j] = 1;
			}
			else
			{
				U[i][j] = A[i][j] / L[i][i];
				for (k = 0; k < i; k++)
				{
					U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
				}
			}
		}
	}
}
