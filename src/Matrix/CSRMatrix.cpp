#include "CSRMatrix.h"

CSRMatrix::CSRMatrix(void)
{

}


CSRMatrix::CSRMatrix(std::initializer_list<double> data, std::initializer_list<int> offsets, std::initializer_list<int> col_indices) : data_(data), offsets_(offsets), col_indices_(col_indices)
{

}


CSRMatrix::CSRMatrix(const vector<double>& data, const vector<int>& offsets, const vector<int>& col_indices) : data_(data), offsets_(offsets), col_indices_(col_indices)
{

}


Vector operator*(const CSRMatrix& A, const Vector& v)
{
	int n_rows = A.rows();

	Vector v_out(n_rows);	

	const vector<double>& data = A.data();
	const vector<int>& offsets = A.offsets();
	const vector<int>& col_indices = A.col_indices();

	for (int i = 0; i < n_rows; i++)
	{
		for (int j = offsets[i]; j < offsets[i+1]; j++)
		{
			int col_index = col_indices[j];
			v_out[i] += data[j] * v[col_index];
		}
	}

	return v_out;
}


ostream& operator<<(ostream& os, const CSRMatrix& A)
{
	const vector<double>& data = A.data();
	const vector<int>& offsets = A.offsets();
	const vector<int>& col_indices = A.col_indices();

	int n_rows = A.rows();

	// Calculate the number of coloumns. This might be wrong in
	// certain rare situations when the last columns of the matrix
	// have no non-zeros
	int n_cols = 0;
	for (int i = 0; i < col_indices.size(); i++)
	{
		if (col_indices[i] > n_cols)
		{
			n_cols = col_indices[i];
		}
	}
	n_cols++;

	int current = 0; // index for data array

	for (int i = 0; i < n_rows; i++) 
	{
		for (int j = 0; j < n_cols; j++) 
		{
			if (current < offsets[i + 1] && j == col_indices[current]) 
			{
				os << data[current] << "\t";
				current++;
			}
			else 
			{
				os << "0\t";
			}
		}
		os << endl << endl;
	}
	
	return os;
}