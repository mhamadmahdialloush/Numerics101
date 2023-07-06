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

	for (int i = 0; i < n_rows; i++)
	{
		// pre-zeros
		for (int jj = 0; jj < col_indices[offsets[i]]; jj++)
		{
			os << 0 << "\t";
		}

		for (int j = offsets[i]; j < offsets[i + 1]-1; j++)
		{
			int col_index = col_indices[j];

			os << data[j] << "\t";

			// in-between zeros
			for (int jj=col_indices[j]+1; jj<col_indices[j+1]; jj++)
			{
				os << 0 << "\t";
			}
		}

		//// post-zeros // BODGEEEE
		//for (int jj = col_indices[offsets[i + 1]]; jj < col_indices[offsets[i + 1] - 1]; jj++)
		//{
		//	os << 0 << "\t";
		//}

		os << data[offsets[i+1]-1] << endl << endl;
	}
	return os;
}