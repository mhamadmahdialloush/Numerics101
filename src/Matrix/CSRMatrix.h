#ifndef CSRMatrix_h
#define CSRMatrix_h


#include <iostream>
#include <vector>

#include "Vector.h"


using namespace std;

class CSRMatrix
{
private:
	vector<double> data_;
	vector<int> offsets_;
	vector<int> col_indices_;

	int rows_;
	int cols_;


public:

	CSRMatrix(void);
	CSRMatrix(std::initializer_list<double> data, std::initializer_list<int> offsets, std::initializer_list<int> col_indices);
	CSRMatrix(const vector<double>& data, const vector<int>& offsets, const vector<int>& col_indices);

	vector<double>& data() { return data_; };
	const vector<double>& data() const { return data_; }

	vector<int>& offsets() { return offsets_; };
	const vector<int>& offsets() const { return offsets_; }

	vector<int>& col_indices() { return col_indices_; };
	const vector<int>& col_indices() const { return col_indices_; }

	const int rows() const { return rows_; }
	const int cols() const { return cols_; }

	void set_rows(int rows) { rows_ = rows; }
	void set_cols(int cols) { cols_ = cols; }
};

Vector operator*(const CSRMatrix& A, const Vector& v);

ostream& operator<<(ostream& os, const CSRMatrix& A);



#endif // CSRMatrix_h

