
#ifndef Matrix_h
#define Matrix_h


#include <iostream>
#include <vector>


using namespace std;



class Vector;



class Matrix
{
private:

    vector<vector<double>> data_;

    double determinant(Matrix& A, int n);
    void adjoint(Matrix& A, Matrix& adj);
    bool inverse(Matrix& A, Matrix& inverse);
    void get_cofactor(Matrix& A, Matrix& temp, int p, int q, int n);

public:

    Matrix(int m, int n);
    Matrix(std::initializer_list<std::initializer_list<double>> list);

    int rows(void) const { return static_cast<int>(data_.size()); };
    int cols(void) const { return static_cast<int>(data_[0].size()); };


    vector<double>& operator[] (int i) { return data_[i]; };
    const vector<double>& operator[] (int i) const { return data_[i]; };


    Matrix& operator*=(const double s);

    double determinant();
    Matrix inverse();

    void lu_decomposition(Matrix& L, Matrix& U);
};


Vector operator*(const Matrix& A, const Vector& v);
Matrix operator*(const Matrix& A1, const Matrix& A2);
Matrix operator*(const Matrix& A, const double s);

ostream& operator<<(ostream& os, const Matrix& A);


#endif // !Matrix_h
