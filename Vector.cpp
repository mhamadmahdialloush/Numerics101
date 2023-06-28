#include "Vector.h"


Vector& Vector::operator=(const Vector& v)
{
    int n = v.size();
    for (int i = 0; i < n; i++)
    {
        (*this)[i] = v[i];
    }
    return (*this);
}


Vector& Vector::operator=(const double val)
{
    int n = (*this).size();
    for (int i = 0; i < n; i++)
    {
        (*this)[i] = val;
    }
    return (*this);
}


Vector operator+(const Vector& v1, const Vector& v2)
{
    int n = v1.size();
    Vector v(n, 0);
    for (int i = 0; i < n; i++)
    {
        v[i] = v1[i] + v2[i];
    }
    return v;
}


Vector operator-(const Vector& v1, const Vector& v2)
{
    int n = v1.size();
    Vector v(n, 0);
    for (int i = 0; i < n; i++)
    {
        v[i] = v1[i] - v2[i];
    }
    return v;
}


ostream& operator<<(ostream& is, const Vector& v)
{
    for (int i = 0; i < v.size(); i++)
    {
        is << v[i] << endl << endl;
    }
    return is;
}

