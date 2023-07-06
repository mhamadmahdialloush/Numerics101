

#ifndef Vector_h
#define Vector_h


#include <iostream>
#include <vector>


using namespace std;


class Vector : public vector<double>
{
private:

public:

    Vector() : vector<double>()
    {
    }

    Vector(const int sz) : vector<double>(sz, 0)
    {
    }

    Vector(const int sz, const double val) : vector<double>(sz, val)
    {
    }

    Vector(const Vector& v) : vector<double>(v.size(), 0)
    {
        for (int i = 0; i < v.size(); i++)
        {
            (*this)[i] = v[i];
        }
    }

    Vector(std::initializer_list<double> list) : vector<double>(list)
    {

    }

    Vector& operator=(const Vector& v);
    Vector& operator=(const double val);

    int size(void) { return static_cast<int>(vector<double>::size()); };
    const int size(void) const { return static_cast<int>(vector<double>::size()); };
};


Vector operator+(const Vector& v1, const Vector& v2);

Vector operator-(const Vector& v1, const Vector& v2);

ostream& operator<<(ostream& os, const Vector& v);


#endif // Vector_h
