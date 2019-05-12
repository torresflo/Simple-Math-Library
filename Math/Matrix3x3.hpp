#ifndef MATH_MATRIX3X3_HPP
#define MATH_MATRIX3X3_HPP

#include <iostream>

#include <Math/Vector.hpp>

namespace Math
{

template <typename T>
class Matrix3x3
{
public:
    Matrix3x3();  // init with identity
    Matrix3x3(const T src[9]);
    Matrix3x3(T m0, T m1, T m2, // 1st column
              T m3, T m4, T m5, // 2nd column
              T m6, T m7, T m8);// 3rd column

    void set(const T src[9]);
    void set(T m0, T m1, T m2,   // 1st column
             T m3, T m4, T m5,   // 2nd column
             T m6, T m7, T m8);  // 3rd column
    void setRow(int index, const T row[3]);
    void setRow(int index, const Vector<T, 3>& v);
    void setColumn(int index, const T col[3]);
    void setColumn(int index, const Vector<T, 3>& v);

    const T* get() const;
    T getDeterminant();

    Matrix3x3& identity();
    Matrix3x3& transpose(); // transpose itself and return reference
    Matrix3x3& invert();

    // operators
    Matrix3x3 operator+(const Matrix3x3& rhs) const;    // add rhs
    Matrix3x3 operator-(const Matrix3x3& rhs) const;    // subtract rhs
    Matrix3x3& operator+=(const Matrix3x3& rhs);        // add rhs and update this object
    Matrix3x3& operator-=(const Matrix3x3& rhs);        // subtract rhs and update this object
    Vector<T, 3> operator*(const Vector<T, 3>& rhs) const;        // multiplication: v' = M * v
    Matrix3x3 operator*(const Matrix3x3& rhs) const;    // multiplication: M3 = M1 * M2
    Matrix3x3& operator*=(const Matrix3x3& rhs);        // multiplication: M1' = M1 * M2
    bool operator==(const Matrix3x3& rhs) const;        // exact compare, no epsilon
    bool operator!=(const Matrix3x3& rhs) const;        // exact compare, no epsilon
    T operator[](int index) const;                      // subscript operator v[0], v[1]
    T& operator[](int index);                           // subscript operator v[0], v[1]

    template <typename U>
    friend Matrix3x3<U> operator-(const Matrix3x3<U>& m);                           // unary operator (-)
    template <typename U>
    friend Matrix3x3<U> operator*(T scalar, const Matrix3x3<U>& m);                 // pre-multiplication
    template <typename U>
    friend Vector<U, 3> operator*(const Vector<U, 3>& vec, const Matrix3x3<U>& m);  // pre-multiplication
    template <typename U>
    friend ::std::ostream& operator<<(::std::ostream& os, const Matrix3x3<U>& m);

private:
    T m[9];
};

} //namespace Math

#include <Math/Matrix3x3.inl>

#endif //MATH_MATRIX3X3_HPP
