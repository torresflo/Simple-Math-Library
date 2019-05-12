#ifndef MATH_MATRIX2X2_HPP
#define MATH_MATRIX2X2_HPP

#include <ostream>

#include <Math/Vector.hpp>

namespace Math
{

template<typename T>
class Matrix2x2
{
public:
    Matrix2x2();  // init with identity
    Matrix2x2(const T src[4]);
    Matrix2x2(T m0, T m1, T m2, T m3);

    void set(const T src[4]);
    void set(T m0, T m1, T m2, T m3);
    void setRow(int index, const T row[2]);
    void setRow(int index, const Vector<T, 2>& v);
    void setColumn(int index, const T col[2]);
    void setColumn(int index, const Vector<T, 2>& v);

    void setAxisX(const T axisX[2]);
    void setAxisX(const Vector<T, 2>& v);
    void setAxisY(const T axisY[2]);
    void setAxisY(const Vector<T, 2>& v);

    void setOrientation(T radians);
    bool isAnOrientation() const;

    Vector<T, 2> getAxisX() const;
    Vector<T, 2> getAxisY() const;
    const T* get() const;
    T getDeterminant() const;

    Matrix2x2& identity();
    Matrix2x2& transpose(); // transpose itself and return reference
    Matrix2x2& invert();

    // operators
    Matrix2x2 operator+(const Matrix2x2& rhs) const;    // add rhs
    Matrix2x2 operator-(const Matrix2x2& rhs) const;    // subtract rhs
    Matrix2x2& operator+=(const Matrix2x2& rhs);        // add rhs and update this object
    Matrix2x2& operator-=(const Matrix2x2& rhs);        // subtract rhs and update this object
    Vector<T, 2> operator*(const Vector<T, 2>& rhs) const;   // multiplication: v' = M * v
    Matrix2x2 operator*(const Matrix2x2& rhs) const;    // multiplication: M3 = M1 * M2
    Matrix2x2& operator*=(const Matrix2x2& rhs);        // multiplication: M1' = M1 * M2
    bool operator==(const Matrix2x2& rhs) const;        // exact compare, no epsilon
    bool operator!=(const Matrix2x2& rhs) const;        // exact compare, no epsilon
    T operator[](int index) const;                      // subscript operator v[0], v[1]
    T& operator[](int index);                           // subscript operator v[0], v[1]

    template<typename U>
    friend Matrix2x2<U> operator-(const Matrix2x2<U>& m);                     // unary operator (-)
    template<typename U>
    friend Matrix2x2<U> operator*(T scalar, const Matrix2x2<U>& m);           // pre-multiplication
    template<typename U>
    friend Vector<U, 2> operator*(const Vector<U, 2>& vec, const Matrix2x2<U>& m);   // pre-multiplication
    template<typename U>
    friend ::std::ostream& operator<<(::std::ostream& os, const Matrix2x2<U>& m);

private:
    T m[4];
};

typedef Matrix2x2<float> Matrix2x2f;

} // namespace Math

#include <Math/Matrix2x2.inl>

#endif //MATH_MATRIX2X2_HPP
