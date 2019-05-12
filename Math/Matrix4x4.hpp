#ifndef MATH_MATRIX4X4_HPP
#define MATH_MATRIX4X4_HPP

#include <ostream>

#include <Math/Vector.hpp>

namespace Math
{

template<typename T>
class Matrix4x4
{
public:
    Matrix4x4();  // init with identity
    Matrix4x4(const T src[16]);
    Matrix4x4(T m00, T m01, T m02, T m03, // 1st column
              T m04, T m05, T m06, T m07, // 2nd column
              T m08, T m09, T m10, T m11, // 3rd column
              T m12, T m13, T m14, T m15);// 4th column
    Matrix4x4(const Matrix4x4<T>& copy);

    void set(const T src[16]);
    void set(T m00, T m01, T m02, T m03, // 1st column
             T m04, T m05, T m06, T m07, // 2nd column
             T m08, T m09, T m10, T m11, // 3rd column
             T m12, T m13, T m14, T m15);// 4th column
    void setRow(int index, const T row[4]);
    void setRow(int index, const Vector<T, 4>& v);
    void setRow(int index, const Vector<T, 3>& v);
    void setColumn(int index, const T col[4]);
    void setColumn(int index, const Vector<T, 4>& v);
    void setColumn(int index, const Vector<T, 3>& v);

    const T* get() const;
    const T* getTranspose(); // return transposed matrix
T  getDeterminant();

    Matrix4x4<T>& identity();
    Matrix4x4<T>& transpose();         // transpose itself and return reference
    Matrix4x4<T>& invert();            // check best inverse method before inverse
    Matrix4x4<T>& invertEuclidean();   // inverse of Euclidean transform matrix
    Matrix4x4<T>& invertAffine();      // inverse of affine transform matrix
    Matrix4x4<T>& invertProjective();  // inverse of projective matrix using partitioning
    Matrix4x4<T>& invertGeneral();     // inverse of generic matrix

    // transform matrix
    Matrix4x4<T>& translate(T x, T y, T z);                    // translation by (x,y,z)
    Matrix4x4<T>& translate(const Vector<T, 3>& v);
    Matrix4x4<T>& rotate(T angle, const Vector<T, 3>& axis);   // rotate angle(degree) along the given axix
    Matrix4x4<T>& rotate(T angle, T x, T y, T z);
    Matrix4x4<T>& rotateX(T angle);                            // rotate on X-axis with degree
    Matrix4x4<T>& rotateY(T angle);                            // rotate on Y-axis with degree
    Matrix4x4<T>& rotateZ(T angle);                            // rotate on Z-axis with degree
    Matrix4x4<T>& scale(T scale);                              // uniform scale
    Matrix4x4<T>& scale(T sx, T sy, T sz);                     // scale by (sx, sy, sz) on each axis

    // operators
    Matrix4x4 operator+(const Matrix4x4<T>& rhs) const;     // add rhs
    Matrix4x4 operator-(const Matrix4x4<T>& rhs) const;     // subtract rhs
    Matrix4x4& operator+=(const Matrix4x4<T>& rhs);         // add rhs and update this object
    Matrix4x4& operator-=(const Matrix4x4<T>& rhs);         // subtract rhs and update this object
    Vector<T, 4> operator*(const Vector<T, 4>& rhs) const;  // multiplication: v' = M * v
    Vector<T, 3> operator*(const Vector<T, 3>& rhs) const;  // multiplication: v' = M * v
    Matrix4x4<T> operator*(const Matrix4x4<T>& rhs) const;  // multiplication: M3 = M1 * M2
    Matrix4x4<T>& operator*=(const Matrix4x4<T>& rhs);      // multiplication: M1' = M1 * M2
    bool operator==(const Matrix4x4<T>& rhs) const;         // exact compare, no epsilon
    bool operator!=(const Matrix4x4<T>& rhs) const;         // exact compare, no epsilon
    T operator[](int index) const;                          // subscript operator v[0], v[1]
    T& operator[](int index);                               // subscript operator v[0], v[1]

    template<typename U>
    friend Matrix4x4<U> operator-(const Matrix4x4<U>& m);                             // unary operator (-)
    template<typename U>
    friend Matrix4x4<U> operator*(U scalar, const Matrix4x4<U>& m);                // pre-multiplication
    template<typename U>
    friend Vector<U, 3> operator*(const Vector<U, 3>& vec, const Matrix4x4<U>& m); // pre-multiplication
    template<typename U>
    friend Vector<U, 4> operator*(const Vector<U, 4>& vec, const Matrix4x4<U>& m); // pre-multiplication
    template<typename U>
    friend ::std::ostream& operator<<(::std::ostream& os, const Matrix4x4<U>& m);

private:
    T getCofactor(T m0, T m1, T m2,
                  T m3, T m4, T m5,
                  T m6, T m7, T m8);

    T m[16];
    T tm[16]; // transpose m
};

typedef Matrix4x4<float> Matrix4x4f;

} //namespace Math

#include <Math/Matrix4x4.inl>

#endif //MATH_MATRIX4X4_HPP
