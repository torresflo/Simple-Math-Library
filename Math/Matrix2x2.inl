#include <cmath>
#include <algorithm>
#include <iomanip>

#include <Math/MathUtils.hpp>

namespace Math
{

template <typename T>
Matrix2x2<T>::Matrix2x2()
{
    // initially identity matrix
    identity();
}

template <typename T>
Matrix2x2<T>::Matrix2x2(const T src[4])
{
    set(src);
}

template <typename T>
Matrix2x2<T>::Matrix2x2(T m0, T m1, T m2, T m3)
{
    set(m0, m1, m2, m3);
}

template <typename T>
void Matrix2x2<T>::set(const T src[4])
{
    m[0] = src[0];  m[1] = src[1];  m[2] = src[2];  m[3] = src[3];
}

template <typename T>
void Matrix2x2<T>::set(T m0, T m1, T m2, T m3)
{
    m[0]= m0;  m[1] = m1;  m[2] = m2;  m[3]= m3;
}

template <typename T>
void Matrix2x2<T>::setRow(int index, const T row[2])
{
    m[index] = row[0];  m[index + 2] = row[1];
}

template <typename T>
void Matrix2x2<T>::setRow(int index, const Vector<T, 2>& v)
{
    m[index] = v.x;  m[index + 2] = v.y;
}

template <typename T>
void Matrix2x2<T>::setColumn(int index, const T col[2])
{
    m[index*2] = col[0];  m[index*2 + 1] = col[1];
}

template <typename T>
void Matrix2x2<T>::setColumn(int index, const Vector<T, 2>& v)
{
    m[index*2] = v.x;  m[index*2 + 1] = v.y;
}

template <typename T>
void Matrix2x2<T>::setAxisX(const T axisX[2])
{
    setRow(0, axisX);
}

template <typename T>
void Matrix2x2<T>::setAxisX(const Vector<T, 2>& v)
{
    setRow(0, v);
}

template <typename T>
void Matrix2x2<T>::setAxisY(const T axisY[2])
{
    setRow(1, axisY);
}

template <typename T>
void Matrix2x2<T>::setAxisY(const Vector<T, 2>& v)
{
    setRow(1, v);
}

template <typename T>
void Matrix2x2<T>::setOrientation(T radians)
{
    T cos = ::std::cos(radians);
    T sin = ::std::sin(radians);

    m[0] = cos;
    m[1] = sin;
    m[2] = -sin;
    m[3] = cos;
}

template <typename T>
bool Matrix2x2<T>::isAnOrientation() const
{
    return ::std::acos(m[0]) == ::std::asin(m[1]) == ::std::asin(-m[2]) == ::std::acos(m[3]); //TODO Epsilon!
}

template <typename T>
Vector<T, 2> Matrix2x2<T>::getAxisX() const
{
    return Vector<T, 2>(m[0], m[1]);
}

template <typename T>
Vector<T, 2> Matrix2x2<T>::getAxisY() const
{
    return Vector<T, 2>(m[2], m[3]);
}

template <typename T>
const T* Matrix2x2<T>::get() const
{
    return m;
}

template <typename T>
T Matrix2x2<T>::getDeterminant() const
{
    return m[0] * m[3] - m[1] * m[2];
}

template <typename T>
Matrix2x2<T>& Matrix2x2<T>::identity()
{
    m[0] = m[3] = 1.0f;
    m[1] = m[2] = 0.0f;
    return *this;
}

template <typename T>
Matrix2x2<T>& Matrix2x2<T>::transpose()
{
    ::std::swap(m[1],  m[2]);
    return *this;
}

template <typename T>
Matrix2x2<T>& Matrix2x2<T>::invert()
{
    T determinant = getDeterminant();
    if(MathUtils::isNullWithEpsilon(fabs(determinant)))
    {
        return identity();
    }

    T tmp = m[0];   // copy the first element
    T invDeterminant = 1.0f / determinant;
    m[0] =  invDeterminant * m[3];
    m[1] = -invDeterminant * m[1];
    m[2] = -invDeterminant * m[2];
    m[3] =  invDeterminant * tmp;

    return *this;
}

template <typename T>
Matrix2x2<T> Matrix2x2<T>::operator+(const Matrix2x2& rhs) const
{
    return Matrix2x2(m[0]+rhs[0], m[1]+rhs[1], m[2]+rhs[2], m[3]+rhs[3]);
}

template <typename T>
Matrix2x2<T> Matrix2x2<T>::operator-(const Matrix2x2& rhs) const
{
    return Matrix2x2(m[0]-rhs[0], m[1]-rhs[1], m[2]-rhs[2], m[3]-rhs[3]);
}

template <typename T>
Matrix2x2<T>& Matrix2x2<T>::operator+=(const Matrix2x2& rhs)
{
    m[0] += rhs[0];  m[1] += rhs[1];  m[2] += rhs[2];  m[3] += rhs[3];
    return *this;
}

template <typename T>
Matrix2x2<T>& Matrix2x2<T>::operator-=(const Matrix2x2& rhs)
{
    m[0] -= rhs[0];  m[1] -= rhs[1];  m[2] -= rhs[2];  m[3] -= rhs[3];
    return *this;
}

template <typename T>
Vector<T, 2> Matrix2x2<T>::operator*(const Vector<T, 2>& rhs) const
{
    return Vector<T, 2>(m[0]*rhs.x + m[2]*rhs.y,  m[1]*rhs.x + m[3]*rhs.y);
}

template <typename T>
Matrix2x2<T> Matrix2x2<T>::operator*(const Matrix2x2& rhs) const
{
    return Matrix2x2<T>(m[0]*rhs[0] + m[2]*rhs[1],  m[1]*rhs[0] + m[3]*rhs[1],
                   m[0]*rhs[2] + m[2]*rhs[3],  m[1]*rhs[2] + m[3]*rhs[3]);
}

template <typename T>
Matrix2x2<T>& Matrix2x2<T>::operator*=(const Matrix2x2& rhs)
{
    *this = *this * rhs;
    return *this;
}

template <typename T>
bool Matrix2x2<T>::operator==(const Matrix2x2<T>& rhs) const
{
    return (m[0] == rhs[0]) && (m[1] == rhs[1]) && (m[2] == rhs[2]) && (m[3] == rhs[3]);
}

template <typename T>
bool Matrix2x2<T>::operator!=(const Matrix2x2<T>& rhs) const
{
    return (m[0] != rhs[0]) || (m[1] != rhs[1]) || (m[2] != rhs[2]) || (m[3] != rhs[3]);
}

template <typename T>
T Matrix2x2<T>::operator[](int index) const
{
    return m[index];
}

template <typename T>
T& Matrix2x2<T>::operator[](int index)
{
    return m[index];
}

template <typename T>
Matrix2x2<T> operator-(const Matrix2x2<T>& rhs)
{
    return Matrix2x2<T>(-rhs[0], -rhs[1], -rhs[2], -rhs[3]);
}

template <typename T>
Matrix2x2<T> operator*(T s, const Matrix2x2<T>& rhs)
{
    return Matrix2x2<T>(s*rhs[0], s*rhs[1], s*rhs[2], s*rhs[3]);
}

template <typename T>
Vector<T, 2> operator*(const Vector<T, 2>& v, const Matrix2x2<T>& rhs)
{
    return Vector<T, 2>(v.x*rhs[0] + v.y*rhs[1],  v.x*rhs[2] + v.y*rhs[3]);
}

template <typename T>
::std::ostream& operator<<(::std::ostream& os, const Matrix2x2<T>& m)
{
    os << ::std::fixed << ::std::setprecision(5);
    os << "[" << ::std::setw(10) << m[0] << " " << ::std::setw(10) << m[2] << "]\n"
       << "[" << ::std::setw(10) << m[1] << " " << ::std::setw(10) << m[3] << "]\n";
    os << ::std::resetiosflags(::std::ios_base::fixed | ::std::ios_base::floatfield);
    return os;
}

} // namespace Math
