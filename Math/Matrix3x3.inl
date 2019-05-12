#include <cmath>
#include <algorithm>
#include <iomanip>

#include <Math/MathUtils.hpp>
#include <Math/Matrix2x2.hpp>

namespace Math
{

template<typename T>
Matrix3x3<T>::Matrix3x3()
{
    // initially identity matrix
    identity();
}

template<typename T>
Matrix3x3<T>::Matrix3x3(const T src[9])
{
    set(src);
}

template<typename T>
Matrix3x3<T>::Matrix3x3(T m0, T m1, T m2,
                        T m3, T m4, T m5,
                        T m6, T m7, T m8)
{
    set(m0, m1, m2, m3, m4, m5, m6, m7, m8);
}

template<typename T>
void Matrix3x3<T>::set(const T src[9])
{
    m[0] = src[0];
    m[1] = src[1];
    m[2] = src[2];
    m[3] = src[3];
    m[4] = src[4];
    m[5] = src[5];
    m[6] = src[6];
    m[7] = src[7];
    m[8] = src[8];
}

template<typename T>
void Matrix3x3<T>::set(T m0, T m1, T m2,
                       T m3, T m4, T m5,
                       T m6, T m7, T m8)
{
    m[0] = m0;
    m[1] = m1;
    m[2] = m2;
    m[3] = m3;
    m[4] = m4;
    m[5] = m5;
    m[6] = m6;
    m[7] = m7;
    m[8] = m8;
}

template<typename T>
void Matrix3x3<T>::setRow(int index, const T row[3])
{
    m[index] = row[0];
    m[index + 3] = row[1];
    m[index + 6] = row[2];
}

template<typename T>
void Matrix3x3<T>::setRow(int index, const Vector<T, 3>& v)
{
    m[index] = v.x;
    m[index + 3] = v.y;
    m[index + 6] = v.z;
}

template<typename T>
void Matrix3x3<T>::setColumn(int index, const T col[3])
{
    m[index * 3] = col[0];
    m[index * 3 + 1] = col[1];
    m[index * 3 + 2] = col[2];
}

template<typename T>
void Matrix3x3<T>::setColumn(int index, const Vector<T, 3>& v)
{
    m[index * 3] = v.x;
    m[index * 3 + 1] = v.y;
    m[index * 3 + 2] = v.z;
}

template<typename T>
const T* Matrix3x3<T>::get() const
{
    return m;
}

template<typename T>
T Matrix3x3<T>::getDeterminant()
{
    return m[0] * (m[4] * m[8] - m[5] * m[7]) -
           m[1] * (m[3] * m[8] - m[5] * m[6]) +
           m[2] * (m[3] * m[7] - m[4] * m[6]);
}

template<typename T>
Matrix3x3 <T>& Matrix3x3<T>::identity()
{
    m[0] = m[4] = m[8] = 1.0f;
    m[1] = m[2] = m[3] = m[5] = m[6] = m[7] = 0.0f;
    return *this;
}

template<typename T>
Matrix3x3 <T>& Matrix3x3<T>::transpose()
{
    ::std::swap(m[1], m[3]);
    ::std::swap(m[2], m[6]);
    ::std::swap(m[5], m[7]);

    return *this;
}

template<typename T>
Matrix3x3 <T>& Matrix3x3<T>::invert()
{
    T determinant, invDeterminant;
    T tmp[9];

    tmp[0] = m[4] * m[8] - m[5] * m[7];
    tmp[1] = m[2] * m[7] - m[1] * m[8];
    tmp[2] = m[1] * m[5] - m[2] * m[4];
    tmp[3] = m[5] * m[6] - m[3] * m[8];
    tmp[4] = m[0] * m[8] - m[2] * m[6];
    tmp[5] = m[2] * m[3] - m[0] * m[5];
    tmp[6] = m[3] * m[7] - m[4] * m[6];
    tmp[7] = m[1] * m[6] - m[0] * m[7];
    tmp[8] = m[0] * m[4] - m[1] * m[3];

    // check determinant if it is 0
    determinant = m[0] * tmp[0] + m[1] * tmp[3] + m[2] * tmp[6];
    if (MathUtils::isNullWithEpsilon(fabs(determinant)))
    {
        return identity(); // cannot inverse, make it idenety matrix
    }

    // divide by the determinant
    invDeterminant = 1.0f / determinant;
    m[0] = invDeterminant * tmp[0];
    m[1] = invDeterminant * tmp[1];
    m[2] = invDeterminant * tmp[2];
    m[3] = invDeterminant * tmp[3];
    m[4] = invDeterminant * tmp[4];
    m[5] = invDeterminant * tmp[5];
    m[6] = invDeterminant * tmp[6];
    m[7] = invDeterminant * tmp[7];
    m[8] = invDeterminant * tmp[8];

    return *this;
}

template<typename T>
Matrix3x3 <T> Matrix3x3<T>::operator+(const Matrix3x3 <T>& rhs) const
{
    return Matrix3x3<T>(m[0] + rhs[0], m[1] + rhs[1], m[2] + rhs[2],
                        m[3] + rhs[3], m[4] + rhs[4], m[5] + rhs[5],
                        m[6] + rhs[6], m[7] + rhs[7], m[8] + rhs[8]);
}

template<typename T>
Matrix3x3 <T> Matrix3x3<T>::operator-(const Matrix3x3 <T>& rhs) const
{
    return Matrix3x3<T>(m[0] - rhs[0], m[1] - rhs[1], m[2] - rhs[2],
                        m[3] - rhs[3], m[4] - rhs[4], m[5] - rhs[5],
                        m[6] - rhs[6], m[7] - rhs[7], m[8] - rhs[8]);
}

template<typename T>
Matrix3x3 <T>& Matrix3x3<T>::operator+=(const Matrix3x3 <T>& rhs)
{
    m[0] += rhs[0];
    m[1] += rhs[1];
    m[2] += rhs[2];
    m[3] += rhs[3];
    m[4] += rhs[4];
    m[5] += rhs[5];
    m[6] += rhs[6];
    m[7] += rhs[7];
    m[8] += rhs[8];
    return *this;
}

template<typename T>
Matrix3x3 <T>& Matrix3x3<T>::operator-=(const Matrix3x3 <T>& rhs)
{
    m[0] -= rhs[0];
    m[1] -= rhs[1];
    m[2] -= rhs[2];
    m[3] -= rhs[3];
    m[4] -= rhs[4];
    m[5] -= rhs[5];
    m[6] -= rhs[6];
    m[7] -= rhs[7];
    m[8] -= rhs[8];
    return *this;
}

template<typename T>
Vector<T, 3> Matrix3x3<T>::operator*(const Vector<T, 3>& rhs) const
{
    return Vector<T, 3>(m[0] * rhs.x + m[3] * rhs.y + m[6] * rhs.z,
                        m[1] * rhs.x + m[4] * rhs.y + m[7] * rhs.z,
                        m[2] * rhs.x + m[5] * rhs.y + m[8] * rhs.z);
}

template<typename T>
Matrix3x3 <T> Matrix3x3<T>::operator*(const Matrix3x3 <T>& rhs) const
{
    return Matrix3x3<T>(m[0] * rhs[0] + m[3] * rhs[1] + m[6] * rhs[2], m[1] * rhs[0] + m[4] * rhs[1] + m[7] * rhs[2],
                        m[2] * rhs[0] + m[5] * rhs[1] + m[8] * rhs[2],
                        m[0] * rhs[3] + m[3] * rhs[4] + m[6] * rhs[5], m[1] * rhs[3] + m[4] * rhs[4] + m[7] * rhs[5],
                        m[2] * rhs[3] + m[5] * rhs[4] + m[8] * rhs[5],
                        m[0] * rhs[6] + m[3] * rhs[7] + m[6] * rhs[8], m[1] * rhs[6] + m[4] * rhs[7] + m[7] * rhs[8],
                        m[2] * rhs[6] + m[5] * rhs[7] + m[8] * rhs[8]);
}

template<typename T>
Matrix3x3 <T>& Matrix3x3<T>::operator*=(const Matrix3x3 <T>& rhs)
{
    *this = *this * rhs;
    return *this;
}

template<typename T>
bool Matrix3x3<T>::operator==(const Matrix3x3 <T>& rhs) const
{
    return (m[0] == rhs[0]) && (m[1] == rhs[1]) && (m[2] == rhs[2]) &&
           (m[3] == rhs[3]) && (m[4] == rhs[4]) && (m[5] == rhs[5]) &&
           (m[6] == rhs[6]) && (m[7] == rhs[7]) && (m[8] == rhs[8]);
}

template<typename T>
bool Matrix3x3<T>::operator!=(const Matrix3x3 <T>& rhs) const
{
    return (m[0] != rhs[0]) || (m[1] != rhs[1]) || (m[2] != rhs[2]) ||
           (m[3] != rhs[3]) || (m[4] != rhs[4]) || (m[5] != rhs[5]) ||
           (m[6] != rhs[6]) || (m[7] != rhs[7]) || (m[8] != rhs[8]);
}

template<typename T>
T Matrix3x3<T>::operator[](int index) const
{
    return m[index];
}

template<typename T>
T& Matrix3x3<T>::operator[](int index)
{
    return m[index];
}

template<typename T>
Matrix3x3 <T> operator-(const Matrix3x3 <T>& rhs)
{
    return Matrix3x3<T>(-rhs[0], -rhs[1], -rhs[2], -rhs[3], -rhs[4], -rhs[5], -rhs[6], -rhs[7], -rhs[8]);
}

template<typename T>
Matrix3x3 <T> operator*(T s, const Matrix3x3 <T>& rhs)
{
    return Matrix3x3<T>(s * rhs[0], s * rhs[1], s * rhs[2], s * rhs[3], s * rhs[4], s * rhs[5], s * rhs[6], s * rhs[7],
                        s * rhs[8]);
}

template<typename T>
Vector<T, 3> operator*(const Vector<T, 3>& v, const Matrix3x3 <T>& m)
{
    return Vector<T, 3>(v.x * m[0] + v.y * m[1] + v.z * m[2], v.x * m[3] + v.y * m[4] + v.z * m[5],
                        v.x * m[6] + v.y * m[7] + v.z * m[8]);
}

template<typename T>
::std::ostream& operator<<(::std::ostream& os, const Matrix3x3 <T>& m)
{
    os << ::std::fixed << ::std::setprecision(5);
    os << "[" << ::std::setw(10) << m[0] << " " << ::std::setw(10) << m[3] << " " << ::std::setw(10) << m[6] << "]\n"
       << "[" << ::std::setw(10) << m[1] << " " << ::std::setw(10) << m[4] << " " << ::std::setw(10) << m[7] << "]\n"
       << "[" << ::std::setw(10) << m[2] << " " << ::std::setw(10) << m[5] << " " << ::std::setw(10) << m[8] << "]\n";
    os << ::std::resetiosflags(::std::ios_base::fixed | ::std::ios_base::floatfield);
    return os;
}

} // namespace Math
