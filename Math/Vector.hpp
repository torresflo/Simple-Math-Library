#ifndef MATH_VECTOR_HPP
#define MATH_VECTOR_HPP

#include <type_traits>
#include <array>
#include <cmath>

namespace Math
{

template <typename T, unsigned int N>
class Vector
{
public:
    Vector(); //Null vector, all values to 0
    Vector(const Vector<T, N>& vector);

    template <typename U>
    Vector(const Vector<U, N>& vector);

    template <unsigned int M>
    Vector(const Vector<T, M>& vector);

    template <typename U, unsigned int M>
    Vector(const Vector<U, M>& vector);

    template<typename ... Values>
    Vector(Values... values);

    T& operator[](::std::size_t index);
    const T& operator[](::std::size_t index) const;

    template <unsigned int M = N>
    typename ::std::enable_if<M >= 1, T>::type x() const { return m_values[0]; };

    template <unsigned int M = N>
    typename ::std::enable_if<M >= 1, void>::type setX(T x) { m_values[0] = x; }

    template <unsigned int M = N>
    typename ::std::enable_if<M >= 2, T>::type y() const { return m_values[1]; }

    template <unsigned int M = N>
    typename ::std::enable_if<M >= 2, void>::type setY(T y) { m_values[1] = y; }

    template <unsigned int M = N>
    typename ::std::enable_if<M >= 3, T>::type z() const { return m_values[2]; }

    template <unsigned int M = N>
    typename ::std::enable_if<M >= 3, void>::type setZ(T z) { m_values[2] = z; }

    template <unsigned int M = N>
    typename ::std::enable_if<M >= 4, T>::type w() const { return m_values[3]; }

    template <unsigned int M = N>
    typename ::std::enable_if<M >= 4, void>::type setW(T w) { m_values[3] = w; }

    template <typename U = T, unsigned int M = N>
    typename ::std::enable_if<M == 2 && std::is_floating_point<U>::value, U>::type getLengthSqr() const
    {
        return x()*x() + y()*y();
    }

    template <typename U = T, unsigned int M = N>
    typename ::std::enable_if<M == 2 && std::is_floating_point<U>::value, U>::type getLength() const
    {
        return sqrt(getLengthSqr());
    }

    template <typename U = T, unsigned int M = N>
    typename ::std::enable_if<M == 2 && std::is_floating_point<U>::value, void>::type rotate(U radians)
    {
        U cos = ::std::cos(radians);
        U sin = ::std::sin(radians);

        U x = m_values[0] * cos - m_values[1] * sin;
        U y = m_values[0] * sin + m_values[1] * cos;

        m_values[0] = x;
        m_values[1] = y;
    }

    template <typename U = T, unsigned int M = N>
    static typename ::std::enable_if<M >= 2, Vector<U, M>>::type getI()
    {
        Vector<U, M> i;
        i.setX(1);
        return i;
    }

    template <typename U = T, unsigned int M = N>
    static typename ::std::enable_if<M >= 2, Vector<U, M>>::type getJ()
    {
        Vector<U, M> j;
        j.setY(1);
        return j;
    }

    template <typename U = T, unsigned int M = N>
    static typename ::std::enable_if<M >= 3, Vector<U, M>>::type getK()
    {
        Vector<U, M> k;
        k.setZ(1);
        return k;
    }

    template <typename U = T, unsigned int M = N>
    static typename ::std::enable_if<M == 2, Vector<U, M>>::type getPerpendicular(const Vector<U, M>& vector)
    {
        return Vector<U, M>(-vector.y(), vector.x());
    }

    template <typename U = T, unsigned int M = N>
    static typename ::std::enable_if<M == 2, U>::type dot(const Vector<U, M>& a, const Vector<U, M>& b)
    {
        return (a.x() * b.x()) + (a.y() * b.y());
    }

    template <typename U = T, unsigned int M = N>
    static typename ::std::enable_if<M == 2, U>::type cross(const Vector<U, M>& a, const Vector<U, M>& b)
    {
        return a.x() * b.y() - a.y() * b.x();
    }

    template <typename U = T, unsigned int M = N>
    static typename ::std::enable_if<M == 2, Vector<U, M>>::type cross(const Vector<U, M>& a, U scalar)
    {
        return Vector<U, M>(scalar * a.y(), -scalar * a.x());
    }

    template <typename U = T, unsigned int M = N>
    static typename ::std::enable_if<M == 2, Vector<U, M>>::type cross(U scalar, const Vector<U, M>& a)
    {
        return Vector<U, M>(-scalar * a.y(), scalar * a.x());
    }

    template <typename U = T, unsigned int M = N>
    static typename ::std::enable_if<M == 2, bool>::type overlap(const Vector<U, M>& a, const Vector<U, M>& b)
    {
        if( (a.x() <= b.x()) && (b.x() <= a.y()) )
            return true;
        if( (b.x() <= a.x()) && (a.x() <= b.y()) )
            return true;
        return false;
    }

    template <typename U = T, unsigned int M = N>
    static typename ::std::enable_if<M == 2, double>::type distance(const Vector<U, M>& p1, const Vector<U, M>& p2)
    {
        return sqrt(((p1.x() - p2.x()) * (p1.x() - p2.x())) + ((p1.y() - p2.y()) * (p1.y() - p2.y())));
    }

    template <typename U = T, unsigned int M = N>
    static typename ::std::enable_if<M == 2, Vector<U, M>>::type getNormalized(const Vector<U, M>& vector)
    {
        U lenght = vector.getLength();
        return Vector<U, M>(vector.x() / lenght, vector.y() / lenght);
    }

public:
    std::array<T, N> m_values;
};

template <typename T, unsigned int N>
Vector<T, N> operator -(const Vector<T, N>& left);

template <typename T, unsigned int N>
Vector<T, N>& operator +=(Vector<T, N>& left, const Vector<T, N>& right);

template <typename T, unsigned int N>
Vector<T, N>& operator -=(Vector<T, N>& left, const Vector<T, N>& right);

template <typename T, unsigned int N>
Vector<T, N> operator +(const Vector<T, N>& left, const Vector<T, N>& right);

template <typename T, unsigned int N>
Vector<T, N> operator -(const Vector<T, N>& left, const Vector<T, N>& right);

template <typename T, unsigned int N>
Vector<T, N> operator *(const Vector<T, N>& left, T right);

template <typename T, unsigned int N>
Vector<T, N> operator *(T left, const Vector<T, N>& right);

template <typename T, unsigned int N>
Vector<T, N>& operator *=(Vector<T, N>& left, T right);

template <typename T, unsigned int N>
Vector<T, N> operator /(const Vector<T, N>& left, T right);

template <typename T, unsigned int N>
Vector<T, N>& operator /=(Vector<T, N>& left, T right);

template <typename T, unsigned int N>
bool operator ==(const Vector<T, N>& left, const Vector<T, N>& right);

template <typename T, unsigned int N>
bool operator !=(const Vector<T, N>& left, const Vector<T, N>& right);

typedef Vector<int, 2> Vector2i;
typedef Vector<float, 2> Vector2f;
typedef Vector<int, 3> Vector3i;
typedef Vector<float, 3> Vector3f;
typedef Vector<int, 4> Vector4i;
typedef Vector<float, 4> Vector4f;

} //namespace Math

#include <Math/Vector.inl>

#endif //MATH_VECTOR_HPP