#include <cstring>

namespace Math
{

template<typename T, unsigned int N>
inline Vector<T, N>::Vector()
{
    m_values.fill(0);
}

template<typename T, unsigned int N>
inline Vector<T, N>::Vector(const Vector <T, N>& vector)
{
    m_values = vector.m_values;
}

template<typename T, unsigned int N>
template<typename U>
inline Vector<T, N>::Vector(const Vector <U, N>& vector)
{
    ::std::copy(vector.m_values.begin(), vector.m_values.end(), m_values.begin());
}

template<typename T, unsigned int N>
template<unsigned int M>
inline Vector<T, N>::Vector(const Vector <T, M>& vector)
{
    for (int i = 0; i < M, i < N; ++i)
        m_values[i] = vector.m_values[i];
    for (int j = N; j < M; ++j)
        m_values[j] = 0;
}

template<typename T, unsigned int N>
template<typename U, unsigned int M>
inline Vector<T, N>::Vector(const Vector <U, M>& vector)
{
    for (int i = 0; i < M, i < N; ++i)
        m_values[i] = static_cast<T>(vector.m_values[i]);
    for (int j = N; j < M; ++j)
        m_values[j] = 0;
};

template<typename T, unsigned int N>
template<typename ... Values>
Vector<T, N>::Vector(Values... values)
        :m_values{static_cast<T>(values)...}
{
    static_assert(sizeof...(Values) == N, "Invalid number of constructor arguments");
}

template<typename T, unsigned int N>
T& Vector<T, N>::operator[](::std::size_t index)
{
    return m_values[index];
}

template<typename T, unsigned int N>
const T& Vector<T, N>::operator[](::std::size_t index) const
{
    return m_values[index];
}

template<typename T, unsigned int N>
inline Vector<T, N> operator-(const Vector<T, N>& left)
{
    Vector<T, N> vector;
    for (int i = 0; i < N; ++i)
        vector.m_values[i] = -left.m_values[i];
    return vector;
}

template<typename T, unsigned int N>
inline Vector<T, N>& operator+=(Vector<T, N>& left, const Vector<T, N>& right)
{
    for (int i = 0; i < N; ++i)
        left.m_values[i] += right.m_values[i];
    return left;
}

template<typename T, unsigned int N>
inline Vector<T, N>& operator-=(Vector<T, N>& left, const Vector<T, N>& right)
{
    for (int i = 0; i < N; ++i)
        left.m_values[i] -= right.m_values[i];
    return left;
}

template<typename T, unsigned int N>
inline Vector<T, N> operator+(const Vector<T, N>& left, const Vector<T, N>& right)
{
    Vector<T, N> vector;
    for (int i = 0; i < N; ++i)
        vector.m_values[i] = left.m_values[i] + right.m_values[i];
    return vector;
}

template<typename T, unsigned int N>
inline Vector<T, N> operator-(const Vector<T, N>& left, const Vector<T, N>& right)
{
    Vector<T, N> vector;
    for (int i = 0; i < N; ++i)
        vector.m_values[i] = left.m_values[i] - right.m_values[i];
    return vector;
}

template<typename T, unsigned int N>
inline Vector<T, N> operator*(const Vector<T, N>& left, T right)
{
    Vector<T, N> vector;
    for (int i = 0; i < N; ++i)
        vector.m_values[i] = left.m_values[i] * right;
    return vector;
}

template<typename T, unsigned int N>
inline Vector<T, N> operator*(T left, const Vector<T, N>& right)
{
    Vector<T, N> vector;
    for (int i = 0; i < N; ++i)
        vector.m_values[i] = right.m_values[i] * left;
    return vector;
}

template<typename T, unsigned int N>
inline Vector<T, N>& operator*=(Vector<T, N>& left, T right)
{
    for (int i = 0; i < N; ++i)
        left.m_values[i] *= right;
    return left;
}

template<typename T, unsigned int N>
inline Vector<T, N> operator/(const Vector<T, N>& left, T right)
{
    Vector<T, N> vector;
    for (int i = 0; i < N; ++i)
        vector.m_values[i] = left.m_values[i] / right;
    return vector;
}

template<typename T, unsigned int N>
inline Vector<T, N>& operator/=(Vector<T, N>& left, T right)
{
    for (int i = 0; i < N; ++i)
        left.m_values[i] /= right;
    return left;
}

template<typename T, unsigned int N>
inline bool operator==(const Vector<T, N>& left, const Vector<T, N>& right)
{
    for (int i = 0; i < N; ++i)
        if (left.m_values[i] != right.m_values[i])
            return false;
    return true;
}

template<typename T, unsigned int N>
inline bool operator!=(const Vector<T, N>& left, const Vector<T, N>& right)
{
    return !(left == right);
}

} // namespace Math
