#include <cassert>
#include "Scalar.hpp"

namespace Math
{

template<typename T, unsigned int N>
inline Scalar<T, N>::Scalar(T initValue)
{
	m_Values.fill(initValue);
}

template<typename T, unsigned int N>
template<typename ...Values>
inline Scalar<T, N>::Scalar(Values ...values)
	:m_Values{ static_cast<T>(values)... }
{
}

template<typename T, unsigned int N>
inline Scalar<T, N>::Scalar(const Scalar & copy)
{
	m_Values = copy.m_Values;
}

template<typename T, unsigned int N>
inline T& Scalar<T, N>::operator[](unsigned int index)
{
	assert(index < m_Values.size());
	return m_Values[index];
}

template<typename T, unsigned int N>
inline const T& Scalar<T, N>::operator[](unsigned int index) const
{
	assert(index < m_Values.size());
	return m_Values[index];
}

}