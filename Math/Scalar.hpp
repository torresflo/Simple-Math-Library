#ifndef MATH_SCALAR_HPP
#define MATH_SCALAR_HPP

#include <array>

namespace Math
{

template <typename T, unsigned int N>
class Scalar
{
public:
	Scalar(T initValue = T());
	template<typename ... Values>
	Scalar(Values... dimensions);
	Scalar(const Scalar& copy);

	T& operator[] (unsigned int index);
	const T& operator[] (unsigned int index) const;

private:
	std::array<T, N> m_Values;
};

using Scalar8UC3 = Scalar<unsigned char, 3>;
using Scalar32FC3 = Scalar<float, 3>;

}

#include <Common\Math\Scalar.inl>

#endif //MATH_SCALAR_HPP