#include <Math\Vector.hpp>

namespace Math
{

const float Utils::PI = 3.14159265359f;
const float Utils::PI2 = 2 * PI;
const float Utils::DEG_TO_RAD = Utils::PI / 180.0f;
const float Utils::EPSILON = 0.00001f;

template<typename T>
inline bool Utils::isNullWithEpsilon(const T& value)
{
    if(value > EPSILON)
        return false;
    return true;
}

template<typename T>
inline bool Utils::isEqualWithEpsilon(T a, T b)
{
    return (abs(a) - abs(b)) <= EPSILON;
}

template<typename T, unsigned int N>
inline bool Utils::isNullWithEpsilon(const Vector<T, N>& vector)
{
	for (unsigned int i = 0; i < N; ++i)
	{
		if(vector[i] > EPSILON)
			return false;
	}
	return true;
}

template<typename T, unsigned int N>
inline Vector<T, N> Utils::getMinVector(const Vector <T, N>& a, const Vector <T, N>& b)
{
    Vector<T, N> vector;
    for(unsigned int i = 0; i < N; ++i)
    {
        vector[i] = min(a[i], b[i]);
    }
    return vector;
}

template<typename T, unsigned int N>
inline Vector<T, N> Utils::getMaxVector(const Vector <T, N>& a, const Vector <T, N>& b)
{
    Vector<T, N> vector;
    for(unsigned int i = 0; i < N; ++i)
    {
        vector[i] = max(a[i], b[i]);
    }
    return vector;
}

template<typename T>
inline T Utils::clamp(T value, T min, T max)
{
    if(value < min)
        return min;
    if(value > max)
        return max;
    return value;
}

inline float Utils::convertAngle(float value, Utils::AngleUnit from, Utils::AngleUnit to)
{
    switch(from)
    {
        case Utils::AngleUnit::Degree:
        {
            switch(to)
            {
                case Utils::AngleUnit::Degree:
                    return value;
                case Utils::AngleUnit::Radian:
                    return value / 180 * Utils::PI;
            }
        }
        case Utils::AngleUnit::Radian:
        {
            switch(to)
            {
                case Utils::AngleUnit::Degree:
                    return value / Utils::PI * 180;
                case Utils::AngleUnit::Radian:
                    return value;
            }
        }
    }
}

inline float Utils::normalizeAngle(float value, Utils::AngleUnit unit)
{
    float min = 0.f;
    float max = 0.f;
    switch (unit)
    {
        case Utils::AngleUnit::Degree:
        {
            max = 360.f;
            break;
        }
        case Utils::AngleUnit::Radian:
        {
            max = Utils::PI2;
        }
    }

    if(value < min || value >= max)
    {
        float i = value / 360.f;
        value -= (i * max);
        if(value < 0.f)
        {
            value += max;
        }
    }
    return value;
}

} //namespace Math