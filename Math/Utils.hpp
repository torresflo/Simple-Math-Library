#ifndef MATH_MATHUTILS_HPP
#define MATH_MATHUTILS_HPP

namespace Math
{

template <typename T, unsigned int N>
class Vector;

class Utils
{
public:
    enum AngleUnit
    {
        Degree,
        Radian
    };

    static const float PI;
    static const float PI2;
    static const float DEG_TO_RAD;
    static const float EPSILON;

    template <typename T>
    static bool isNullWithEpsilon(const T& value);

    template <typename T>
    static bool isEqualWithEpsilon(T a, T b);

    template <typename T, unsigned int N>
    static bool isNullWithEpsilon(const Vector<T, N>& vector);

    template <typename T, unsigned int N>
    static Vector<T, N> getMinVector(const Vector <T, N>& a, const Vector <T, N>& b);

    template <typename T, unsigned int N>
    static Vector<T, N> getMaxVector(const Vector <T, N>& a, const Vector <T, N>& b);

    template <typename T>
    static T clamp(T value, T min, T max);

    static float convertAngle(float value, AngleUnit from, AngleUnit to);
    static float normalizeAngle(float value, AngleUnit unit);
};

} //namespace Math

#include <Math/Utils.inl>

#endif //MATH_MATHUTILS_HPP
