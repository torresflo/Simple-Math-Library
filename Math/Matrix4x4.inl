#include <cmath>
#include <algorithm>

#include <Math/MathUtils.hpp>
#include <Math/Matrix2x2.hpp>
#include <Math/Matrix3x3.hpp>

namespace Math
{

template<typename T>
Matrix4x4<T>::Matrix4x4()
{
    // initially identity matrix
    identity();
}

template<typename T>
Matrix4x4<T>::Matrix4x4(const T src[16])
{
    set(src);
}

template<typename T>
Matrix4x4<T>::Matrix4x4(T m00, T m01, T m02, T m03,
                        T m04, T m05, T m06, T m07,
                        T m08, T m09, T m10, T m11,
                        T m12, T m13, T m14, T m15)
{
    set(m00, m01, m02, m03,  m04, m05, m06, m07,  m08, m09, m10, m11,  m12, m13, m14, m15);
}

template<typename T>
void Matrix4x4<T>::set(const T src[16])
{
    m[0] = src[0];  m[1] = src[1];  m[2] = src[2];  m[3] = src[3];
    m[4] = src[4];  m[5] = src[5];  m[6] = src[6];  m[7] = src[7];
    m[8] = src[8];  m[9] = src[9];  m[10]= src[10]; m[11]= src[11];
    m[12]= src[12]; m[13]= src[13]; m[14]= src[14]; m[15]= src[15];
}

template<typename T>
void Matrix4x4<T>::set(T m00, T m01, T m02, T m03,
                         T m04, T m05, T m06, T m07,
                         T m08, T m09, T m10, T m11,
                         T m12, T m13, T m14, T m15)
{
    m[0] = m00;  m[1] = m01;  m[2] = m02;  m[3] = m03;
    m[4] = m04;  m[5] = m05;  m[6] = m06;  m[7] = m07;
    m[8] = m08;  m[9] = m09;  m[10]= m10;  m[11]= m11;
    m[12]= m12;  m[13]= m13;  m[14]= m14;  m[15]= m15;
}

template<typename T>
void Matrix4x4<T>::setRow(int index, const T row[4])
{
    m[index] = row[0];  m[index + 4] = row[1];  m[index + 8] = row[2];  m[index + 12] = row[3];
}

template<typename T>
void Matrix4x4<T>::setRow(int index, const Vector<T, 4>& v)
{
    m[index] = v.x;  m[index + 4] = v.y;  m[index + 8] = v.z;  m[index + 12] = v.w;
}

template<typename T>
void Matrix4x4<T>::setRow(int index, const Vector<T, 3>& v)
{
    m[index] = v.x;  m[index + 4] = v.y;  m[index + 8] = v.z;
}

template<typename T>
void Matrix4x4<T>::setColumn(int index, const T col[4])
{
    m[index*4] = col[0];  m[index*4 + 1] = col[1];  m[index*4 + 2] = col[2];  m[index*4 + 3] = col[3];
}

template<typename T>
void Matrix4x4<T>::setColumn(int index, const Vector<T, 4>& v)
{
    m[index*4] = v.x;  m[index*4 + 1] = v.y;  m[index*4 + 2] = v.z;  m[index*4 + 3] = v.w;
}

template<typename T>
void Matrix4x4<T>::setColumn(int index, const Vector<T, 3>& v)
{
    m[index*4] = v.x;  m[index*4 + 1] = v.y;  m[index*4 + 2] = v.z;
}

template<typename T>
const T* Matrix4x4<T>::get() const
{
    return m;
}

template<typename T>
const T* Matrix4x4<T>::getTranspose()
{
    tm[0] = m[0];   tm[1] = m[4];   tm[2] = m[8];   tm[3] = m[12];
    tm[4] = m[1];   tm[5] = m[5];   tm[6] = m[9];   tm[7] = m[13];
    tm[8] = m[2];   tm[9] = m[6];   tm[10]= m[10];  tm[11]= m[14];
    tm[12]= m[3];   tm[13]= m[7];   tm[14]= m[11];  tm[15]= m[15];
    return tm;
}

template<typename T>
T Matrix4x4<T>::getDeterminant()
{
    return m[0] * getCofactor(m[5],m[6],m[7], m[9],m[10],m[11], m[13],m[14],m[15]) -
           m[1] * getCofactor(m[4],m[6],m[7], m[8],m[10],m[11], m[12],m[14],m[15]) +
           m[2] * getCofactor(m[4],m[5],m[7], m[8],m[9], m[11], m[12],m[13],m[15]) -
           m[3] * getCofactor(m[4],m[5],m[6], m[8],m[9], m[10], m[12],m[13],m[14]);
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::identity()
{
    m[0] = m[5] = m[10] = m[15] = 1.0f;
    m[1] = m[2] = m[3] = m[4] = m[6] = m[7] = m[8] = m[9] = m[11] = m[12] = m[13] = m[14] = 0.0f;
    return *this;
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::transpose()
{
    ::std::swap(m[1],  m[4]);
    ::std::swap(m[2],  m[8]);
    ::std::swap(m[3],  m[12]);
    ::std::swap(m[6],  m[9]);
    ::std::swap(m[7],  m[13]);
    ::std::swap(m[11], m[14]);

    return *this;
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::invert()
{
    // If the 4th row is [0,0,0,1] then it is affine matrix and
    // it has no projective transformation.
    if(m[3] == 0 && m[7] == 0 && m[11] == 0 && m[15] == 1)
        this->invertAffine();
    else
    {
        this->invertGeneral();
        /*@@ invertProjective() is not optimized (slower than generic one)
        if(fabs(m[0]*m[5] - m[1]*m[4]) > MathUtils::EPSILON)
            this->invertProjective();   // inverse using matrix partition
        else
            this->invertGeneral();      // generalized inverse
        */
    }

    return *this;
}

///////////////////////////////////////////////////////////////////////////////
// compute the inverse of 4x4 Euclidean transformation matrix
//
// Euclidean transformation is translation, rotation, and reflection.
// With Euclidean transform, only the position and orientation of the object
// will be changed. Euclidean transform does not change the shape of an object
// (no scaling). Length and angle are reserved.
//
// Use inverseAffine() if the matrix has scale and shear transformation.
//
// M = [ R | T ]
//     [ --+-- ]    (R denotes 3x3 rotation/reflection matrix)
//     [ 0 | 1 ]    (T denotes 1x3 translation matrix)
//
// y = M*x  ->  y = R*x + T  ->  x = R^-1*(y - T)  ->  x = R^T*y - R^T*T
// (R is orthogonal,  R^-1 = R^T)
//
//  [ R | T ]-1    [ R^T | -R^T * T ]    (R denotes 3x3 rotation matrix)
//  [ --+-- ]   =  [ ----+--------- ]    (T denotes 1x3 translation)
//  [ 0 | 1 ]      [  0  |     1    ]    (R^T denotes R-transpose)
///////////////////////////////////////////////////////////////////////////////
template<typename T>
Matrix4x4<T>& Matrix4x4<T>::invertEuclidean()
{
    // transpose 3x3 rotation matrix part
    // | R^T | 0 |
    // | ----+-- |
    // |  0  | 1 |
    T tmp;
    tmp = m[1];  m[1] = m[4];  m[4] = tmp;
    tmp = m[2];  m[2] = m[8];  m[8] = tmp;
    tmp = m[6];  m[6] = m[9];  m[9] = tmp;

    // compute translation part -R^T * T
    // | 0 | -R^T x |
    // | --+------- |
    // | 0 |   0    |
    T x = m[12];
    T y = m[13];
    T z = m[14];
    m[12] = -(m[0] * x + m[4] * y + m[8] * z);
    m[13] = -(m[1] * x + m[5] * y + m[9] * z);
    m[14] = -(m[2] * x + m[6] * y + m[10]* z);

    // last row should be unchanged (0,0,0,1)

    return *this;
}

///////////////////////////////////////////////////////////////////////////////
// compute the inverse of a 4x4 affine transformation matrix
//
// Affine transformations are generalizations of Euclidean transformations.
// Affine transformation includes translation, rotation, reflection, scaling,
// and shearing. Length and angle are NOT preserved.
// M = [ R | T ]
//     [ --+-- ]    (R denotes 3x3 rotation/scale/shear matrix)
//     [ 0 | 1 ]    (T denotes 1x3 translation matrix)
//
// y = M*x  ->  y = R*x + T  ->  x = R^-1*(y - T)  ->  x = R^-1*y - R^-1*T
//
//  [ R | T ]-1   [ R^-1 | -R^-1 * T ]
//  [ --+-- ]   = [ -----+---------- ]
//  [ 0 | 1 ]     [  0   +     1     ]
///////////////////////////////////////////////////////////////////////////////
template<typename T>
Matrix4x4<T>& Matrix4x4<T>::invertAffine()
{
    // R^-1
    Matrix3x3<T> r(m[0],m[1],m[2], m[4],m[5],m[6], m[8],m[9],m[10]);
    r.invert();
    m[0] = r[0];  m[1] = r[1];  m[2] = r[2];
    m[4] = r[3];  m[5] = r[4];  m[6] = r[5];
    m[8] = r[6];  m[9] = r[7];  m[10]= r[8];

    // -R^-1 * T
    T x = m[12];
    T y = m[13];
    T z = m[14];
    m[12] = -(r[0] * x + r[3] * y + r[6] * z);
    m[13] = -(r[1] * x + r[4] * y + r[7] * z);
    m[14] = -(r[2] * x + r[5] * y + r[8] * z);

    // last row should be unchanged (0,0,0,1)
    //m[3] = m[7] = m[11] = 0.0f;
    //m[15] = 1.0f;

    return * this;
}

///////////////////////////////////////////////////////////////////////////////
// inverse matrix using matrix partitioning (blockwise inverse)
// It devides a 4x4 matrix into 4 of 2x2 matrices. It works in case of where
// det(A) != 0. If not, use the generic inverse method
// inverse formula.
// M = [ A | B ]    A, B, C, D are 2x2 matrix blocks
//     [ --+-- ]    det(M) = |A| * |D - ((C * A^-1) * B)|
//     [ C | D ]
//
// M^-1 = [ A' | B' ]   A' = A^-1 - (A^-1 * B) * C'
//        [ ---+--- ]   B' = (A^-1 * B) * -D'
//        [ C' | D' ]   C' = -D' * (C * A^-1)
//                      D' = (D - ((C * A^-1) * B))^-1
//
// NOTE: I wrap with () if it it used more than once.
//       The matrix is invertable even if det(A)=0, so must check det(A) before
//       calling this function, and use invertGeneric() instead.
///////////////////////////////////////////////////////////////////////////////
template<typename T>
Matrix4x4<T>& Matrix4x4<T>::invertProjective()
{
    // partition
    Matrix2x2<T> a(m[0], m[1], m[4], m[5]);
    Matrix2x2<T> b(m[8], m[9], m[12], m[13]);
    Matrix2x2<T> c(m[2], m[3], m[6], m[7]);
    Matrix2x2<T> d(m[10], m[11], m[14], m[15]);

    // pre-compute repeated parts
    a.invert();             // A^-1
    Matrix2x2<T> ab = a * b;     // A^-1 * B
    Matrix2x2<T> ca = c * a;     // C * A^-1
    Matrix2x2<T> cab = ca * b;   // C * A^-1 * B
    Matrix2x2<T> dcab = d - cab; // D - C * A^-1 * B

    // check determinant if |D - C * A^-1 * B| = 0
    //NOTE: this function assumes det(A) is already checked. if |A|=0 then,
    //      cannot use this function.
    T determinant = dcab[0] * dcab[3] - dcab[1] * dcab[2];
    if(MathUtils::isNullWithEpsilon(fabs(determinant)))
    {
        return identity();
    }

    // compute D' and -D'
    Matrix2x2<T> d1 = dcab;      //  (D - C * A^-1 * B)
    d1.invert();            //  (D - C * A^-1 * B)^-1
    Matrix2x2<T> d2 = -d1;       // -(D - C * A^-1 * B)^-1

    // compute C'
    Matrix2x2<T> c1 = d2 * ca;   // -D' * (C * A^-1)

    // compute B'
    Matrix2x2<T> b1 = ab * d2;   // (A^-1 * B) * -D'

    // compute A'
    Matrix2x2<T> a1 = a - (ab * c1); // A^-1 - (A^-1 * B) * C'

    // assemble inverse matrix
    m[0] = a1[0];  m[4] = a1[2]; /*|*/ m[8] = b1[0];  m[12]= b1[2];
    m[1] = a1[1];  m[5] = a1[3]; /*|*/ m[9] = b1[1];  m[13]= b1[3];
    /*-----------------------------+-----------------------------*/
    m[2] = c1[0];  m[6] = c1[2]; /*|*/ m[10]= d1[0];  m[14]= d1[2];
    m[3] = c1[1];  m[7] = c1[3]; /*|*/ m[11]= d1[1];  m[15]= d1[3];

    return *this;
}

///////////////////////////////////////////////////////////////////////////////
// compute the inverse of a general 4x4 matrix using Cramer's Rule
// If cannot find inverse, return indentity matrix
// M^-1 = adj(M) / det(M)
///////////////////////////////////////////////////////////////////////////////
template<typename T>
Matrix4x4<T>& Matrix4x4<T>::invertGeneral()
{
    // get cofactors of minor matrices
    T cofactor0 = getCofactor(m[5],m[6],m[7], m[9],m[10],m[11], m[13],m[14],m[15]);
    T cofactor1 = getCofactor(m[4],m[6],m[7], m[8],m[10],m[11], m[12],m[14],m[15]);
    T cofactor2 = getCofactor(m[4],m[5],m[7], m[8],m[9], m[11], m[12],m[13],m[15]);
    T cofactor3 = getCofactor(m[4],m[5],m[6], m[8],m[9], m[10], m[12],m[13],m[14]);

    // get determinant
    T determinant = m[0] * cofactor0 - m[1] * cofactor1 + m[2] * cofactor2 - m[3] * cofactor3;
    if(MathUtils::isNullWithEpsilon(fabs(determinant)))
    {
        return identity();
    }

    // get rest of cofactors for adj(M)
    T cofactor4 = getCofactor(m[1],m[2],m[3], m[9],m[10],m[11], m[13],m[14],m[15]);
    T cofactor5 = getCofactor(m[0],m[2],m[3], m[8],m[10],m[11], m[12],m[14],m[15]);
    T cofactor6 = getCofactor(m[0],m[1],m[3], m[8],m[9], m[11], m[12],m[13],m[15]);
    T cofactor7 = getCofactor(m[0],m[1],m[2], m[8],m[9], m[10], m[12],m[13],m[14]);

    T cofactor8 = getCofactor(m[1],m[2],m[3], m[5],m[6], m[7],  m[13],m[14],m[15]);
    T cofactor9 = getCofactor(m[0],m[2],m[3], m[4],m[6], m[7],  m[12],m[14],m[15]);
    T cofactor10= getCofactor(m[0],m[1],m[3], m[4],m[5], m[7],  m[12],m[13],m[15]);
    T cofactor11= getCofactor(m[0],m[1],m[2], m[4],m[5], m[6],  m[12],m[13],m[14]);

    T cofactor12= getCofactor(m[1],m[2],m[3], m[5],m[6], m[7],  m[9], m[10],m[11]);
    T cofactor13= getCofactor(m[0],m[2],m[3], m[4],m[6], m[7],  m[8], m[10],m[11]);
    T cofactor14= getCofactor(m[0],m[1],m[3], m[4],m[5], m[7],  m[8], m[9], m[11]);
    T cofactor15= getCofactor(m[0],m[1],m[2], m[4],m[5], m[6],  m[8], m[9], m[10]);

    // build inverse matrix = adj(M) / det(M)
    // adjugate of M is the transpose of the cofactor matrix of M
    T invDeterminant = 1.0f / determinant;
    m[0] =  invDeterminant * cofactor0;
    m[1] = -invDeterminant * cofactor4;
    m[2] =  invDeterminant * cofactor8;
    m[3] = -invDeterminant * cofactor12;

    m[4] = -invDeterminant * cofactor1;
    m[5] =  invDeterminant * cofactor5;
    m[6] = -invDeterminant * cofactor9;
    m[7] =  invDeterminant * cofactor13;

    m[8] =  invDeterminant * cofactor2;
    m[9] = -invDeterminant * cofactor6;
    m[10]=  invDeterminant * cofactor10;
    m[11]= -invDeterminant * cofactor14;

    m[12]= -invDeterminant * cofactor3;
    m[13]=  invDeterminant * cofactor7;
    m[14]= -invDeterminant * cofactor11;
    m[15]=  invDeterminant * cofactor15;

    return *this;
}

///////////////////////////////////////////////////////////////////////////////
// translate this matrix by (x, y, z)
///////////////////////////////////////////////////////////////////////////////
template<typename T>
Matrix4x4<T>& Matrix4x4<T>::translate(const Vector<T, 3>& v)
{
    return translate(v.x, v.y, v.z);
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::translate(T x, T y, T z)
{
    m[0] += m[3] * x;   m[4] += m[7] * x;   m[8] += m[11]* x;   m[12]+= m[15]* x;
    m[1] += m[3] * y;   m[5] += m[7] * y;   m[9] += m[11]* y;   m[13]+= m[15]* y;
    m[2] += m[3] * z;   m[6] += m[7] * z;   m[10]+= m[11]* z;   m[14]+= m[15]* z;

    return *this;
}

///////////////////////////////////////////////////////////////////////////////
// uniform scale
///////////////////////////////////////////////////////////////////////////////
template<typename T>
Matrix4x4<T>& Matrix4x4<T>::scale(T s)
{
    return scale(s, s, s);
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::scale(T x, T y, T z)
{
    m[0] *= x;   m[4] *= x;   m[8] *= x;   m[12] *= x;
    m[1] *= y;   m[5] *= y;   m[9] *= y;   m[13] *= y;
    m[2] *= z;   m[6] *= z;   m[10]*= z;   m[14] *= z;
    return *this;
}

///////////////////////////////////////////////////////////////////////////////
// build a rotation matrix with given angle(degree) and rotation axis, then
// multiply it with this object
///////////////////////////////////////////////////////////////////////////////
template<typename T>
Matrix4x4<T>& Matrix4x4<T>::rotate(T angle, const Vector<T, 3>& axis)
{
    return rotate(angle, axis.x, axis.y, axis.z);
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::rotate(T angle, T x, T y, T z)
{
    T c = cosf(angle * Math::MathUtils::DEG_TO_RAD);    // cosine
    T s = sinf(angle * Math::MathUtils::DEG_TO_RAD);    // sine
    T c1 = 1.0f - c;                // 1 - c
    T m0 = m[0],  m4 = m[4],  m8 = m[8],  m12= m[12],
            m1 = m[1],  m5 = m[5],  m9 = m[9],  m13= m[13],
            m2 = m[2],  m6 = m[6],  m10= m[10], m14= m[14];

    // build rotation matrix
    T r0 = x * x * c1 + c;
    T r1 = x * y * c1 + z * s;
    T r2 = x * z * c1 - y * s;
    T r4 = x * y * c1 - z * s;
    T r5 = y * y * c1 + c;
    T r6 = y * z * c1 + x * s;
    T r8 = x * z * c1 + y * s;
    T r9 = y * z * c1 - x * s;
    T r10= z * z * c1 + c;

    // multiply rotation matrix
    m[0] = r0 * m0 + r4 * m1 + r8 * m2;
    m[1] = r1 * m0 + r5 * m1 + r9 * m2;
    m[2] = r2 * m0 + r6 * m1 + r10* m2;
    m[4] = r0 * m4 + r4 * m5 + r8 * m6;
    m[5] = r1 * m4 + r5 * m5 + r9 * m6;
    m[6] = r2 * m4 + r6 * m5 + r10* m6;
    m[8] = r0 * m8 + r4 * m9 + r8 * m10;
    m[9] = r1 * m8 + r5 * m9 + r9 * m10;
    m[10]= r2 * m8 + r6 * m9 + r10* m10;
    m[12]= r0 * m12+ r4 * m13+ r8 * m14;
    m[13]= r1 * m12+ r5 * m13+ r9 * m14;
    m[14]= r2 * m12+ r6 * m13+ r10* m14;

    return *this;
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::rotateX(T angle)
{
    T c = cosf(angle * Math::MathUtils::DEG_TO_RAD);
    T s = sinf(angle * Math::MathUtils::DEG_TO_RAD);
    T m1 = m[1],  m2 = m[2],
            m5 = m[5],  m6 = m[6],
            m9 = m[9],  m10= m[10],
            m13= m[13], m14= m[14];

    m[1] = m1 * c + m2 *-s;
    m[2] = m1 * s + m2 * c;
    m[5] = m5 * c + m6 *-s;
    m[6] = m5 * s + m6 * c;
    m[9] = m9 * c + m10*-s;
    m[10]= m9 * s + m10* c;
    m[13]= m13* c + m14*-s;
    m[14]= m13* s + m14* c;

    return *this;
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::rotateY(T angle)
{
    T c = cosf(angle * Math::MathUtils::DEG_TO_RAD);
    T s = sinf(angle * Math::MathUtils::DEG_TO_RAD);
    T m0 = m[0],  m2 = m[2],
            m4 = m[4],  m6 = m[6],
            m8 = m[8],  m10= m[10],
            m12= m[12], m14= m[14];

    m[0] = m0 * c + m2 * s;
    m[2] = m0 *-s + m2 * c;
    m[4] = m4 * c + m6 * s;
    m[6] = m4 *-s + m6 * c;
    m[8] = m8 * c + m10* s;
    m[10]= m8 *-s + m10* c;
    m[12]= m12* c + m14* s;
    m[14]= m12*-s + m14* c;

    return *this;
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::rotateZ(T angle)
{
    T c = cosf(angle * Math::MathUtils::DEG_TO_RAD);
    T s = sinf(angle * Math::MathUtils::DEG_TO_RAD);
    T m0 = m[0],  m1 = m[1],
            m4 = m[4],  m5 = m[5],
            m8 = m[8],  m9 = m[9],
            m12= m[12], m13= m[13];

    m[0] = m0 * c + m1 *-s;
    m[1] = m0 * s + m1 * c;
    m[4] = m4 * c + m5 *-s;
    m[5] = m4 * s + m5 * c;
    m[8] = m8 * c + m9 *-s;
    m[9] = m8 * s + m9 * c;
    m[12]= m12* c + m13*-s;
    m[13]= m12* s + m13* c;

    return *this;
}

template<typename T>
Matrix4x4<T> Matrix4x4<T>::operator+(const Matrix4x4<T>& rhs) const
{
    return Matrix4x4<T>(m[0]+rhs[0],   m[1]+rhs[1],   m[2]+rhs[2],   m[3]+rhs[3],
                   m[4]+rhs[4],   m[5]+rhs[5],   m[6]+rhs[6],   m[7]+rhs[7],
                   m[8]+rhs[8],   m[9]+rhs[9],   m[10]+rhs[10], m[11]+rhs[11],
                   m[12]+rhs[12], m[13]+rhs[13], m[14]+rhs[14], m[15]+rhs[15]);
}

template<typename T>
Matrix4x4<T> Matrix4x4<T>::operator-(const Matrix4x4<T>& rhs) const
{
    return Matrix4x4<T>(m[0]-rhs[0],   m[1]-rhs[1],   m[2]-rhs[2],   m[3]-rhs[3],
                   m[4]-rhs[4],   m[5]-rhs[5],   m[6]-rhs[6],   m[7]-rhs[7],
                   m[8]-rhs[8],   m[9]-rhs[9],   m[10]-rhs[10], m[11]-rhs[11],
                   m[12]-rhs[12], m[13]-rhs[13], m[14]-rhs[14], m[15]-rhs[15]);
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::operator+=(const Matrix4x4<T>& rhs)
{
    m[0] += rhs[0];   m[1] += rhs[1];   m[2] += rhs[2];   m[3] += rhs[3];
    m[4] += rhs[4];   m[5] += rhs[5];   m[6] += rhs[6];   m[7] += rhs[7];
    m[8] += rhs[8];   m[9] += rhs[9];   m[10]+= rhs[10];  m[11]+= rhs[11];
    m[12]+= rhs[12];  m[13]+= rhs[13];  m[14]+= rhs[14];  m[15]+= rhs[15];
    return *this;
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::operator-=(const Matrix4x4<T>& rhs)
{
    m[0] -= rhs[0];   m[1] -= rhs[1];   m[2] -= rhs[2];   m[3] -= rhs[3];
    m[4] -= rhs[4];   m[5] -= rhs[5];   m[6] -= rhs[6];   m[7] -= rhs[7];
    m[8] -= rhs[8];   m[9] -= rhs[9];   m[10]-= rhs[10];  m[11]-= rhs[11];
    m[12]-= rhs[12];  m[13]-= rhs[13];  m[14]-= rhs[14];  m[15]-= rhs[15];
    return *this;
}

template<typename T>
Vector<T, 4> Matrix4x4<T>::operator*(const Vector<T, 4>& rhs) const
{
    return Vector<T, 4>(m[0]*rhs.x + m[4]*rhs.y + m[8]*rhs.z  + m[12]*rhs.w,
                   m[1]*rhs.x + m[5]*rhs.y + m[9]*rhs.z  + m[13]*rhs.w,
                   m[2]*rhs.x + m[6]*rhs.y + m[10]*rhs.z + m[14]*rhs.w,
                   m[3]*rhs.x + m[7]*rhs.y + m[11]*rhs.z + m[15]*rhs.w);
}

template<typename T>
Vector<T, 3> Matrix4x4<T>::operator*(const Vector<T, 3>& rhs) const
{
    return Vector<T, 3>(m[0]*rhs.x + m[4]*rhs.y + m[8]*rhs.z,
                   m[1]*rhs.x + m[5]*rhs.y + m[9]*rhs.z,
                   m[2]*rhs.x + m[6]*rhs.y + m[10]*rhs.z);
}

template<typename T>
Matrix4x4<T> Matrix4x4<T>::operator*(const Matrix4x4<T>& n) const
{
    return Matrix4x4<T>(m[0]*n[0]  + m[4]*n[1]  + m[8]*n[2]  + m[12]*n[3],   m[1]*n[0]  + m[5]*n[1]  + m[9]*n[2]  + m[13]*n[3],   m[2]*n[0]  + m[6]*n[1]  + m[10]*n[2]  + m[14]*n[3],   m[3]*n[0]  + m[7]*n[1]  + m[11]*n[2]  + m[15]*n[3],
                   m[0]*n[4]  + m[4]*n[5]  + m[8]*n[6]  + m[12]*n[7],   m[1]*n[4]  + m[5]*n[5]  + m[9]*n[6]  + m[13]*n[7],   m[2]*n[4]  + m[6]*n[5]  + m[10]*n[6]  + m[14]*n[7],   m[3]*n[4]  + m[7]*n[5]  + m[11]*n[6]  + m[15]*n[7],
                   m[0]*n[8]  + m[4]*n[9]  + m[8]*n[10] + m[12]*n[11],  m[1]*n[8]  + m[5]*n[9]  + m[9]*n[10] + m[13]*n[11],  m[2]*n[8]  + m[6]*n[9]  + m[10]*n[10] + m[14]*n[11],  m[3]*n[8]  + m[7]*n[9]  + m[11]*n[10] + m[15]*n[11],
                   m[0]*n[12] + m[4]*n[13] + m[8]*n[14] + m[12]*n[15],  m[1]*n[12] + m[5]*n[13] + m[9]*n[14] + m[13]*n[15],  m[2]*n[12] + m[6]*n[13] + m[10]*n[14] + m[14]*n[15],  m[3]*n[12] + m[7]*n[13] + m[11]*n[14] + m[15]*n[15]);
}

template<typename T>
Matrix4x4<T>& Matrix4x4<T>::operator*=(const Matrix4x4<T>& rhs)
{
    *this = *this * rhs;
    return *this;
}

template<typename T>
bool Matrix4x4<T>::operator==(const Matrix4x4<T>& n) const
{
    return (m[0] == n[0])  && (m[1] == n[1])  && (m[2] == n[2])  && (m[3] == n[3])  &&
           (m[4] == n[4])  && (m[5] == n[5])  && (m[6] == n[6])  && (m[7] == n[7])  &&
           (m[8] == n[8])  && (m[9] == n[9])  && (m[10]== n[10]) && (m[11]== n[11]) &&
           (m[12]== n[12]) && (m[13]== n[13]) && (m[14]== n[14]) && (m[15]== n[15]);
}

template<typename T>
bool Matrix4x4<T>::operator!=(const Matrix4x4<T>& n) const
{
    return (m[0] != n[0])  || (m[1] != n[1])  || (m[2] != n[2])  || (m[3] != n[3])  ||
           (m[4] != n[4])  || (m[5] != n[5])  || (m[6] != n[6])  || (m[7] != n[7])  ||
           (m[8] != n[8])  || (m[9] != n[9])  || (m[10]!= n[10]) || (m[11]!= n[11]) ||
           (m[12]!= n[12]) || (m[13]!= n[13]) || (m[14]!= n[14]) || (m[15]!= n[15]);
}

template<typename T>
T Matrix4x4<T>::operator[](int index) const
{
    return m[index];
}

template<typename T>
T& Matrix4x4<T>::operator[](int index)
{
    return m[index];
}

template<typename T>
Matrix4x4<T> operator-(const Matrix4x4<T>& rhs)
{
    return Matrix4x4<T>(-rhs[0], -rhs[1], -rhs[2], -rhs[3], -rhs[4], -rhs[5], -rhs[6], -rhs[7], -rhs[8], -rhs[9], -rhs[10], -rhs[11], -rhs[12], -rhs[13], -rhs[14], -rhs[15]);
}

template<typename T>
Matrix4x4<T> operator*(T s, const Matrix4x4<T>& rhs)
{
    return Matrix4x4<T>(s*rhs[0], s*rhs[1], s*rhs[2], s*rhs[3], s*rhs[4], s*rhs[5], s*rhs[6], s*rhs[7], s*rhs[8], s*rhs[9], s*rhs[10], s*rhs[11], s*rhs[12], s*rhs[13], s*rhs[14], s*rhs[15]);
}

template<typename T>
Vector<T, 4> operator*(const Vector<T, 4>& v, const Matrix4x4<T>& m)
{
    return Vector<T, 4>(v.x*m[0] + v.y*m[1] + v.z*m[2] + v.w*m[3],  v.x*m[4] + v.y*m[5] + v.z*m[6] + v.w*m[7],  v.x*m[8] + v.y*m[9] + v.z*m[10] + v.w*m[11], v.x*m[12] + v.y*m[13] + v.z*m[14] + v.w*m[15]);
}

template<typename T>
Vector<T, 3> operator*(const Vector<T, 3>& v, const Matrix4x4<T>& m)
{
    return Vector<T, 3>(v.x*m[0] + v.y*m[1] + v.z*m[2],  v.x*m[4] + v.y*m[5] + v.z*m[6],  v.x*m[8] + v.y*m[9] + v.z*m[10]);
}

template<typename T>
::std::ostream& operator<<(::std::ostream& os, const Matrix4x4<T>& m)
{
    os << ::std::fixed << ::std::setprecision(5);
    os << "[" << ::std::setw(10) << m[0] << " " << ::std::setw(10) << m[4] << " " << ::std::setw(10) << m[8]  <<  " " << ::std::setw(10) << m[12] << "]\n"
    << "[" << ::std::setw(10) << m[1] << " " << ::std::setw(10) << m[5] << " " << ::std::setw(10) << m[9]  <<  " " << ::std::setw(10) << m[13] << "]\n"
    << "[" << ::std::setw(10) << m[2] << " " << ::std::setw(10) << m[6] << " " << ::std::setw(10) << m[10] <<  " " << ::std::setw(10) << m[14] << "]\n"
    << "[" << ::std::setw(10) << m[3] << " " << ::std::setw(10) << m[7] << " " << ::std::setw(10) << m[11] <<  " " << ::std::setw(10) << m[15] << "]\n";
    os << ::std::resetiosflags(::std::ios_base::fixed | ::std::ios_base::floatfield);
    return os;
}

///////////////////////////////////////////////////////////////////////////////
// compute cofactor of 3x3 minor matrix without sign
// input params are 9 elements of the minor matrix
// NOTE: The caller must know its sign.
///////////////////////////////////////////////////////////////////////////////
template<typename T>
T Matrix4x4<T>::getCofactor(T m0, T m1, T m2,
                           T m3, T m4, T m5,
                           T m6, T m7, T m8)
{
    return m0 * (m4 * m8 - m5 * m7) -
           m1 * (m3 * m8 - m5 * m6) +
           m2 * (m3 * m7 - m4 * m6);
}

} //namespace Math