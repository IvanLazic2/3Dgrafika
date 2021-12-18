#pragma once
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <cstdlib> 
#include <cstdio> 
#include <iomanip> 

// glavni predlozak koji je parametriziran s dimnezijom i tipom
template <size_t DIM, typename T>
struct vec // definiraj vektor
{
    // konstruktor
    vec()
    {
        for (size_t i = DIM; i--; data_[i] = T())
            ;
    }
    // operator pristupa
    T &operator[](const size_t i)
    {
        assert(i < DIM);
        return data_[i];
    }
    const T &operator[](const size_t i) const
    {
        assert(i < DIM);
        return data_[i];
    }

private:
    T data_[DIM];
};

// specijalizacija tipa ovisno o float ili int i po velicine dimenzije
typedef vec<2, float> Vec2f;
typedef vec<3, float> Vec3f;
typedef vec<3, int> Vec3d;
typedef vec<4, float> Vec4f;

template <typename T>
struct vec<2, T>
{
    // konstruktori
    vec() : x(T()), y(T()) {}
    vec(T X, T Y) : x(X), y(Y) {}
    template <class U>
    vec<2, T>(const vec<2, U> &v);
    // operator pristupa
    T &operator[](const size_t i)
    {
        assert(i < 2);
        return i <= 0 ? x : y;
    }
    const T &operator[](const size_t i) const
    {
        assert(i < 2);
        return i <= 0 ? x : y;
    }
    // komponente
    T x, y;
};

template <typename T>
struct vec<3, T>
{
    // konstruktori
    vec() : x(T()), y(T()), z(T()) {}
    vec(T X, T Y, T Z) : x(X), y(Y), z(Z) {}
    // operatori pristupa
    T &operator[](const size_t i)
    {
        assert(i < 3);
        return i <= 0 ? x : (1 == i ? y : z);
    }
    const T &operator[](const size_t i) const
    {
        assert(i < 3);
        return i <= 0 ? x : (1 == i ? y : z);
    }

    // normalizacija
    float norm() { return std::sqrt(x * x + y * y + z * z); }
    float norm() const { return std::sqrt(x * x + y * y + z * z); }
    vec<3, T> &normalize(T l = 1)
    {
        *this = (*this) * (l / norm());
        return *this;
    }
    // komponente
    T x, y, z;
};

template <typename T>
struct vec<4, T>
{
    // konstruktori
    vec() : x(T()), y(T()), z(T()), w(T()) {}
    vec(T X, T Y, T Z, T W) : x(X), y(Y), z(Z), w(W) {}
    // operator pristupa
    T &operator[](const size_t i)
    {
        assert(i < 4);
        return i <= 0 ? x : (1 == i ? y : (2 == i ? z : w));
    }
    const T &operator[](const size_t i) const
    {
        assert(i < 4);
        return i <= 0 ? x : (1 == i ? y : (2 == i ? z : w));
    }

    // komponente
    T x, y, z, w;
};

// implementacija skalarnog produkta
template <size_t DIM, typename T>
T operator*(const vec<DIM, T> &lhs, const vec<DIM, T> &rhs)
{
    T ret = T();
    for (size_t i = DIM; i--; ret += lhs[i] * rhs[i])
        ;
    return ret;
}

// implementacija zbrajanja vektora
template <size_t DIM, typename T>
vec<DIM, T> operator+(vec<DIM, T> lhs, const vec<DIM, T> &rhs)
{
    for (size_t i = DIM; i--; lhs[i] += rhs[i])
        ;
    return lhs;
}
// implementacija oduzimanja vektora rhs
template <size_t DIM, typename T>
vec<DIM, T> operator-(vec<DIM, T> lhs, const vec<DIM, T> &rhs)
{
    for (size_t i = DIM; i--; lhs[i] -= rhs[i])
        ;
    return lhs;
}

// implementacija produkta po kompoenentama
template <size_t DIM, typename T, typename U>
vec<DIM, T> operator*(const vec<DIM, T> &lhs, const U &rhs)
{
    vec<DIM, T> ret;
    for (size_t i = DIM; i--; ret[i] = lhs[i] * rhs)
        ;
    return ret;
}

template <size_t DIM, typename T, typename U>
vec<DIM, T> operator/(const vec<DIM, T> &lhs, const U &rhs)
{
    vec<DIM, T> ret;
    for (size_t i = DIM; i--; ret[i] = lhs[i] / rhs)
        ;
    return ret;
}

// implementacija oduzimanja vektora lhs
template <size_t DIM, typename T>
vec<DIM, T> operator-(const vec<DIM, T> &lhs)
{
    return lhs * T(-1);
}

// vektorski produkt
template <typename T>
vec<3, T> cross(vec<3, T> v1, vec<3, T> v2)
{
    return vec<3, T>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

// determinanata za 3x3 matrice
template <typename T>
float determinant(const vec<3, T> &v1, const vec<3, T> &v2, const vec<3, T> &v3)
{
    return (v1.x * (v2.y * v3.z - v2.z * v3.y) - v2.x * (v1.y * v3.z - v1.z * v3.y) + v3.x * (v1.y * v2.z - v1.z * v2.y));
}

// ispis vektora
template <size_t DIM, typename T>
std::ostream &operator<<(std::ostream &out, const vec<DIM, T> &v)
{
    for (unsigned int i = 0; i < DIM; i++)
    {
        out << v[i] << " ";
    }
    return out;
}

//Matrix 4x4

template <typename T>
class Matrix44
{
public:
    T x[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

    Matrix44() {}

    Matrix44(T a, T b, T c, T d, T e, T f, T g, T h,
             T i, T j, T k, T l, T m, T n, T o, T p)
    {
        x[0][0] = a;
        x[0][1] = b;
        x[0][2] = c;
        x[0][3] = d;
        x[1][0] = e;
        x[1][1] = f;
        x[1][2] = g;
        x[1][3] = h;
        x[2][0] = i;
        x[2][1] = j;
        x[2][2] = k;
        x[2][3] = l;
        x[3][0] = m;
        x[3][1] = n;
        x[3][2] = o;
        x[3][3] = p;
    }

    const T *operator[](uint8_t i) const { return x[i]; }
    T *operator[](uint8_t i) { return x[i]; }

    // Multiply the current matrix with another matrix (rhs)
    Matrix44 operator*(const Matrix44 &v) const
    {
        Matrix44 tmp;
        multiply(*this, v, tmp);

        return tmp;
    }

    static void multiply(const Matrix44<T> &a, const Matrix44 &b, Matrix44 &c)
    {
#if 0 
        for (uint8_t i = 0; i < 4; ++i) { 
            for (uint8_t j = 0; j < 4; ++j) { 
                c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + 
                    a[i][2] * b[2][j] + a[i][3] * b[3][j]; 
            } 
        }
#else
        // A restric qualified pointer (or reference) is basically a promise
        // to the compiler that for the scope of the pointer, the target of the
        // pointer will only be accessed through that pointer (and pointers
        // copied from it.
        const T *__restrict ap = &a.x[0][0];
        const T *__restrict bp = &b.x[0][0];
        T *__restrict cp = &c.x[0][0];

        T a0, a1, a2, a3;

        a0 = ap[0];
        a1 = ap[1];
        a2 = ap[2];
        a3 = ap[3];

        cp[0] = a0 * bp[0] + a1 * bp[4] + a2 * bp[8] + a3 * bp[12];
        cp[1] = a0 * bp[1] + a1 * bp[5] + a2 * bp[9] + a3 * bp[13];
        cp[2] = a0 * bp[2] + a1 * bp[6] + a2 * bp[10] + a3 * bp[14];
        cp[3] = a0 * bp[3] + a1 * bp[7] + a2 * bp[11] + a3 * bp[15];

        a0 = ap[4];
        a1 = ap[5];
        a2 = ap[6];
        a3 = ap[7];

        cp[4] = a0 * bp[0] + a1 * bp[4] + a2 * bp[8] + a3 * bp[12];
        cp[5] = a0 * bp[1] + a1 * bp[5] + a2 * bp[9] + a3 * bp[13];
        cp[6] = a0 * bp[2] + a1 * bp[6] + a2 * bp[10] + a3 * bp[14];
        cp[7] = a0 * bp[3] + a1 * bp[7] + a2 * bp[11] + a3 * bp[15];

        a0 = ap[8];
        a1 = ap[9];
        a2 = ap[10];
        a3 = ap[11];

        cp[8] = a0 * bp[0] + a1 * bp[4] + a2 * bp[8] + a3 * bp[12];
        cp[9] = a0 * bp[1] + a1 * bp[5] + a2 * bp[9] + a3 * bp[13];
        cp[10] = a0 * bp[2] + a1 * bp[6] + a2 * bp[10] + a3 * bp[14];
        cp[11] = a0 * bp[3] + a1 * bp[7] + a2 * bp[11] + a3 * bp[15];

        a0 = ap[12];
        a1 = ap[13];
        a2 = ap[14];
        a3 = ap[15];

        cp[12] = a0 * bp[0] + a1 * bp[4] + a2 * bp[8] + a3 * bp[12];
        cp[13] = a0 * bp[1] + a1 * bp[5] + a2 * bp[9] + a3 * bp[13];
        cp[14] = a0 * bp[2] + a1 * bp[6] + a2 * bp[10] + a3 * bp[14];
        cp[15] = a0 * bp[3] + a1 * bp[7] + a2 * bp[11] + a3 * bp[15];
#endif
    }

    // \brief return a transposed copy of the current matrix as a new matrix
    Matrix44 transposed() const
    {
#if 0 
        Matrix44 t; 
        for (uint8_t i = 0; i < 4; ++i) { 
            for (uint8_t j = 0; j < 4; ++j) { 
                t[i][j] = x[j][i]; 
            } 
        } 
 
        return t;
#else
        return Matrix44(x[0][0],
                        x[1][0],
                        x[2][0],
                        x[3][0],
                        x[0][1],
                        x[1][1],
                        x[2][1],
                        x[3][1],
                        x[0][2],
                        x[1][2],
                        x[2][2],
                        x[3][2],
                        x[0][3],
                        x[1][3],
                        x[2][3],
                        x[3][3]);
#endif
    }

    // \brief transpose itself
    Matrix44 &transpose()
    {
        Matrix44 tmp(x[0][0],
                     x[1][0],
                     x[2][0],
                     x[3][0],
                     x[0][1],
                     x[1][1],
                     x[2][1],
                     x[3][1],
                     x[0][2],
                     x[1][2],
                     x[2][2],
                     x[3][2],
                     x[0][3],
                     x[1][3],
                     x[2][3],
                     x[3][3]);
        *this = tmp;

        return *this;
    }

    template <typename S>
    void multVecMatrix(const vec<3, S> &src, vec<3, S> &dst) const
    {
        S a, b, c, w;

        a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0] + x[3][0];
        b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1] + x[3][1];
        c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2] + x[3][2];
        w = src[0] * x[0][3] + src[1] * x[1][3] + src[2] * x[2][3] + x[3][3];

        dst.x = a / w;
        dst.y = b / w;
        dst.z = c / w;
    }

    template<typename S> 
    void multDirMatrix(const vec<3, S> &src, vec<3, S> &dst) const 
    { 
        S a, b, c; 
 
        a = src[0] * x[0][0] + src[1] * x[1][0] + src[2] * x[2][0]; 
        b = src[0] * x[0][1] + src[1] * x[1][1] + src[2] * x[2][1]; 
        c = src[0] * x[0][2] + src[1] * x[1][2] + src[2] * x[2][2]; 
 
        dst.x = a; 
        dst.y = b; 
        dst.z = c; 
    } 

    Matrix44 inverse() const
    {
        int i, j, k;
        Matrix44 s;
        Matrix44 t(*this);

        // Forward elimination
        for (i = 0; i < 3; i++)
        {
            int pivot = i;

            T pivotsize = t[i][i];

            if (pivotsize < 0)
                pivotsize = -pivotsize;

            for (j = i + 1; j < 4; j++)
            {
                T tmp = t[j][i];

                if (tmp < 0)
                    tmp = -tmp;

                if (tmp > pivotsize)
                {
                    pivot = j;
                    pivotsize = tmp;
                }
            }

            if (pivotsize == 0)
            {
                // Cannot invert singular matrix
                return Matrix44();
            }

            if (pivot != i)
            {
                for (j = 0; j < 4; j++)
                {
                    T tmp;

                    tmp = t[i][j];
                    t[i][j] = t[pivot][j];
                    t[pivot][j] = tmp;

                    tmp = s[i][j];
                    s[i][j] = s[pivot][j];
                    s[pivot][j] = tmp;
                }
            }

            for (j = i + 1; j < 4; j++)
            {
                T f = t[j][i] / t[i][i];

                for (k = 0; k < 4; k++)
                {
                    t[j][k] -= f * t[i][k];
                    s[j][k] -= f * s[i][k];
                }
            }
        }

        // Backward substitution
        for (i = 3; i >= 0; --i)
        {
            T f;

            if ((f = t[i][i]) == 0)
            {
                // Cannot invert singular matrix
                return Matrix44();
            }

            for (j = 0; j < 4; j++)
            {
                t[i][j] /= f;
                s[i][j] /= f;
            }

            for (j = 0; j < i; j++)
            {
                f = t[j][i];

                for (k = 0; k < 4; k++)
                {
                    t[j][k] -= f * t[i][k];
                    s[j][k] -= f * s[i][k];
                }
            }
        }

        return s;
    }

    // \brief set current matrix to its inverse
    const Matrix44<T> &invert()
    {
        *this = inverse();
        return *this;
    }

    friend std::ostream &operator<<(std::ostream &s, const Matrix44 &m)
    {
        std::ios_base::fmtflags oldFlags = s.flags();
        int width = 12; // total with of the displayed number
        s.precision(5); // control the number of displayed decimals
        s.setf(std::ios_base::fixed);

        s << "[" << std::setw(width) << m[0][0] << " " << std::setw(width) << m[0][1] << " " << std::setw(width) << m[0][2] << " " << std::setw(width) << m[0][3] << "\n"
          <<

            " " << std::setw(width) << m[1][0] << " " << std::setw(width) << m[1][1] << " " << std::setw(width) << m[1][2] << " " << std::setw(width) << m[1][3] << "\n"
          <<

            " " << std::setw(width) << m[2][0] << " " << std::setw(width) << m[2][1] << " " << std::setw(width) << m[2][2] << " " << std::setw(width) << m[2][3] << "\n"
          <<

            " " << std::setw(width) << m[3][0] << " " << std::setw(width) << m[3][1] << " " << std::setw(width) << m[3][2] << " " << std::setw(width) << m[3][3] << "]";

        s.flags(oldFlags);
        return s;
    }
};

typedef Matrix44<float> Matrix44f;