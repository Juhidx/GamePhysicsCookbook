﻿#include "matrices.h"
#include <cmath>
#include <cfloat>

#define CMP(x, y)                       \
(fabsf((x)+(y)) <= FLT_EPSILON *    \
fmaxf(1.0f,                     \
fmaxf(fabsf(x), fabsf(y)))      \
)

mat2 operator*(const mat2& matrix, const float scalar)
{
    mat2 result;
    for (int i = 0; i < 4; ++i)
    {
        result.asArray[i] = matrix.asArray[i] * scalar;
    }
    return result;
}

mat3 operator*(const mat3& matrix, const float scalar)
{
    mat3 result;
    for (int i = 0; i < 9; ++i)
    {
        result.asArray[i] = matrix.asArray[i] * scalar;
    }
    return result;
}

mat4 operator*(const mat4& matrix, const float scalar)
{
    mat4 result;
    for (int i = 0; i < 16; ++i)
    {
        result.asArray[i] = matrix.asArray[i] * scalar;
    }
    return result;
}

mat2 operator*(const mat2& matA, const mat2& matB)
{
    mat2 res;
    Multiply(res.asArray, matA.asArray, 2, 2, matB.asArray, 2, 2);
    return res;
}

mat3 operator*(const mat3& matA, const mat3& matB)
{
    mat3 res;
    Multiply(res.asArray, matA.asArray, 3, 3, matB.asArray, 3, 3);
    return res;
}

mat4 operator*(const mat4& matA, const mat4& matB)
{
    mat4 res;
    Multiply(res.asArray, matA.asArray, 4, 4, matB.asArray, 4, 4);
    return res;
}


void Transpose(const float* srcMat, float* dstMat, int srcRows, int srcCols)
{
    int row;
    int col;
    for (int i = 0; i < srcRows * srcCols; ++i)
    {
        row = i / srcRows;
        col = i % srcRows;
        dstMat[i] = srcMat[srcCols *  col + row];
    }
}

mat2 Transpose(const mat2& matrix)
{
    mat2 result;
    Transpose(matrix.asArray, result.asArray, 2, 2);
    return result;
}

mat3 Transpose(const mat3& matrix)
{
    mat3 result;
    Transpose(matrix.asArray, result.asArray, 3, 3);
    return result;
}

mat4 Transpose(const mat4& matrix)
{
    mat4 result;
    Transpose(matrix.asArray, result.asArray, 4, 4);
    return result;
}

bool Multiply(float* out, const float* matA, int aRows, int aCols, const float* matB, int bRows, int bCols)
{
    if(aCols != bRows)
    {
        return false;
    }

    int a;
    int b;

    for (int i = 0; i < aRows; ++i)
    {
        for (int j = 0; j < bCols; ++j)
        {
            out[bCols * i + j] = 0.0f;
            for (int k = 0; k < bRows; ++k)
            {
                a = aCols * i + k;
                b = bCols * k + j;
                out[bCols * i + j] += matA[a] * matB[b];
            }
        }
    }

    return true;
}