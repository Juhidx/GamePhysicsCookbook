#include "matrices.h"
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

mat2 Transpose(const mat2& mat)
{
    mat2 result;
    Transpose(mat.asArray, result.asArray, 2, 2);
    return result;
}

mat3 Transpose(const mat3& mat)
{
    mat3 result;
    Transpose(mat.asArray, result.asArray, 3, 3);
    return result;
}

mat4 Transpose(const mat4& mat)
{
    mat4 result;
    Transpose(mat.asArray, result.asArray, 4, 4);
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

float Determinant(const mat2& mat)
{
    return mat._11 * mat._22 - mat._12 * mat._21;
}

float Determinant(const mat3& mat)
{
    float result = 0.0f;
    int index;
    mat3 cofactor = Cofactor(mat);

    for (int i = 0; i < 3; ++i)
    {
        index = 3 * 0 + i;
        result += mat.asArray[index] * cofactor[0][i];
    }
    
    return result;
}

float Determinant(const mat4& mat)
{
    float result = 0.0f;
    mat4 cofactor = Cofactor(mat);

    for (int i = 0; i < 4; ++i)
    {
        result += mat.asArray[4 * 0 + i] * cofactor[0][i];
    }

    return result;
}


mat2 Cut(const mat3& mat, int row, int col)
{
    mat2 result;
    int index = 0;
    int target;
    int source;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if(i == row || j == col)
            {
                continue;
            }

            target = index++;
            source = 3 * i + j;
            result.asArray[target] = mat.asArray[source];
        }
    }

    return result;
}

mat3 Cut(const mat4& mat, int row, int col)
{
    mat3 result;
    int index = 0;
    int target;
    int source;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            if(i == row || j == col)
            {
                continue;
            }

            target = index++;
            source = 4 * i + j;
            result.asArray[target] = mat.asArray[source];
        }
    }
}


mat2 Minor(const mat2& mat)
{
    return mat2(
        mat._22, mat._21,
        mat._12, mat._11);
}


mat3 Minor(const mat3& mat)
{
    mat3 result;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            result[i][j] = Determinant(Cut(mat, i, j));
        }
    }

    return result;
}

mat4 Minor(const mat4& mat)
{
    mat4 result;

    for(int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            result[i][j] = Determinant(Cut(mat, i , j));
        }
    }
    
    return result;
}


void Cofactor(float* out, const float* minor, int rows, int cols)
{
    int t;
    int s;
    float sign;
    
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            s = t = cols * j + i;
            out[t] = minor[s] * powf(-1.0f, i + j); // + or -
        }
    }
}

mat2 Cofactor(const mat2& mat)
{
    mat2 result;
    Cofactor(result.asArray, Minor(mat).asArray, 2, 2);
    return result;
}

mat3 Cofactor(const mat3& mat)
{
    mat3 result;
    Cofactor(result.asArray, Minor(mat).asArray, 3, 3);
    return result;
}

mat4 Cofactor(const mat4& mat)
{
    mat4 result;
    Cofactor(result.asArray, Minor(mat).asArray, 4, 4);
    return result;
}

mat2 Adjugate(const mat2& mat)
{
    return Transpose(Cofactor(mat));
}

mat3 Adjugate(const mat3& mat)
{
    return Transpose(Cofactor(mat));
}

mat4 Adjugate(const mat4& mat)
{
    return Transpose(Cofactor(mat));
}

mat2 Inverse(const mat2& mat)
{
    float det = mat._11 * mat._22 - mat._12 * mat._21;;
    
    if(CMP(det, 0.0f)) {return mat2();}
    mat2 result;
    float inv_det = 1.0f / det;
    
    result._11 = mat._22 * inv_det;
    result._12 = -mat._12 * inv_det;
    result._21 = -mat._21 * inv_det;
    result._22 = mat._11 * inv_det;
    
    return result;
}

mat3 Inverse(const mat3& mat)
{
    float det = mat._11 * (mat._22 * mat._33 - mat._23 * mat._32) -
                mat._12 * (mat._21 * mat._33 - mat._23 * mat._31) +
                mat._13 * (mat._21 * mat._32 - mat._22 * mat._31);
    
    if(CMP(det, 0.0f)) {return mat3();}

    mat3 result;
    float inv_det = 1.0f / det;

    result._11 = (mat._22 * mat._33 - mat._23 * mat._32) * inv_det;
    result._12 = -(mat._12 * mat._33 - mat._13 * mat._32) * inv_det;
    result._13 = (mat._12 * mat._23 - mat._13 * mat._22) * inv_det;
    result._21 = -(mat._21 * mat._33 - mat._23 * mat._31) * inv_det;
    result._22 = (mat._11 * mat._33 - mat._13 * mat._31) * inv_det;
    result._23 = -(mat._11 * mat._23 - mat._13 * mat._21) * inv_det;
    result._31 = (mat._21 * mat._32 - mat._22 * mat._31) * inv_det;
    result._32 = -(mat._11 * mat._32 - mat._12 * mat._31) * inv_det;
    result._33 = (mat._11 * mat._22 - mat._12 * mat._21) * inv_det;

    return result;
}

mat4 Inverse(const mat4& mat)
{
    float det = 
        mat._11 * (mat._22 * mat._33 * mat._44 - 
                    mat._23 * mat._34 * mat._42 + 
                    mat._24 * mat._32 * mat._43 - 
                    mat._24 * mat._33 * mat._42 + 
                    mat._23 * mat._34 * mat._41 - 
                    mat._22 * mat._33 * mat._44);

    if(CMP(det, 0.0f)) {return mat4();}

    float inv_det = 1.0f / det;
    
    mat4 result;
    
    // First row
    result._11 = (mat._22 * mat._33 * mat._44 - 
                  mat._23 * mat._34 * mat._42 + 
                  mat._24 * mat._32 * mat._43 - 
                  mat._24 * mat._33 * mat._42 + 
                  mat._23 * mat._34 * mat._41 - 
                  mat._22 * mat._33 * mat._44) * inv_det;
    result._12 = -(mat._21 * mat._33 * mat._44 - 
                  mat._23 * mat._34 * mat._41 + 
                  mat._24 * mat._31 * mat._43 - 
                  mat._24 * mat._33 * mat._41 + 
                  mat._21 * mat._34 * mat._42 - 
                  mat._22 * mat._33 * mat._44) * inv_det;
    result._13 = (mat._21 * mat._32 * mat._44 - 
                 mat._22 * mat._34 * mat._41 + 
                 mat._24 * mat._31 * mat._42 - 
                 mat._23 * mat._34 * mat._41 + 
                 mat._22 * mat._33 * mat._44 - 
                 mat._21 * mat._33 * mat._44) * inv_det;
    result._14 = -(mat._21 * mat._32 * mat._43 - 
                  mat._22 * mat._33 * mat._41 + 
                  mat._23 * mat._31 * mat._42 - 
                  mat._24 * mat._31 * mat._33 + 
                  mat._24 * mat._32 * mat._31 - 
                  mat._23 * mat._32 * mat._34) * inv_det;

    // Second row  
    result._21 = -(mat._12 * mat._33 * mat._44 - 
                  mat._13 * mat._34 * mat._41 + 
                  mat._14 * mat._31 * mat._43 - 
                  mat._14 * mat._33 * mat._41 + 
                  mat._12 * mat._34 * mat._31 - 
                  mat._11 * mat._33 * mat._44) * inv_det;
    result._22 = (mat._11 * mat._34 * mat._44 - 
                 mat._13 * mat._32 * mat._44 + 
                 mat._14 * mat._21 * mat._43 - 
                 mat._14 * mat._33 * mat._41 + 
                 mat._11 * mat._33 * mat._42 - 
                 mat._12 * mat._33 * mat._44) * inv_det;
    result._23 = -(mat._11 * mat._23 * mat._44 - 
                  mat._12 * mat._24 * mat._41 + 
                  mat._14 * mat._21 * mat._43 - 
                  mat._14 * mat._23 * mat._41 + 
                  mat._12 * mat._24 * mat._31 - 
                  mat._11 * mat._23 * mat._44) * inv_det;
    result._24 = (mat._11 * mat._23 * mat._34 - 
                 mat._12 * mat._24 * mat._31 + 
                 mat._13 * mat._21 * mat._33 - 
                 mat._13 * mat._23 * mat._31 + 
                 mat._12 * mat._24 * mat._21 - 
                 mat._11 * mat._23 * mat._34) * inv_det;

    // Third row  
    result._31 = (mat._12 * mat._23 * mat._44 - 
                 mat._13 * mat._24 * mat._31 + 
                 mat._14 * mat._21 * mat._32 - 
                 mat._14 * mat._23 * mat._31 + 
                 mat._13 * mat._24 * mat._21 - 
                 mat._11 * mat._23 * mat._34) * inv_det;
    result._32 = -(mat._11 * mat._24 * mat._34 - 
                  mat._12 * mat._23 * mat._31 + 
                  mat._13 * mat._21 * mat._32 - 
                  mat._13 * mat._22 * mat._31 + 
                  mat._12 * mat._23 * mat._21 - 
                  mat._11 * mat._23 * mat._34) * inv_det;
    result._33 = (mat._11 * mat._22 * mat._34 - 
                 mat._12 * mat._23 * mat._31 + 
                 mat._13 * mat._21 * mat._32 - 
                 mat._13 * mat._22 * mat._31 + 
                 mat._12 * mat._23 * mat._21 - 
                 mat._11 * mat._22 * mat._33) * inv_det;
    result._34 = -(mat._11 * mat._22 * mat._33 - 
                  mat._12 * mat._23 * mat._31 + 
                  mat._13 * mat._21 * mat._32 - 
                  mat._13 * mat._22 * mat._31 + 
                  mat._12 * mat._23 * mat._21 - 
                  mat._11 * mat._22 * mat._33) * inv_det;

    // Fourth row  
    result._41 = -(mat._13 * mat._22 * mat._34 - 
                  mat._12 * mat._23 * mat._31 + 
                  mat._13 * mat._21 * mat._32 - 
                  mat._13 * mat._22 * mat._31 + 
                  mat._12 * mat._23 * mat._21 - 
                  mat._11 * mat._23 * mat._34) * inv_det;
    result._42 = (mat._11 * mat._24 * mat._34 - 
                 mat._12 * mat._23 * mat._31 + 
                 mat._13 * mat._21 * mat._32 - 
                 mat._13 * mat._22 * mat._31 + 
                 mat._12 * mat._23 * mat._21 - 
                 mat._11 * mat._23 * mat._34) * inv_det;
    result._43 = -(mat._11 * mat._22 * mat._33 - 
                  mat._12 * mat._23 * mat._31 + 
                  mat._13 * mat._21 * mat._32 - 
                  mat._13 * mat._22 * mat._31 + 
                  mat._12 * mat._23 * mat._21 - 
                  mat._11 * mat._22 * mat._33) * inv_det;
    result._44 = (mat._11 * mat._22 * mat._33 - 
                 mat._12 * mat._23 * mat._31 + 
                 mat._13 * mat._21 * mat._32 - 
                 mat._13 * mat._22 * mat._31 + 
                 mat._12 * mat._23 * mat._21 - 
                 mat._11 * mat._22 * mat._33) * inv_det;

    return result;
}

mat4 Translation(float x, float y, float z)
{
    return {1.0f, 0.0f, 0.0f, 0.0f,
                0.0f, 1.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 1.0f, 0.0f,
                x, y, z, 1.0f};
}

mat4 Translation(const vec3& pos)
{
    return {1.0f, 0.0f, 0.0f, 0.0f,
                0.0f, 1.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 1.0f, 0.0f,
                pos.x, pos.y, pos.z, 1.0f};
}

vec3 GetTranslation(const mat4& mat)
{
    return vec3(mat._41, mat._42, mat._43);
}

mat4 Scale(float x, float y, float z)
{
    return {x, 0.0f, 0.0f, 0.0f,
                0.0f, y, 0.0f, 0.0f,
                0.0f, 0.0f, z, 0.0f,
                0.0f, 0.0f, 0.0f, 1.0f};
}

mat4 Scale(const vec3& vec)
{
    return {vec.x, 0.0f, 0.0f, 0.0f,
                0.0f, vec.y, 0.0f, 0.0f,
                0.0f, 0.0f, vec.z, 0.0f,
                0.0f, 0.0f, 0.0f, 1.0f};
}

vec3 GetScale(const mat4& mat)
{
    return vec3(mat._11, mat._22, mat._33);
}

mat3 Rotation3x3(float pitch, float yaw, float roll)
{
    return ZRotation3x3(roll) * XRotation3x3(pitch) * YRotation3x3(yaw);
}

mat4 Rotation(float pitch, float yaw, float roll)
{
    return ZRotation(roll) * XRotation(pitch) * YRotation(yaw);
}

mat3 ZRotation3x3(float angle)
{
    angle = DEG2RAD(angle);
    return {cosf(angle), sinf(angle), 0.0f,
                -sinf(angle), cosf(angle), 0.0f,
                0.0f, 0.0f, 1.0f};
}

mat4 ZRotation(float angle)
{
    angle = DEG2RAD(angle);
    return {cosf(angle), sinf(angle), 0.0f, 0.0f,
                -sinf(angle), cosf(angle), 0.0f, 0.0f,
                0.0f, 0.0f, 1.0f, 0.0f,
                0.0f, 0.0f, 0.0f, 1.0f};
}

mat3 XRotation3x3(float angle) {
    angle = DEG2RAD(angle);
    return mat3(
    1.0f, 0.0f, 0.0f,
    0.0f, cosf(angle), sinf(angle),
    0.0f, -sinf(angle), cos(angle)
    );
}

mat4 XRotation(float angle) {
    angle = DEG2RAD(angle);
    return mat4(
    1.0f, 0.0f, 0.0f, 0.0f,
    0.0f, cosf(angle), sinf(angle), 0.0f,
    0.0f, -sinf(angle), cos(angle), 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f
    );
}

mat3 YRotation3x3(float angle) {
    angle = DEG2RAD(angle);
    return mat3(
    cosf(angle), 0.0f, -sinf(angle),
    0.0f, 1.0f, 0.0f,
    sinf(angle), 0.0f, cosf(angle)
    );
}

mat4 YRotation(float angle) {
    angle = DEG2RAD(angle);
    return mat4(
    cosf(angle), 0.0f, -sinf(angle), 0.0f,
    0.0f, 1.0f, 0.0f, 0.0f,
    sinf(angle), 0.0f, cosf(angle), 0.0f,
    0.0f, 0.0f, 0.0f, 1.0f
    );
}

mat3 AxisAngle3x3(const vec3& axis, float angle)
{
    angle = DEG2RAD(angle);
    float c = cosf(angle);
    float s = sinf(angle);
    float t = 1.0 - cosf(angle);
    
    float x = axis.x;
    float y = axis.y;
    float z = axis.z;

    if(!CMP(MagnitudeSq(axis), 1.0f))
    {
        float inv_len = 1.0f / Magnitude(axis);
        x *= inv_len;
        y *= inv_len;
        z *= inv_len;
    }

    return {t * (x * x) + c,t * x * y + s * z,t * x * z - s * y,
                t * x * y - s * z,t * (y * y) + c,t * y * z + s * x,
                t * x * z + s * y,t * y * z - s * x,t * (z * z) + c};
}

mat4 AxisAngle(const vec3& axis, float angle)
{
    angle = DEG2RAD(angle);
    float c = cosf(angle);
    float s = sinf(angle);
    float t = 1.0 - cosf(angle);
    
    float x = axis.x;
    float y = axis.y;
    float z = axis.z;

    if(!CMP(MagnitudeSq(axis), 1.0f))
    {
        float inv_len = 1.0f / Magnitude(axis);
        x *= inv_len;
        y *= inv_len;
        z *= inv_len;
    }

    return {t*(x*x) + c, t*x*y + s*z, t*x*z - s*y, 0.0f,
                t*x*y - s*z, t*(y*y) + c, t*y*z + s*x, 0.0f,
                t*x*z + s*y, t*y*z - s*x, t*(z*z) + c, 0.0f,
                0.0f, 0.0f, 0.0f, 1.0f};
}
