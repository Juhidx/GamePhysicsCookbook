﻿#include "vectors.h"
#include <cmath>
#include <cfloat>

#define CMP(x, y)                       \
    (fabsf((x)+(y)) <= FLT_EPSILON *    \
        fmaxf(1.0f,                     \
        fmaxf(fabsf(x), fabsf(y)))      \
    )

vec2 operator+(const vec2& l, const vec2& r)
{
    return vec2{l.x + r.x, l.y + r.y};
}

vec2 operator-(const vec2& l, const vec2& r)
{
    return vec2{l.x - r.x, l.y - r.y};
}

vec2 operator*(const vec2& l, const vec2& r)
{
    return vec2{l.x * r.x, l.y * r.y};
}

vec2 operator*(const vec2& l, const float s)
{
    return vec2{l.x * s, l.y * s};
}

bool operator==(const vec2& l, const vec2& r)
{
    return CMP(l.x, r.x) && CMP(l.y, r.y);
}

bool operator!=(const vec2& l, const vec2& r)
{
    return !(operator==(l, r));
}

vec3 operator+(const vec3& l, const vec3& r)
{
    return vec3{l.x + r.x, l.y + r.y, l.z + r.z};
}

vec3 operator-(const vec3& l, const vec3& r)
{
    return vec3{l.x - r.x, l.y - r.y, l.z - r.z};
}

vec3 operator*(const vec3& l, const vec3& r)
{
    return vec3{l.x * r.x, l.y * r.y, l.z * r.z};
}

vec3 operator*(const vec3& l, const float s)
{
    return vec3{l.x * s, l.y * s, l.z * s};
}

bool operator==(const vec3& l, const vec3& r)
{
    return CMP(l.x, r.x) && CMP(l.y, r.y) && CMP(l.z, r.z);
}

bool operator!=(const vec3& l, const vec3& r)
{
    return !(operator==(l, r));
}

float Dot(const vec2& l, const vec2& r)
{
    return l.x * r.x + l.y * r.y;
}

float Dot(const vec3& l, const vec3& r)
{
    return l.x * r.x + l.y * r.y + l.z * r.z;
}

float Magnitude(const vec2& v)
{
    return sqrtf(Dot(v, v));
}

float Magnitude(const vec3& v)
{
    return sqrtf(Dot(v, v));
}

float MagnitudeSq(const vec2& v)
{
    return Dot(v, v);
}

float MagnitudeSq(const vec3& v)
{
    return Dot(v, v);
}

float Distance(const vec2& p1, const vec2& p2)
{
    vec2 t = operator-(p1, p2);
    return Magnitude(t);
}

float Distance(const vec3& p1, const vec3& p2)
{
    vec3 t = operator-(p1, p2);
    return Magnitude(t);
}

float DistanceSq(const vec2& p1, const vec2& p2)
{
    vec2 t = operator-(p1, p2);
    return MagnitudeSq(t);
}

float DistanceSq(const vec3& p1, const vec3& p2)
{
    vec3 t = operator-(p1, p2);
    return MagnitudeSq(t);
}

void Normalize(vec2& v)
{
    v = operator*(v, (1.0f / Magnitude(v)));
}

void Normalize(vec3& v)
{
    v = operator*(v, (1.0f / Magnitude(v)));
}

vec2 Normalized(const vec2& v)
{
    return operator*(v, (1.0f / Magnitude(v)));
}

vec3 Normalized(const vec3& v)
{
    return operator*(v, (1.0f / Magnitude(v)));
}

vec3 Cross(const vec3& l, const vec3& r)
{
    vec3 result;
    result.x = l.y * r.z - l.z * r.y;
    result.y = l.z * r.x - l.x * r.z;
    result.z = l.x * r.y - l.y * r.x;
    return result;
}

float Angle(const vec2& l, const vec2& r)
{
    float m = sqrtf(MagnitudeSq(l) * MagnitudeSq(r));
    return acos(Dot(l, r) / m);
}

float Angle(const vec3& l, const vec3& r)
{
    float m = sqrtf(MagnitudeSq(l) * MagnitudeSq(r));
    return acos(Dot(l, r) / m);
}

vec2 Project(const vec2& length, const vec2& direction)
{
    float dot = Dot(length, direction);
    float maqSq = MagnitudeSq(direction);

    if(maqSq == 0)
    {
        return {0,0};
    }
    
    return direction * (dot/maqSq);
}

vec3 Project(const vec3& length, const vec3& direction)
{
    float dot = Dot(length, direction);
    float maqSq = MagnitudeSq(direction);

    if(maqSq == 0)
    {
        return {0,0,0};
    }
    
    return direction * (dot/maqSq);
}

vec2 Perpendicular(const vec2& len, const vec2& dir)
{
    return len - Project(len, dir);
}

vec3 Perpendicular(const vec3& len, const vec3& dir)
{
    return len - Project(len, dir);
}

vec2 Reflection(const vec2& vec, const vec2& normal)
{
    float d = Dot(vec, normal);
    return vec - normal * (d * 2.0f);
}

vec3 Reflection(const vec3& vec, const vec3& normal)
{
    float d = Dot(vec, normal);
    return vec - normal * (d * 2.0f);
}
