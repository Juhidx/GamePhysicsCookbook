﻿#pragma once

#ifndef _H_MATH_VECTORS_
#define _H_MATH_VECTORS_
#define RAD2DEG(x) ((x) * 57.295754f)
#define DEG2RAD(x) ((x) * 0.0174533f)
#include <assert.h>

typedef struct vec2
{
    union
    {
        struct {
            float x;
            float y;
        };

        float asArray[2];
    };

    const float& operator[] (int i)
    {
        assert(i >= 0 && i < 2);
        return asArray[i];
    }

    vec2() : x(0.0f), y(0.0f) {}
    vec2(float _x, float _y) : x(_x), y(_x) {};
} vec2;

typedef struct vec3
{
    union 
    {
        struct
        {
            float x;
            float y;
            float z;
        };

        float asArray[3];
    };

    const float& operator[] (int i)
    {
        assert(i >= 0 && i < 3);
        return asArray[i];
    }
    
    vec3() : x(0.0f), y(0.0f), z(0.0f) { }
    vec3(float _x, float _y, float _z) : x(_x), y(_y), z(_z) { }
} vec3;

vec2 operator+ (const vec2& l, const vec2& r);
vec2 operator- (const vec2& l, const vec2& r);
vec2 operator* (const vec2& l, const vec2& r);
vec2 operator* (const vec2& l, const float s);
bool operator== (const vec2& l, const vec2& r);
bool operator!= (const vec2& l, const vec2& r);

vec3 operator+ (const vec3& l, const vec3& r);
vec3 operator- (const vec3& l, const vec3& r);
vec3 operator* (const vec3& l, const vec3& r);
vec3 operator* (const vec3& l, const float s);
bool operator== (const vec3& l, const vec3& r);
bool operator!= (const vec3& l, const vec3& r);

float Dot(const vec2& l, const vec2& r);
float Dot(const vec3& l, const vec3& r);
float Magnitude(const vec2& v);
float Magnitude(const vec3& v);
float MagnitudeSq(const vec2& v);
float MagnitudeSq(const vec3& v);
float Distance(const vec2& p1, const vec2& p2);
float Distance(const vec3& p1, const vec3& p2);
float DistanceSq(const vec2& p1, const vec2& p2);
float DistanceSq(const vec3& p1, const vec3& p2);
void Normalize(vec2& v);
void Normalize(vec3& v);
vec2 Normalized(const vec2& v);
vec3 Normalized(const vec3& v);
vec3 Cross(const vec3& l, const vec3& r);
float Angle(const vec2& l, const vec2& r);
float Angle(const vec3& l, const vec3& r);
vec2 Project(const vec2& length, const vec2& direction);
vec3 Project(const vec3& length, const vec3& direction);
vec2 Perpendicular(const vec2& len, const vec2& dir);
vec3 Perpendicular(const vec3& len, const vec3& dir);
vec2 Reflection(const vec2& vec, const vec2& normal);
vec3 Reflection(const vec3& vec, const vec3& normal);

#endif

class vectors
{
public:
};