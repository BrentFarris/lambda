#ifndef LAMBDA_MATH_VECTOR4_H
#define LAMBDA_MATH_VECTOR4_H

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "hardware_vector.h"

struct Vector4 {
	union {
		float linear[4];
		struct {
			float x;
			float y;
			float z;
			float w;
		};
	};

	Vector4() : x(0.0F), y(0.0F), z(0.0F), w(0.0F) {}
	Vector4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
	Vector4(float xyzw) : x(xyzw), y(xyzw), z(xyzw), w(xyzw) {}
	Vector4(float* xyzw) : x(xyzw[0]), y(xyzw[1]), z(xyzw[2]), w(xyzw[3]) {}

	static Vector4 zero() {
		return Vector4(0.0F);
	}

	static Vector4 one() {
		return Vector4(1.0F);
	}
	
	static Vector4 identity() {
		return Vector4(0.0F, 0.0F, 0.0F, 1.0F);
	}

	static float length_linear(float vec[4]) {
#if defined(USE_SIMD)
		vf32x4 v = vf32x4_ld((float*)vec);
		vf32x4 res = vf32x4_mul(v, v);
		float m = vf32x4_sum(res);
#else
		float m = (x * x) + (y * y) + (z * z) + (w * w);
#endif
		return sqrt(m);
	}

	inline float length() {
		return length_linear(linear);
	}

	static void normalize_linear(float vec[4]) {
		float mag = length_linear(vec);
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(vec);
		vf32x4 r = vf32x4_set(mag, mag, mag, mag);
		vf32x4_st(vec, vf32x4_div(l, r));
#else
		vec[0] /= mag;
		vec[1] /= mag;
		vec[2] /= mag;
		vec[3] /= mag;
#endif
	}

	inline void normalize() {
		normalize_linear(linear);
	}

	Vector4 normalized() {
		float mag = length();
#if defined(USE_SIMD)
		Vector4 ret = *this;
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_set(mag, mag, mag, mag);
		vf32x4_st(ret.linear, vf32x4_div(l, r));
		return ret;
#else
		return Vector4(x / mag, y / mag, z / mag, w / mag);
#endif
	}

	Vector4 abs() {
		return Vector4(std::abs(x), std::abs(y), std::abs(z), std::abs(w));
	}

	std::string to_string() const {
		std::stringstream ss;
		ss << std::to_string(x) << ", " << std::to_string(y) << ", " << std::to_string(z) << ", " << std::to_string(w);
		return ss.str();
	}

	void from_string(const std::string& str) {
		char skip;
		std::stringstream ss;
		ss << str;
		ss >> x >> skip >> y >> skip >> z >> skip >> w;
	}

	void print() {
		std::cout << "Vector4<" << x << ", " << y << ", " << z << ", " << w << ">";
	}

	static float dot_linear(const float lhs[4], const float rhs[4]) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(lhs);
		vf32x4 r = vf32x4_ld(rhs);
		vf32x4 v = vf32x4_mul(l, r);
		return vf32x4_sum(v);
#else
		return x * other.x + y * other.y + z * other.z + w * other.w;
#endif
	}

	static inline float dot(const Vector4& lhs, const Vector4& rhs) {
		return dot_linear(lhs.linear, rhs.linear);
	}

	static Vector4 min(const Vector4& a, const Vector4& b) {
		return Vector4(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z), std::min(a.w, b.w));
	}

	static Vector4 max(const Vector4& a, const Vector4& b) {
		return Vector4(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z), std::max(a.w, b.w));
	}

	static Vector4 lerp(const Vector4& from, const Vector4& to, float t) {
#if defined(USE_SIMD)
		Vector4 ret = to;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_ld(from.linear);
		vf32x4 vt = vf32x4_set(t, t, t, t);
		vf32x4_st(ret.linear, vf32x4_add(l, vf32x4_mul(vf32x4_sub(l, r), vt)));
		return ret;
#else
		return Vector4(from.x + (to.x - from.x) * t,
			from.y + (to.y - from.y) * t,
			from.z + (to.z - from.z) * t,
			from.w + (to.w - from.w) * t);
#endif
	}

	Vector4 operator+(const Vector4& other) const {
#if defined(USE_SIMD)
		Vector4 ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(ret.linear, vf32x4_add(l, r));
		return ret;
#else
		return Vector4(x + other.x, y + other.y, z + other.z, w + other.w);
#endif
	}

	Vector4 operator-(const Vector4& other) const {
#if defined(USE_SIMD)
		Vector4 ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(ret.linear, vf32x4_sub(l, r));
		return ret;
#else
		return Vector4(x - other.x, y - other.y, z - other.z, w - other.w);
#endif
	}

	Vector4 operator*(float scalar) const {
#if defined(USE_SIMD)
		Vector4 ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
		vf32x4_st(ret.linear, vf32x4_mul(l, r));
		return ret;
#else
		return Vector4(x * scalar, y * scalar, z * scalar, w * scalar);
#endif
	}

	Vector4 operator*(const Vector4& other) const {
#if defined(USE_SIMD)
		Vector4 ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(ret.linear, vf32x4_mul(l, r));
		return ret;
#else
		return Vector4(x * other.x, y * other.y, z * other.z, w * other.w);
#endif
	}

	Vector4 operator/(float scalar) const {
#if defined(USE_SIMD)
		Vector4 ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
		vf32x4_st(ret.linear, vf32x4_div(l, r));
		return ret;
#else
		return Vector4(x / scalar, y / scalar, z / scalar, w / scalar);
#endif
	}

	Vector4 operator/(const Vector4& other) const {
#if defined(USE_SIMD)
		Vector4 ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(ret.linear, vf32x4_div(l, r));
		return ret;
#else
		return Vector4(x / other.x, y / other.y, z / other.z, w / other.w);
#endif
	}

	void operator+=(const Vector4& other) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(linear, vf32x4_add(l, r));
#else
		x += other.x;
		y += other.y;
		z += other.z;
		w += other.w;
#endif
	}

	void operator-=(const Vector4& other) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(linear, vf32x4_sub(l, r));
#else
		x -= other.x;
		y -= other.y;
		z -= other.z;
		w -= other.w;
#endif
	}

	void operator*=(float scalar) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
		vf32x4_st(linear, vf32x4_mul(l, r));
#else
		x *= scalar;
		y *= scalar;
		z *= scalar;
		w *= scalar;
#endif
	}

	void operator*=(const Vector4& other) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(linear, vf32x4_add(l, r));
#else
		x *= other.x;
		y *= other.y;
		z *= other.z;
		w *= other.w;
#endif
	}

	void operator/=(float scalar) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
		vf32x4_st(linear, vf32x4_div(l, r));
#else
		x /= scalar;
		y /= scalar;
		z /= scalar;
		w /= scalar;
#endif
	}

	void operator/=(const Vector4& other) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(linear, vf32x4_div(l, r));
#else
		x /= other.x;
		y /= other.y;
		z /= other.z;
		w /= other.w;
#endif
	}

	bool operator==(const Vector4& other) const {
		float epsilon = std::numeric_limits<float>::epsilon();
		float d[4];
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(d, vf32x4_sub(l, r));
#else
		d[0] = x - other.x;
		d[1] = y - other.y;
		d[2] = z - other.z;
		d[3] = w - other.w;
#endif
		// TODO:  Could probably use the abs of the sum rather than each
		return std::abs(d[0]) < epsilon
			&& std::abs(d[1]) < epsilon
			&& std::abs(d[2]) < epsilon
			&& std::abs(d[3]) < epsilon;
	}

	bool operator!=(const Vector4& other) const {
		float epsilon = std::numeric_limits<float>::epsilon();
		float d[4];
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(d, vf32x4_sub(l, r));
#else
		d[0] = x - other.x;
		d[1] = y - other.y;
		d[2] = z - other.z;
		d[3] = w - other.w;
#endif
		// TODO:  Could probably use the abs of the sum rather than each
		return std::abs(d[0]) > epsilon
			&& std::abs(d[1]) > epsilon
			&& std::abs(d[2]) > epsilon
			&& std::abs(d[3]) > epsilon;
	}
};

#endif
