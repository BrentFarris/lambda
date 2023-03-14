#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "hardware_vector.h"

struct Vector4 {
	union {
		double linear[4];
		struct {
			double x;
			double y;
			double z;
			double w;
		};
	};

	Vector4() {
		x = 0.0;
		y = 0.0;
		z = 0.0;
		w = 0.0;
	}

	Vector4(double x, double y, double z, double w) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}

	Vector4(double xyz) {
		this->x = xyz;
		this->y = xyz;
		this->z = xyz;
		this->w = xyz;
	}

	static Vector4 zero() {
		return Vector4(0.0);
	}

	static Vector4 one() {
		return Vector4(1.0);
	}
	
	static Vector4 identity() {
		return Vector4(0.0, 0.0, 0.0, 1.0);
	}

	double length() {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_ld(linear);
		vf64x4 v = vf64x4_mul(l, r);
		double m = vf64x4_sum(v);
#else
		double m = (x * x) + (y * y) + (z * z) + (w * w);
#endif
		return sqrt(m);
	}

	void normalize() {
		double mag = length();
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_set(mag, mag, mag, mag);
		vf64x4_st(linear, vf64x4_div(l, r));
#else
		x /= mag;
		y /= mag;
		z /= mag;
		w /= mag;
#endif
	}

	Vector4 normalized() {
		double mag = length();
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		Vector4 ret = *this;
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_set(mag, mag, mag, mag);
		vf64x4_st(ret.linear, vf64x4_div(l, r));
		return ret;
#else
		return Vector4(x / mag, y / mag, z / mag, w / mag);
#endif
	}

	double dot(Vector4& other) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_ld(other.linear);
		vf64x4 v = vf64x4_mul(l, r);
		return vf64x4_sum(v);
#else
		return x * other.x + y * other.y + z * other.z + w * other.w;
#endif
	}

	Vector4 abs() {
		return Vector4(std::abs(x), std::abs(y), std::abs(z), std::abs(w));
	}

	std::string to_string() {
		std::stringstream ss;
		ss << std::to_string(x) << ", " << std::to_string(y) << ", " << std::to_string(z) << ", " << std::to_string(w);
		return ss.str();
	}

	void from_string(std::string& str) {
		char skip;
		std::stringstream ss;
		ss << str;
		ss >> x >> skip >> y >> skip >> z >> skip >> w;
	}

	void print() {
		std::cout << "Vector4<" << x << ", " << y << ", " << z << ", " << w << ">";
	}

	static Vector4 min(Vector4& a, Vector4& b) {
		return Vector4(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z), std::min(a.w, b.w));
	}

	static Vector4 max(Vector4& a, Vector4& b) {
		return Vector4(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z), std::max(a.w, b.w));
	}

	static Vector4 lerp(Vector4& from, Vector4& to, double t) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		Vector4 ret = to;
		vf64x4 l = vf64x4_ld(ret.linear);
		vf64x4 r = vf64x4_ld(from.linear);
		vf64x4 vt = vf64x4_set(t, t, t, t);
		vf64x4_st(ret.linear, vf64x4_add(l, vf64x4_mul(vf64x4_sub(l, r), vt)));
		return ret;
#else
		return Vector4(from.x + (to.x - from.x) * t,
			from.y + (to.y - from.y) * t,
			from.z + (to.z - from.z) * t,
			from.w + (to.w - from.w) * t);
#endif
	}

	Vector4 operator+(Vector4& other) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		Vector4 ret = *this;
		vf64x4 l = vf64x4_ld(ret.linear);
		vf64x4 r = vf64x4_ld(other.linear);
		vf64x4_st(ret.linear, vf64x4_add(l, r));
		return ret;
#else
		return Vector4(x + other.x, y + other.y, z + other.z, w + other.w);
#endif
	}

	Vector4 operator-(Vector4& other) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		Vector4 ret = *this;
		vf64x4 l = vf64x4_ld(ret.linear);
		vf64x4 r = vf64x4_ld(other.linear);
		vf64x4_st(ret.linear, vf64x4_sub(l, r));
		return ret;
#else
		return Vector4(x - other.x, y - other.y, z - other.z, w - other.w);
#endif
	}

	Vector4 operator*(double scalar) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		Vector4 ret = *this;
		vf64x4 l = vf64x4_ld(ret.linear);
		vf64x4 r = vf64x4_set(scalar, scalar, scalar, scalar);
		vf64x4_st(ret.linear, vf64x4_mul(l, r));
		return ret;
#else
		return Vector4(x * scalar, y * scalar, z * scalar, w * scalar);
#endif
	}

	Vector4 operator*(Vector4& other) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		Vector4 ret = *this;
		vf64x4 l = vf64x4_ld(ret.linear);
		vf64x4 r = vf64x4_ld(other.linear);
		vf64x4_st(ret.linear, vf64x4_mul(l, r));
		return ret;
#else
		return Vector4(x * other.x, y * other.y, z * other.z, w * other.w);
#endif
	}

	Vector4 operator/(double scalar) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		Vector4 ret = *this;
		vf64x4 l = vf64x4_ld(ret.linear);
		vf64x4 r = vf64x4_set(scalar, scalar, scalar, scalar);
		vf64x4_st(ret.linear, vf64x4_div(l, r));
		return ret;
#else
		return Vector4(x / scalar, y / scalar, z / scalar, w / scalar);
#endif
	}

	Vector4 operator/(Vector4 other) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		Vector4 ret = *this;
		vf64x4 l = vf64x4_ld(ret.linear);
		vf64x4 r = vf64x4_ld(other.linear);
		vf64x4_st(ret.linear, vf64x4_div(l, r));
		return ret;
#else
		return Vector4(x / other.x, y / other.y, z / other.z, w / other.w);
#endif
	}

	void operator+=(Vector4& other) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_ld(other.linear);
		vf64x4_st(linear, vf64x4_add(l, r));
#else
		x += other.x;
		y += other.y;
		z += other.z;
		w += other.w;
#endif
	}

	void operator-=(Vector4& other) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_ld(other.linear);
		vf64x4_st(linear, vf64x4_sub(l, r));
#else
		x -= other.x;
		y -= other.y;
		z -= other.z;
		w -= other.w;
#endif
	}

	void operator*=(double scalar) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_set(scalar, scalar, scalar, scalar);
		vf64x4_st(linear, vf64x4_mul(l, r));
#else
		x *= scalar;
		y *= scalar;
		z *= scalar;
		w *= scalar;
#endif
	}

	void operator*=(Vector4& other) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_ld(other.linear);
		vf64x4_st(linear, vf64x4_add(l, r));
#else
		x *= other.x;
		y *= other.y;
		z *= other.z;
		w *= other.w;
#endif
	}

	void operator/=(double scalar) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_set(scalar, scalar, scalar, scalar);
		vf64x4_st(linear, vf64x4_div(l, r));
#else
		x /= scalar;
		y /= scalar;
		z /= scalar;
		w /= scalar;
#endif
	}

	void operator/=(Vector4& other) {
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_ld(other.linear);
		vf64x4_st(linear, vf64x4_div(l, r));
#else
		x /= other.x;
		y /= other.y;
		z /= other.z;
		w /= other.w;
#endif
	}

	bool operator==(Vector4& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		double d[4];
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_ld(other.linear);
		vf64x4_st(d, vf64x4_sub(l, r));
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

	bool operator!=(Vector4& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		double d[4];
#if defined(USE_SIMD) && !defined(SIMD_NEON)
		vf64x4 l = vf64x4_ld(linear);
		vf64x4 r = vf64x4_ld(other.linear);
		vf64x4_st(d, vf64x4_sub(l, r));
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
