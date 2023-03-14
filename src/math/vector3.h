#ifndef LAMBDA_MATH_VECTOR3_H
#define LAMBDA_MATH_VECTOR3_H

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

struct Vector3 {
	union {
		double linear[3];
		struct {
			double x;
			double y;
			double z;
		};
	};

	Vector3() : x(0.0), y(0.0), z(0.0) {}
	Vector3(double x, double y, double z) : x(x), y(y), z(z) {}
	Vector3(double xyz) : x(xyz), y(xyz), z(xyz) {}
	Vector3(double* xyz) : x(xyz[0]), y(xyz[1]), z(xyz[2]) {}

	static Vector3 zero() {
		return Vector3(0.0);
	}

	static Vector3 one() {
		return Vector3(1.0);
	}

	static Vector3 up() {
		return Vector3(0.0, 1.0, 0.0);
	}

	static Vector3 down() {
		return Vector3(0.0, -1.0, 0.0);
	}

	static Vector3 left() {
		return Vector3(-1.0, 0.0, 0.0);
	}

	static Vector3 right() {
		return Vector3(1.0, 0.0, 0.0);
	}

	static Vector3 forward() {
		return Vector3(0.0, 0.0, -1.0);
	}

	static Vector3 backward() {
		return Vector3(0.0, 0.0, 1.0);
	}

	double length() const {
		return sqrt(x * x + y * y + z * z);
	}

	double magnitude() const {
		return sqrt((x * x) + (y * y) + (z * z));
	}

	void normalize() {
		double mag = magnitude();
		x /= mag;
		y /= mag;
		z /= mag;
	}

	Vector3 normal() const {
		double mag = magnitude();
		return Vector3(x / mag, y / mag, z / mag);
	}

	Vector3 abs() {
		return Vector3(std::abs(x), std::abs(y), std::abs(z));
	}

	Vector3 orthogonal() const {
		double tx = std::abs(x);
		double ty = std::abs(y);
		double tz = std::abs(z);
		Vector3 other = tx < ty
			? (tx < tz ? right() : forward())
			: (ty < tz ? up() : forward());
		return cross(*this, other);
	}

	Vector3 negative() const {
		return { -x, -y, -z };
	}

	std::string to_string() {
		std::stringstream ss;
		ss << std::to_string(x) << ", " << std::to_string(y) << ", " << std::to_string(z);
		return ss.str();
	}

	void from_string(const std::string& str) {
		char skip;
		std::stringstream ss;
		ss << str;
		ss >> x >> skip >> y >> skip >> z;
	}

	void print() {
		std::cout << "Vector3<" << x << ", " << y << ", " << z << ">";
	}

	static double distance(const Vector3& from, const Vector3& to) {
		return sqrt((to.x - from.x) * (to.x - from.x)
			+ (to.y - from.y) * (to.y - from.y)
			+ (to.z - from.z) * (to.z - from.z));
	}

	static double dot(const Vector3& lhs, const Vector3& rhs) {
		return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
	}

	static Vector3 cross(const Vector3& lhs, const Vector3& rhs) {
		return Vector3(lhs.y * rhs.z - lhs.z * rhs.y,
			lhs.z * rhs.x - lhs.x * rhs.z,
			lhs.x * rhs.y - lhs.y * rhs.x);
	}

	static Vector3 min(const Vector3& a, const Vector3& b) {
		return Vector3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
	}

	static Vector3 max(const Vector3& a, const Vector3& b) {
		return Vector3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
	}

	static Vector3 max_abs(const Vector3& a, const Vector3& b) {
		return Vector3(std::max(std::abs(a.x), std::abs(b.x)),
			std::max(std::abs(a.y), std::abs(b.y)),
			std::max(std::abs(a.z), std::abs(b.z)));
	}

	static Vector3 lerp(const Vector3& from, const Vector3& to, double t) {
		return Vector3(from.x + (to.x - from.x) * t,
			from.y + (to.y - from.y) * t,
			from.z + (to.z - from.z) * t);
	}

	Vector3 operator+(const Vector3& other) const {
		return Vector3(x + other.x, y + other.y, z + other.z);
	}

	Vector3 operator-(const Vector3& other) const {
		return Vector3(x - other.x, y - other.y, z - other.z);
	}

	Vector3 operator*(double scalar) const {
		return Vector3(x * scalar, y * scalar, z * scalar);
	}

	Vector3 operator*(const Vector3& other) const {
		return Vector3(x * other.x, y * other.y, z * other.z);
	}

	Vector3 operator/(double scalar) const {
		return Vector3(x / scalar, y / scalar, z / scalar);
	}

	Vector3 operator/(const Vector3& other) const {
		return Vector3(x / other.x, y / other.y, z / other.z);
	}

	void operator+=(const Vector3& other) {
		x += other.x;
		y += other.y;
		z += other.z;
	}

	void operator-=(const Vector3& other) {
		x -= other.x;
		y -= other.y;
		z -= other.z;
	}

	void operator*=(double scalar) {
		x *= scalar;
		y *= scalar;
		z *= scalar;
	}

	void operator*=(const Vector3& other) {
		x *= other.x;
		y *= other.y;
		z *= other.z;
	}

	void operator/=(double scalar) {
		x /= scalar;
		y /= scalar;
		z /= scalar;
	}

	void operator/=(const Vector3& other) {
		x /= other.x;
		y /= other.y;
		z /= other.z;
	}

	bool operator==(const Vector3& other) const {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) < epsilon
			&& std::abs(y - other.y) < epsilon
			&& std::abs(z - other.z) < epsilon;
	}

	bool operator!=(const Vector3& other) const {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) > epsilon
			|| std::abs(y - other.y) > epsilon
			|| std::abs(z - other.z) > epsilon;
	}
};

#endif
