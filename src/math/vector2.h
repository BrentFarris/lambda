#ifndef LAMBDA_MATH_VECTOR2_H
#define LAMBDA_MATH_VECTOR2_H

#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

struct Vector2 {
	union {
		float linear[2];
		struct {
			float x;
			float y;
		};
	};

	Vector2() : x(0.0F), y(0.0F) {}
	Vector2(float x, float y) : x(x), y(y) {}
	Vector2(float xy) : x(xy), y(xy) {}
	Vector2(float* xy) : x(xy[0]), y(xy[1]) {}

	static Vector2 zero() {
		return Vector2(0.0F, 0.0F);
	}

	static Vector2 one() {
		return Vector2(1.0F, 1.0F);
	}

	static Vector2 up() {
		return Vector2(0.0F, 1.0F);
	}

	static Vector2 down() {
		return Vector2(0.0F, -1.0F);
	}

	static Vector2 left() {
		return Vector2(-1.0F, 0.0F);
	}

	static Vector2 right() {
		return Vector2(1.0F, 0.0F);
	}

	float length() {
		return sqrt(x * x + y * y);
	}

	float magnitude() {
		return sqrt((x * x) + (y * y));
	}

	void normalize() {
		float mag = magnitude();
		x /= mag;
		y /= mag;
	}

	Vector2 normalized() {
		float mag = magnitude();
		return Vector2(x / mag, y / mag);
	}

	float dot(const Vector2& other) {
		return x * other.x + y * other.y;
	}

	float distance(const Vector2& target) {
		return sqrt(((target.x - x) * (target.x - x)) + ((target.y - y) * (target.y - y)));
	}

	float angle(const Vector2& target) {
		return (atan2(x - target.x, y - target.y) * (180.0F / 3.14159265358979323846F)) + 180.0F;
	}

	Vector2 min(const Vector2& a, const Vector2& b) {
		return Vector2(std::min(a.x, b.x), std::min(a.y, b.y));
	}

	Vector2 max(const Vector2& a, const Vector2& b) {
		return Vector2(std::max(a.x, b.x), std::max(a.y, b.y));
	}

	static Vector2 lerp(const Vector2& from, const Vector2& to, float t) {
		return Vector2(from.x + (to.x - from.x) * t,
			from.y + (to.y - from.y) * t);
	}

	std::string to_string() const {
		std::stringstream ss;
		ss << std::to_string(x) << ", " << std::to_string(y);
		return ss.str();
	}

	void from_string(const std::string& str) {
		char skip;
		std::stringstream ss;
		ss << str;
		ss >> x >> skip >> y;
	}

	void print() {
		std::cout << "Vector2<" << x << ", " << y << ">";
	}

	Vector2 operator+(const Vector2& other) const {
		return Vector2(x + other.x, y + other.y);
	}

	Vector2 operator-(const Vector2& other) const {
		return Vector2(x - other.x, y - other.y);
	}

	Vector2 operator*(float scalar) const {
		return Vector2(x * scalar, y * scalar);
	}

	Vector2 operator*(const Vector2& other) const {
		return Vector2(x * other.x, y * other.y);
	}

	Vector2 operator/(float scalar) const {
		return Vector2(x / scalar, y / scalar);
	}

	Vector2 operator/(const Vector2& other) const {
		return Vector2(x / other.x, y / other.y);
	}

	void operator+=(const Vector2& other) {
		x += other.x;
		y += other.y;
	}

	void operator-=(const Vector2& other) {
		x -= other.x;
		y -= other.y;
	}

	void operator*=(float scalar) {
		x *= scalar;
		y *= scalar;
	}

	void operator*=(const Vector2& other) {
		x *= other.x;
		y *= other.y;
	}

	void operator/=(float scalar) {
		x /= scalar;
		y /= scalar;
	}

	void operator/=(const Vector2& other) {
		x /= other.x;
		y /= other.y;
	}

	bool operator==(const Vector2& other) const {
		float epsilon = std::numeric_limits<float>::epsilon();
		return std::abs(x - other.x) < epsilon && std::abs(y - other.y) < epsilon;
	}

	bool operator!=(const Vector2& other) const {
		float epsilon = std::numeric_limits<float>::epsilon();
		return std::abs(x - other.x) > epsilon || std::abs(y - other.y) > epsilon;
	}
};

#endif
