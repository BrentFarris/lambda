#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

struct Vector2 {
	union {
		double linear[2];
		struct {
			double x;
			double y;
		};
	};

	Vector2() {
		x = 0.0;
		y = 0.0;
	}

	Vector2(double x, double y) {
		this->x = x;
		this->y = y;
	}

	Vector2(double xy) {
		this->x = xy;
		this->y = xy;
	}

	static Vector2 zero() {
		return Vector2(0.0, 0.0);
	}

	static Vector2 one() {
		return Vector2(1.0, 1.0);
	}

	static Vector2 up() {
		return Vector2(0.0, 1.0);
	}

	static Vector2 down() {
		return Vector2(0.0, -1.0);
	}

	static Vector2 left() {
		return Vector2(-1.0, 0.0);
	}

	static Vector2 right() {
		return Vector2(1.0, 0.0);
	}

	double length() {
		return sqrt(x * x + y * y);
	}

	double magnitude() {
		return sqrt((x * x) + (y * y));
	}

	void normalize() {
		double mag = magnitude();
		x /= mag;
		y /= mag;
	}

	Vector2 normalized() {
		double mag = magnitude();
		return Vector2(x / mag, y / mag);
	}

	double dot(Vector2& other) {
		return x * other.x + y * other.y;
	}

	double distance(Vector2& target) {
		return sqrt(((target.x - x) * (target.x - x)) + ((target.y - y) * (target.y - y)));
	}

	double angle(Vector2& target) {
		return (atan2(x - target.x, y - target.y) * (180.0F / 3.14159265358979323846)) + 180.0;
	}

	Vector2 min(Vector2& a, Vector2& b) {
		return Vector2(std::min(a.x, b.x), std::min(a.y, b.y));
	}

	Vector2 max(Vector2& a, Vector2& b) {
		return Vector2(std::max(a.x, b.x), std::max(a.y, b.y));
	}

	static Vector2 lerp(Vector2& from, Vector2& to, double t) {
		return Vector2(from.x + (to.x - from.x) * t,
			from.y + (to.y - from.y) * t);
	}

	std::string to_string() {
		std::stringstream ss;
		ss << std::to_string(x) << ", " << std::to_string(y);
		return ss.str();
	}

	void from_string(std::string& str) {
		char skip;
		std::stringstream ss;
		ss << str;
		ss >> x >> skip >> y;
	}

	void print() {
		std::cout << "Vector2<" << x << ", " << y << ">";
	}

	Vector2 operator+(Vector2& other) {
		return Vector2(x + other.x, y + other.y);
	}

	Vector2 operator-(Vector2& other) {
		return Vector2(x - other.x, y - other.y);
	}

	Vector2 operator*(double scalar) {
		return Vector2(x * scalar, y * scalar);
	}

	Vector2 operator*(Vector2& other) {
		return Vector2(x * other.x, y * other.y);
	}

	Vector2 operator/(double scalar) {
		return Vector2(x / scalar, y / scalar);
	}

	Vector2 operator/(Vector2 other) {
		return Vector2(x / other.x, y / other.y);
	}

	void operator+=(Vector2& other) {
		x += other.x;
		y += other.y;
	}

	void operator-=(Vector2& other) {
		x -= other.x;
		y -= other.y;
	}

	void operator*=(double scalar) {
		x *= scalar;
		y *= scalar;
	}

	void operator*=(Vector2& other) {
		x *= other.x;
		y *= other.y;
	}

	void operator/=(double scalar) {
		x /= scalar;
		y /= scalar;
	}

	void operator/=(Vector2& other) {
		x /= other.x;
		y /= other.y;
	}

	bool operator==(Vector2& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) < epsilon && std::abs(y - other.y) < epsilon;
	}

	bool operator!=(Vector2& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) > epsilon || std::abs(y - other.y) > epsilon;
	}
};

