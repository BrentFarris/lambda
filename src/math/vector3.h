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

	Vector3() {
		x = 0.0;
		y = 0.0;
		z = 0.0;
	}

	Vector3(double x, double y, double z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	Vector3(double xyz) {
		this->x = xyz;
		this->y = xyz;
		this->z = xyz;
	}

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

	double length() {
		return sqrt(x * x + y * y + z * z);
	}

	double magnitude() {
		return sqrt((x * x) + (y * y) + (z * z));
	}

	void normalize() {
		double mag = magnitude();
		x /= mag;
		y /= mag;
		z /= mag;
	}

	Vector3 normalized() {
		double mag = magnitude();
		return Vector3(x / mag, y / mag, z / mag);
	}

	double dot(Vector3& other) {
		return x * other.x + y * other.y + z * other.z;
	}

	Vector3 cross(Vector3& other) {
		return Vector3(y * other.z - z * other.y,
			z * other.x - x * other.z,
			x * other.y - y * other.x);
	}

	double distance(Vector3& target) {
		return sqrt((target.x - x) * (target.x - x)
			+ (target.y - y) * (target.y - y)
			+ (target.z - z) * (target.z - z));
	}

	Vector3 abs() {
		return Vector3(std::abs(x), std::abs(y), std::abs(z));
	}

	std::string to_string() {
		std::stringstream ss;
		ss << std::to_string(x) << ", " << std::to_string(y) << ", " << std::to_string(z);
		return ss.str();
	}

	void from_string(std::string& str) {
		char skip;
		std::stringstream ss;
		ss << str;
		ss >> x >> skip >> y >> skip >> z;
	}

	void print() {
		std::cout << "Vector3<" << x << ", " << y << ", " << z << ">";
	}

	static Vector3 min(Vector3& a, Vector3& b) {
		return Vector3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
	}

	static Vector3 max(Vector3& a, Vector3& b) {
		return Vector3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
	}

	static Vector3 max_abs(Vector3& a, Vector3& b) {
		return Vector3(std::max(std::abs(a.x), std::abs(b.x)),
			std::max(std::abs(a.y), std::abs(b.y)),
			std::max(std::abs(a.z), std::abs(b.z)));
	}

	static Vector3 lerp(Vector3& from, Vector3& to, double t) {
		return Vector3(from.x + (to.x - from.x) * t,
			from.y + (to.y - from.y) * t,
			from.z + (to.z - from.z) * t);
	}

	Vector3 operator+(Vector3& other) {
		return Vector3(x + other.x, y + other.y, z + other.z);
	}

	Vector3 operator-(Vector3& other) {
		return Vector3(x - other.x, y - other.y, z - other.z);
	}

	Vector3 operator*(double scalar) {
		return Vector3(x * scalar, y * scalar, z * scalar);
	}

	Vector3 operator*(Vector3& other) {
		return Vector3(x * other.x, y * other.y, z * other.z);
	}

	Vector3 operator/(double scalar) {
		return Vector3(x / scalar, y / scalar, z / scalar);
	}

	Vector3 operator/(Vector3 other) {
		return Vector3(x / other.x, y / other.y, z / other.z);
	}

	void operator+=(Vector3& other) {
		x += other.x;
		y += other.y;
		z += other.z;
	}

	void operator-=(Vector3& other) {
		x -= other.x;
		y -= other.y;
		z -= other.z;
	}

	void operator*=(double scalar) {
		x *= scalar;
		y *= scalar;
		z *= scalar;
	}

	void operator*=(Vector3& other) {
		x *= other.x;
		y *= other.y;
		z *= other.z;
	}

	void operator/=(double scalar) {
		x /= scalar;
		y /= scalar;
		z /= scalar;
	}

	void operator/=(Vector3& other) {
		x /= other.x;
		y /= other.y;
		z /= other.z;
	}

	bool operator==(Vector3& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) < epsilon
			&& std::abs(y - other.y) < epsilon
			&& std::abs(z - other.z) < epsilon;
	}

	bool operator!=(Vector3& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) > epsilon
			|| std::abs(y - other.y) > epsilon
			|| std::abs(z - other.z) > epsilon;
	}
};
