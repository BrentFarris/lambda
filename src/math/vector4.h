#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

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
		return sqrt(x * x + y * y + z * z + w * w);
	}

	double magnitude() {
		return sqrt((x * x) + (y * y) + (z * z) + (w * w));
	}

	void normalize() {
		double mag = magnitude();
		x /= mag;
		y /= mag;
		z /= mag;
		w /= mag;
	}

	Vector4 normalized() {
		double mag = magnitude();
		return Vector4(x / mag, y / mag, z / mag, w / mag);
	}

	double dot(Vector4& other) {
		return x * other.x + y * other.y + z * other.z + w * other.w;
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
		return Vector4(from.x + (to.x - from.x) * t,
			from.y + (to.y - from.y) * t,
			from.z + (to.z - from.z) * t,
			from.w + (to.w - from.w) * t);
	}

	Vector4 operator+(Vector4& other) {
		return Vector4(x + other.x, y + other.y, z + other.z, w + other.w);
	}

	Vector4 operator-(Vector4& other) {
		return Vector4(x - other.x, y - other.y, z - other.z, w - other.w);
	}

	Vector4 operator*(double scalar) {
		return Vector4(x * scalar, y * scalar, z * scalar, w * scalar);
	}

	Vector4 operator*(Vector4& other) {
		return Vector4(x * other.x, y * other.y, z * other.z, w * other.w);
	}

	Vector4 operator/(double scalar) {
		return Vector4(x / scalar, y / scalar, z / scalar, w / scalar);
	}

	Vector4 operator/(Vector4 other) {
		return Vector4(x / other.x, y / other.y, z / other.z, w / other.w);
	}

	void operator+=(Vector4& other) {
		x += other.x;
		y += other.y;
		z += other.z;
		w += other.w;
	}

	void operator-=(Vector4& other) {
		x -= other.x;
		y -= other.y;
		z -= other.z;
		w -= other.w;
	}

	void operator*=(double scalar) {
		x *= scalar;
		y *= scalar;
		z *= scalar;
		w *= scalar;
	}

	void operator*=(Vector4& other) {
		x *= other.x;
		y *= other.y;
		z *= other.z;
		w *= other.w;
	}

	void operator/=(double scalar) {
		x /= scalar;
		y /= scalar;
		z /= scalar;
		w /= scalar;
	}

	void operator/=(Vector4& other) {
		x /= other.x;
		y /= other.y;
		z /= other.z;
		w /= other.w;
	}

	bool operator==(Vector4& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) < epsilon
			&& std::abs(y - other.y) < epsilon
			&& std::abs(z - other.z) < epsilon
			&& std::abs(w - other.w) < epsilon;
	}

	bool operator!=(Vector4& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) > epsilon
			|| std::abs(y - other.y) > epsilon
			|| std::abs(z - other.z) > epsilon
			|| std::abs(w - other.w) > epsilon;
	}
};
