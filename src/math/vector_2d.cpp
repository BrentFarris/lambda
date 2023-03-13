export module Vector2D;

import <cmath>;
import <string>;
import <sstream>;
import <iostream>;
import <algorithm>;

export struct Vector2D {
	double x;
	double y;

	Vector2D() {
		x = 0.0;
		y = 0.0;
	}

	Vector2D(double x, double y) {
		this->x = x;
		this->y = y;
	}

	Vector2D(double xy) {
		this->x = xy;
		this->y = xy;
	}

	static Vector2D zero() {
		return Vector2D(0.0, 0.0);
	}

	static Vector2D one() {
		return Vector2D(1.0, 1.0);
	}

	static Vector2D up() {
		return Vector2D(0.0, 1.0);
	}

	static Vector2D down() {
		return Vector2D(0.0, -1.0);
	}

	static Vector2D left() {
		return Vector2D(-1.0, 0.0);
	}

	static Vector2D right() {
		return Vector2D(1.0, 0.0);
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

	Vector2D normalized() {
		double mag = magnitude();
		return Vector2D(x / mag, y / mag);
	}

	double dot(Vector2D& other) {
		return x * other.x + y * other.y;
	}

	double distance(Vector2D& target) {
		return sqrt(((target.x - x) * (target.x - x)) + ((target.y - y) * (target.y - y)));
	}

	double angle(Vector2D& target) {
		return (atan2(x - target.x, y - target.y) * (180.0F / 3.14159265358979323846)) + 180.0;
	}

	Vector2D min(Vector2D& a, Vector2D& b) {
		return Vector2D(std::min(a.x, b.x), std::min(a.y, b.y));
	}

	Vector2D max(Vector2D& a, Vector2D& b) {
		return Vector2D(std::max(a.x, b.x), std::max(a.y, b.y));
	}

	Vector2D lerp(Vector2D& other, double t) {
		return Vector2D(x + (other.x - x) * t, y + (other.y - y) * t);
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
		std::cout << "Vector2D<" << x << ", " << y << ">";
	}

	Vector2D operator+(Vector2D& other) {
		return Vector2D(x + other.x, y + other.y);
	}

	Vector2D operator-(Vector2D& other) {
		return Vector2D(x - other.x, y - other.y);
	}

	Vector2D operator*(double scalar) {
		return Vector2D(x * scalar, y * scalar);
	}

	Vector2D operator*(Vector2D& other) {
		return Vector2D(x * other.x, y * other.y);
	}

	Vector2D operator/(double scalar) {
		return Vector2D(x / scalar, y / scalar);
	}

	Vector2D operator/(Vector2D other) {
		return Vector2D(x / other.x, y / other.y);
	}

	void operator+=(Vector2D& other) {
		x += other.x;
		y += other.y;
	}

	void operator-=(Vector2D& other) {
		x -= other.x;
		y -= other.y;
	}

	void operator*=(double scalar) {
		x *= scalar;
		y *= scalar;
	}

	void operator*=(Vector2D& other) {
		x *= other.x;
		y *= other.y;
	}

	void operator/=(double scalar) {
		x /= scalar;
		y /= scalar;
	}

	void operator/=(Vector2D& other) {
		x /= other.x;
		y /= other.y;
	}

	bool operator==(Vector2D& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) < epsilon && std::abs(y - other.y) < epsilon;
	}

	bool operator!=(Vector2D& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) > epsilon || std::abs(y - other.y) > epsilon;
	}
};

