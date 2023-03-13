export module Point2D;

import <limits>;
import <string>;
import <sstream>;
import <iostream>;
import <algorithm>;

export struct Point2D {
	double x;
	double y;

	Point2D() {
		x = 0.0;
		y = 0.0;
	}

	Point2D(double x, double y) {
		this->x = x;
		this->y = y;
	}

	Point2D(double xy) {
		this->x = xy;
		this->y = xy;
	}

	static Point2D zero() {
		return Point2D(0.0, 0.0);
	}

	static Point2D one() {
		return Point2D(1.0, 1.0);
	}

	static Point2D up() {
		return Point2D(0.0, 1.0);
	}

	static Point2D down() {
		return Point2D(0.0, -1.0);
	}

	static Point2D left() {
		return Point2D(-1.0, 0.0);
	}

	static Point2D right() {
		return Point2D(1.0, 0.0);
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

	Point2D normalized() {
		double mag = magnitude();
		return Point2D(x / mag, y / mag);
	}

	double dot(Point2D& other) {
		return x * other.x + y * other.y;
	}

	double distance(Point2D& target) {
		return sqrt(((target.x - x) * (target.x - x)) + ((target.y - y) * (target.y - y)));
	}

	double angle(Point2D& target) {
		return (atan2(x - target.x, y - target.y) * (180.0F / 3.14159265358979323846)) + 180.0;
	}

	Point2D min(Point2D& a, Point2D& b) {
		return Point2D(std::min(a.x, b.x), std::min(a.y, b.y));
	}

	Point2D max(Point2D& a, Point2D& b) {
		return Point2D(std::max(a.x, b.x), std::max(a.y, b.y));
	}

	Point2D lerp(Point2D& other, double t) {
		return Point2D(x + (other.x - x) * t, y + (other.y - y) * t);
	}

	std::string to_string() {
		std::stringstream ss;
		ss << std::to_string(x) << ", " << std::to_string(y);
		return std::move(ss.str());
	}

	void from_string(std::string& str) {
		char skip;
		std::stringstream ss;
		ss << str;
		ss >> x >> skip >> y;
	}

	void print() {
		std::cout << "Point2D<" << x << ", " << y << ">";
	}

	Point2D operator+(Point2D& other) {
		return Point2D(x + other.x, y + other.y);
	}

	Point2D operator-(Point2D& other) {
		return Point2D(x - other.x, y - other.y);
	}

	Point2D operator*(double scalar) {
		return Point2D(x * scalar, y * scalar);
	}

	Point2D operator*(Point2D& other) {
		return Point2D(x * other.x, y * other.y);
	}

	Point2D operator/(double scalar) {
		return Point2D(x / scalar, y / scalar);
	}

	Point2D operator/(Point2D other) {
		return Point2D(x / other.x, y / other.y);
	}

	void operator+=(Point2D& other) {
		x += other.x;
		y += other.y;
	}

	void operator-=(Point2D& other) {
		x -= other.x;
		y -= other.y;
	}

	void operator*=(double scalar) {
		x *= scalar;
		y *= scalar;
	}

	void operator*=(Point2D& other) {
		x *= other.x;
		y *= other.y;
	}

	void operator/=(double scalar) {
		x /= scalar;
		y /= scalar;
	}

	void operator/=(Point2D& other) {
		x /= other.x;
		y /= other.y;
	}

	bool operator==(Point2D& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) < epsilon && std::abs(y - other.y) < epsilon;
	}

	bool operator!=(Point2D& other) {
		double epsilon = std::numeric_limits<double>::epsilon();
		return std::abs(x - other.x) > epsilon || std::abs(y - other.y) > epsilon;
	}
};

