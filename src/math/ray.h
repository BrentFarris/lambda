#ifndef LAMBDA_MATH_RAY_H
#define LAMBDA_MATH_RAY_H

#include "vector3.h"

struct Ray {
	Vector3 origin;
	Vector3 direction;

	Ray() = default;
	Ray(const Vector3& origin, const Vector3& direction)
		: origin(origin), direction(direction) {}

	bool triangle_hit(double rayLen, const Vector3& a,
		const Vector3& b, const Vector3& c) const;

	bool plane_hit(const Vector3& planePosition,
		const Vector3& planeNormal, Vector3& outHit) const
	{
		double d = Vector3::dot(planeNormal, direction);
		double epsilon = std::numeric_limits<double>::epsilon();
		if (std::abs(d) < epsilon)
			return false;
		const Vector3 diff = planePosition - origin;
		const double distance = Vector3::dot(diff, planeNormal) / d;
		if (distance <= 0.0)
			return false;
		outHit = point(distance);
		return true;
	}

	bool sphere_hit(const Vector3& center, double radius, double maxLen) const {
		Vector3 delta = center - origin;
		double len = Vector3::dot(direction, delta);
		if (len < 0.0 || len > (maxLen + radius))
			return false;
		double d2 = Vector3::dot(delta, delta) - (len * len);
		return d2 <= (radius * radius);
	}

	Vector3 point(double distance) const {
		Vector3 s = direction * distance;
		s += origin;
		return s;
	}
};

#endif
