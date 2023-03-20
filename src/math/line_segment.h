#ifndef LAMBDA_MATH_LINE_SEGMENT_H
#define LAMBDA_MATH_LINE_SEGMENT_H

#include "ray.h"
#include "vector3.h"

namespace lambda {
	struct LineSegment {
		Vector3 a;
		Vector3 b;

		LineSegment() = default;
		LineSegment(const Vector3& a, const Vector3& b) : a(a), b(b) {}

		bool triangle_hit(const Vector3& p0, const Vector3& p1, const Vector3& p2) const {
			const Vector3 p = a;
			const Vector3 q = b;
			const Vector3 ab = p1 - p0;
			const Vector3 ac = p2 - p0;
			const Vector3 qp = p - q;

			// Compute triangle normal. Can be pre-calculated or cached if
			// intersecting multiple segments against the same triangle
			Vector3 n = Vector3::cross(ab, ac);

			// Compute denominator d. If d <= 0, segment is parallel to or points
			// away from triangle, so exit early
			float d = Vector3::dot(qp, n);
			if (d <= 0.0)
				return false;

			// Compute intersection t value of pq with plane of triangle. A ray
			// intersects iff 0 <= t. Segment intersects iff 0 <= t <= 1. Delay
			// dividing by d until intersection has been found to pierce triangle
			Vector3 ap = p - p0;
			float t = Vector3::dot(ap, n);
			if (t < 0.0)
				return false;
			//if (t > d)  // For segment; exclude this code line for a ray test
			//	return false;

			// Compute barycentric coordinate components and test if within bounds
			Vector3 e = Vector3::cross(qp, ap);
			float v = Vector3::dot(ac, e);
			if (v < 0.0 || v > d)
				return false;
			float w = -Vector3::dot(ab, e);
			if (w < 0.0 || v + w > d)
				return false;

			// Segment/ray intersects triangle. Perform delayed division and
			// compute the last barycentric coordinate component
			// This can be uncommented if we wish to get u,v,w and t as outs to this func
			//float ood = 1.0F / d;
			//t *= ood;
			//v *= ood;
			//w *= ood;
			//float u = 1.0F - v - w;
			return true;
		}

		static LineSegment from_ray(const Ray& r, float len) {
			return LineSegment(r.origin, r.point(len));
		}
	};
}

#endif
