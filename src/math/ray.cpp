#include "ray.h"
#include "line_segment.h"

bool Ray::triangle_hit(double rayLen, const Vector3& a,
	const Vector3& b, const Vector3& c) const
{
	LineSegment s = { origin, r.point(rayLen) };
	return s.triangle_hit(a, b, c);
}