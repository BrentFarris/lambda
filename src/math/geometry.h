#ifndef LAMBDA_MATH_GEOMETRY_H
#define LAMBDA_MATH_GEOMETRY_H

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

template <class T>
static T rad2deg(T radian) {
	return radian * (180.0 / M_PI);
}

template <class T>
static T deg2rad(T degree) {
	return degree * (M_PI / 180.0);
}

#endif
