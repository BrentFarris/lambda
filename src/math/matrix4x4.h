#ifndef LAMBDA_MATH_MATRIX4X4_H
#define LAMBDA_MATH_MATRIX4X4_H

#include "vector3.h"
#include "vector4.h"
#include "geometry.h"
#include "hardware_vector.h"

enum class Mat4Row {
	ROW_0 = 0,
	ROW_1 = 1,
	ROW_2 = 2,
	ROW_3 = 3
};

enum class Mat4Col {
	COL_0 = 0,
	COL_1 = 1,
	COL_2 = 2,
	COL_3 = 3
};

struct Matrix4x4 {
	union {
		double linear[16];
		struct {
			double x0y0;
			double x1y0;
			double x2y0;
			double x3y0;
			double x0y1;
			double x1y1;
			double x2y1;
			double x3y1;
			double x0y2;
			double x1y2;
			double x2y2;
			double x3y2;
			double x0y3;
			double x1y3;
			double x2y3;
			double x3y3;
		};
	};

	Matrix4x4() :
		x0y0(1.0), x1y0(0.0), x2y0(0.0), x3y0(0.0),
		x0y1(0.0), x1y1(1.0), x2y1(0.0), x3y1(0.0),
		x0y2(0.0), x1y2(0.0), x2y2(1.0), x3y2(0.0),
		x0y3(0.0), x1y3(0.0), x2y3(0.0), x3y3(1.0) {}
	Matrix4x4(double setAll) :
		x0y0(setAll), x1y0(setAll), x2y0(setAll), x3y0(setAll),
		x0y1(setAll), x1y1(setAll), x2y1(setAll), x3y1(setAll),
		x0y2(setAll), x1y2(setAll), x2y2(setAll), x3y2(setAll),
		x0y3(setAll), x1y3(setAll), x2y3(setAll), x3y3(setAll) {}
	Matrix4x4(double x0y0, double x1y0, double x2y0, double x3y0,
		double x0y1, double x1y1, double x2y1, double x3y1,
		double x0y2, double x1y2, double x2y2, double x3y2,
		double x0y3, double x1y3, double x2y3, double x3y3) :
			x0y0(x0y0), x1y0(x1y0), x2y0(x2y0), x3y0(x3y0),
			x0y1(x0y1), x1y1(x1y1), x2y1(x2y1), x3y1(x3y1),
			x0y2(x0y2), x1y2(x1y2), x2y2(x2y2), x3y2(x3y2),
		x0y3(x0y3), x1y3(x1y3), x2y3(x2y3), x3y3(x3y3) {}
	Matrix4x4(double* other) :
			x0y0(other[0]), x1y0(other[1]), x2y0(other[2]), x3y0(other[3]),
			x0y1(other[4]), x1y1(other[5]), x2y1(other[6]), x3y1(other[7]),
			x0y2(other[8]), x1y2(other[9]), x2y2(other[10]), x3y2(other[11]),
			x0y3(other[12]), x1y3(other[13]), x2y3(other[14]), x3y3(other[15]) {}

	double* at(Mat4Row rowIndex, Mat4Col colIndex) const {
		return (double*)linear + (int(rowIndex) * (size_t)4) + int(colIndex);
	}

	void transpose() {
		Matrix4x4 result = *this;
		x0y1 = result.x1y0;
		x0y2 = result.x2y0;
		x0y3 = result.x3y0;
		x1y0 = result.x0y1;
		x1y2 = result.x2y1;
		x1y3 = result.x3y1;
		x2y0 = result.x0y2;
		x2y1 = result.x1y2;
		x2y3 = result.x3y2;
		x3y0 = result.x0y3;
		x3y1 = result.x1y3;
		x3y2 = result.x2y3;
	}

	Matrix4x4 transposed() {
		Matrix4x4 m;
		m.x0y0 = x0y0;
		m.x1y0 = x0y1;
		m.x2y0 = x0y2;
		m.x3y0 = x0y3;
		m.x0y1 = x1y0;
		m.x1y1 = x1y1;
		m.x2y1 = x1y2;
		m.x3y1 = x1y3;
		m.x0y2 = x2y0;
		m.x1y2 = x2y1;
		m.x2y2 = x2y2;
		m.x3y2 = x2y3;
		m.x0y3 = x3y0;
		m.x1y3 = x3y1;
		m.x2y3 = x3y2;
		m.x3y3 = x3y3;
		return m;
	}

	void reset() {
		x0y0 = 1.0F;
		x1y0 = 0.0F;
		x2y0 = 0.0F;
		x3y0 = 0.0F;
		x0y1 = 0.0F;
		x1y1 = 1.0F;
		x2y1 = 0.0F;
		x3y1 = 0.0F;
		x0y2 = 0.0F;
		x1y2 = 0.0F;
		x2y2 = 1.0F;
		x3y2 = 0.0F;
		x0y3 = 0.0F;
		x1y3 = 0.0F;
		x2y3 = 0.0F;
		x3y3 = 1.0F;
	}

	void add(const Matrix4x4& rhs) {
#ifdef USE_SIMD
		vf64x4 a0 = vf64x4_set(x0y0, x1y0, x2y0, x3y0);
		vf64x4 a1 = vf64x4_set(x0y1, x1y1, x2y1, x3y1);
		vf64x4 a2 = vf64x4_set(x0y2, x1y2, x2y2, x3y2);
		vf64x4 a3 = vf64x4_set(x0y3, x1y3, x2y3, x3y3);
		vf64x4 b0 = vf64x4_set(rhs.x0y0, rhs.x1y0, rhs.x2y0, rhs.x3y0);
		vf64x4 b1 = vf64x4_set(rhs.x0y1, rhs.x1y1, rhs.x2y1, rhs.x3y1);
		vf64x4 b2 = vf64x4_set(rhs.x0y2, rhs.x1y2, rhs.x2y2, rhs.x3y2);
		vf64x4 b3 = vf64x4_set(rhs.x0y3, rhs.x1y3, rhs.x2y3, rhs.x3y3);
		vf64x4n_add(resA0B0, a0, b0);
		vf64x4_st((double*)linear, resA0B0);
		vf64x4n_add(resA1B1, a1, b1);
		vf64x4_st((double*)linear + 4, resA1B1);
		vf64x4n_add(resA2B2, a2, b2);
		vf64x4_st((double*)linear + 8, resA2B2);
		vf64x4n_add(resA3B3, a3, b3);
		vf64x4_st((double*)linear + 12, resA3B3);
#else
		x0y0 += rhs->x0y0;
		x1y0 += rhs->x1y0;
		x2y0 += rhs->x2y0;
		x3y0 += rhs->x3y0;
		x0y1 += rhs->x0y1;
		x1y1 += rhs->x1y1;
		x2y1 += rhs->x2y1;
		x3y1 += rhs->x3y1;
		x0y2 += rhs->x0y2;
		x1y2 += rhs->x1y2;
		x2y2 += rhs->x2y2;
		x3y2 += rhs->x3y2;
		x0y3 += rhs->x0y3;
		x1y3 += rhs->x1y3;
		x2y3 += rhs->x2y3;
		x3y3 += rhs->x3y3;
#endif
	}

	void subtract(const Matrix4x4& rhs) {
#ifdef USE_SIMD
		vf64x4 a0 = vf64x4_set(x0y0, x1y0, x2y0, x3y0);
		vf64x4 a1 = vf64x4_set(x0y1, x1y1, x2y1, x3y1);
		vf64x4 a2 = vf64x4_set(x0y2, x1y2, x2y2, x3y2);
		vf64x4 a3 = vf64x4_set(x0y3, x1y3, x2y3, x3y3);
		vf64x4 b0 = vf64x4_set(rhs.x0y0, rhs.x1y0, rhs.x2y0, rhs.x3y0);
		vf64x4 b1 = vf64x4_set(rhs.x0y1, rhs.x1y1, rhs.x2y1, rhs.x3y1);
		vf64x4 b2 = vf64x4_set(rhs.x0y2, rhs.x1y2, rhs.x2y2, rhs.x3y2);
		vf64x4 b3 = vf64x4_set(rhs.x0y3, rhs.x1y3, rhs.x2y3, rhs.x3y3);
		vf64x4n_sub(resA0B0, a0, b0);
		vf64x4_st((double*)linear, resA0B0);
		vf64x4n_sub(resA1B1, a1, b1);
		vf64x4_st((double*)linear + 4, resA1B1);
		vf64x4n_sub(resA2B2, a2, b2);
		vf64x4_st((double*)linear + 8, resA2B2);
		vf64x4n_sub(resA3B3, a3, b3);
		vf64x4_st((double*)linear + 12, resA3B3);
#else
		x0y0 -= rhs->x0y0;
		x1y0 -= rhs->x1y0;
		x2y0 -= rhs->x2y0;
		x3y0 -= rhs->x3y0;
		x0y1 -= rhs->x0y1;
		x1y1 -= rhs->x1y1;
		x2y1 -= rhs->x2y1;
		x3y1 -= rhs->x3y1;
		x0y2 -= rhs->x0y2;
		x1y2 -= rhs->x1y2;
		x2y2 -= rhs->x2y2;
		x3y2 -= rhs->x3y2;
		x0y3 -= rhs->x0y3;
		x1y3 -= rhs->x1y3;
		x2y3 -= rhs->x2y3;
		x3y3 -= rhs->x3y3;
#endif
	}

	void negate() {
#ifdef USE_SIMD
		vf64x4 a0 = vf64x4_set(x0y0, x1y0, x2y0, x3y0);
		vf64x4 a1 = vf64x4_set(x0y1, x1y1, x2y1, x3y1);
		vf64x4 a2 = vf64x4_set(x0y2, x1y2, x2y2, x3y2);
		vf64x4 a3 = vf64x4_set(x0y3, x1y3, x2y3, x3y3);
		vf64x4 n = vf64x4_set(-1.0, -1.0, -1.0, -1.0);
		vf64x4n_mul(resA0B0, a0, n);
		vf64x4_st((double*)linear, resA0B0);
		vf64x4n_mul(resA1B1, a1, n);
		vf64x4_st((double*)linear + 4, resA1B1);
		vf64x4n_mul(resA2B2, a2, n);
		vf64x4_st((double*)linear + 8, resA2B2);
		vf64x4n_mul(resA3B3, a3, n);
		vf64x4_st((double*)linear + 12, resA3B3);
#else
		x0y0 *= -1.0;
		x1y0 *= -1.0;
		x2y0 *= -1.0;
		x3y0 *= -1.0;
		x0y1 *= -1.0;
		x1y1 *= -1.0;
		x2y1 *= -1.0;
		x3y1 *= -1.0;
		x0y2 *= -1.0;
		x1y2 *= -1.0;
		x2y2 *= -1.0;
		x3y2 *= -1.0;
		x0y3 *= -1.0;
		x1y3 *= -1.0;
		x2y3 *= -1.0;
		x3y3 *= -1.0;
#endif
	}

	Vector3 project(const Matrix4x4& mvp, const Vector3& pos, const Vector4& viewport) const {
		Vector4 pos4 = Vector4(pos.x, pos.y, pos.z, 1.0);
		pos4 = mvp * pos4;
		double z = pos4.z;
		pos4 /= pos4.w;
		pos4 *= 0.5;
		pos4 += Vector4(0.5);
		return Vector3(pos4.x * viewport.z + viewport.x,
			pos4.y * viewport.w + viewport.y, z);
	}

	Vector3 unproject(const Matrix4x4& invMat, const Vector3& pos, const Vector4& viewport) {
		Vector4 v;
		v.x = 2.0 * (pos.x - viewport.x) / viewport.z - 1.0;
		v.y = 2.0 * (pos.y - viewport.y) / viewport.w - 1.0;
		v.z = 2.0 * pos.z - 1.0;
		v.w = 1.0;
		v = invMat * v;
		v *= 1.0 / v.w;
		return Vector3(v.x, v.y, v.z);
	}

	Vector3 position() const {
		return Vector3(x0y3, x1y3, x2y3);
	}

	void translate(const Vector3& translation) {
#ifdef USE_SIMD
		vf64x4 a = vf64x4_ld(&x0y3);
		vf64x4 b = vf64x4_set(translation.x, translation.y, translation.z, 0.0);
		vf64x4 sum = vf64x4_add(a, b);
		vf64x4_st(&x0y3, sum);
#else
		matrix4x4->x0y3 += translation->x;
		matrix4x4->x1y3 += translation->y;
		matrix4x4->x2y3 += translation->z;
#endif
	}

	void set_translation(const Vector3& translation) {
		x0y3 = translation.x;
		x1y3 = translation.y;
		x2y3 = translation.z;
	}

	void scale(const Vector3& scale) {
		x0y0 *= scale.x;
		x1y1 *= scale.y;
		x2y2 *= scale.z;
	}

	void look_at(const Vector3& eye, const Vector3& center, const Vector3& up) {
		Vector3 f = eye - center;
		f.normalize();
		Vector3 s = Vector3::cross(up, f);
		s.normalize();
		Vector3 u = Vector3::cross(f, s);
		Vector3 ns = s.negative();
		Vector3 nu = u.negative();
		Vector3 nf = f.negative();
		*this = {
			s.x, u.x, f.x, 0.0,
			s.y, u.y, f.y, 0.0,
			s.z, u.z, f.z, 0.0,
			Vector3::dot(ns, eye), Vector3::dot(nu, eye), Vector3::dot(nf, eye), 1.0
		};
	}

	void rotate(const Vector3& rotate);

	void rotate_x(double angles) {
		Matrix4x4 rot{};
		rot.x1y1 = cos(deg2rad(angles));
		rot.x2y1 = -sin(deg2rad(angles));
		rot.x1y2 = sin(deg2rad(angles));
		rot.x2y2 = cos(deg2rad(angles));
		*this *= rot;
	}

	void rotate_y(double angles) {
		Matrix4x4 rot{};
		rot.x0y0 = cos(deg2rad(angles));
		rot.x2y0 = sin(deg2rad(angles));
		rot.x0y2 = -sin(deg2rad(angles));
		rot.x2y2 = cos(deg2rad(angles));
		*this *= rot;
	}

	void rotate_z(double angles) {
		Matrix4x4 rot{};
		rot.x0y0 = cos(deg2rad(angles));
		rot.x1y0 = -sin(deg2rad(angles));
		rot.x0y1 = sin(deg2rad(angles));
		rot.x1y1 = cos(deg2rad(angles));
		*this *= rot;
	}

	void rotate_angles(double angle, const Vector3& axis) {
		double a = angle;
		double c = cos(a);
		double s = sin(a);
		Vector3 axisNorm = axis.normal();
		Vector3 temp = axisNorm * (1.0 - c);
		//vec<3, T, Q> temp((T(1) - c) * axis);
		Matrix4x4 rot{};
		rot.x0y0 = c + temp.x * axisNorm.x;
		rot.x0y1 = temp.x * axisNorm.y + s * axisNorm.z;
		rot.x0y2 = temp.x * axisNorm.z - s * axisNorm.y;
		rot.x1y0 = temp.y * axisNorm.x - s * axisNorm.z;
		rot.x1y1 = c + temp.y * axisNorm.y;
		rot.x1y2 = temp.y * axisNorm.z + s * axisNorm.x;
		rot.x2y0 = temp.z * axisNorm.x + s * axisNorm.y;
		rot.x2y1 = temp.z * axisNorm.y - s * axisNorm.x;
		rot.x2y2 = c + temp.z * axisNorm.z;
		const Vector4* matCols = (Vector4*)linear;
		Matrix4x4 res;
		const Vector4 c0x0y0 = matCols[0] * &rot.x0y0;
		const Vector4 c0x1y0 = matCols[0] * &rot.x1y0;
		const Vector4 c0x2y0 = matCols[0] * &rot.x2y0;
		const Vector4 c1x0y1 = matCols[1] * &rot.x0y1;
		const Vector4 c1x1y1 = matCols[1] * &rot.x1y1;
		const Vector4 c1x2y1 = matCols[1] * &rot.x2y1;
		const Vector4 c2x0y2 = matCols[2] * &rot.x0y2;
		const Vector4 c2x1y2 = matCols[2] * &rot.x1y2;
		const Vector4 c2x2y2 = matCols[2] * &rot.x2y2;
		const Vector4 r0 = c0x0y0 + c1x0y1;
		const Vector4 r1 = c0x1y0 + c1x1y1;
		const Vector4 r2 = c0x2y0 + c1x2y1;
		((Vector4*)&res)[0] = r0 + c2x0y2;
		((Vector4*)&res)[1] = r1 + c2x1y2;
		((Vector4*)&res)[2] = r2 + c2x2y2;
		((Vector4*)&res)[3] = matCols[3];
		*this = res;
	}

	void inverse() {
		double t[6] = { 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F };
		double a = x0y0, b = x0y1, c = x0y2, d = x0y3,
			e = x1y0, f = x1y1, g = x1y2, h = x1y3,
			i = x2y0, j = x2y1, k = x2y2, l = x2y3,
			m = x3y0, n = x3y1, o = x3y2, p = x3y3;
		t[0] = k * p - o * l; t[1] = j * p - n * l; t[2] = j * o - n * k;
		t[3] = i * p - m * l; t[4] = i * o - m * k; t[5] = i * n - m * j;
		x0y0 = f * t[0] - g * t[1] + h * t[2];
		x1y0 = -(e * t[0] - g * t[3] + h * t[4]);
		x2y0 = e * t[1] - f * t[3] + h * t[5];
		x3y0 = -(e * t[2] - f * t[4] + g * t[5]);
		x0y1 = -(b * t[0] - c * t[1] + d * t[2]);
		x1y1 = a * t[0] - c * t[3] + d * t[4];
		x2y1 = -(a * t[1] - b * t[3] + d * t[5]);
		x3y1 = a * t[2] - b * t[4] + c * t[5];
		t[0] = g * p - o * h; t[1] = f * p - n * h; t[2] = f * o - n * g;
		t[3] = e * p - m * h; t[4] = e * o - m * g; t[5] = e * n - m * f;
		x0y2 = b * t[0] - c * t[1] + d * t[2];
		x1y2 = -(a * t[0] - c * t[3] + d * t[4]);
		x2y2 = a * t[1] - b * t[3] + d * t[5];
		x3y2 = -(a * t[2] - b * t[4] + c * t[5]);
		t[0] = g * l - k * h; t[1] = f * l - j * h; t[2] = f * k - j * g;
		t[3] = e * l - i * h; t[4] = e * k - i * g; t[5] = e * j - i * f;
		x0y3 = -(b * t[0] - c * t[1] + d * t[2]);
		x1y3 = a * t[0] - c * t[3] + d * t[4];
		x2y3 = -(a * t[1] - b * t[3] + d * t[5]);
		x3y3 = a * t[2] - b * t[4] + c * t[5];
		double det = 1.0 / (a * x0y0 + b * x1y0 + c * x2y0 + d * x3y0);
		x0y0 *= det; x0y1 *= det; x0y2 *= det; x0y3 *= det;
		x1y0 *= det; x1y1 *= det; x1y2 *= det; x1y3 *= det;
		x2y0 *= det; x2y1 *= det; x2y2 *= det; x2y3 *= det;
		x3y0 *= det; x3y1 *= det; x3y2 *= det; x3y3 *= det;
	}

	void transform_point(Vector3& point) const {
		Vector4 pt0 = Vector4(point.x, point.y, point.z, 1.0);
		Vector4 res = *this * pt0;
		point.x = res.x;
		point.y = res.y;
		point.z = res.z;
		point /= res.w;
	}

	Vector3 right() const {
		return Vector3(x0y0, x1y0, x2y0);
	}

	Vector3 up() const {
		return Vector3(x0y1, x1y1, x2y1);
	}

	Vector3 forward() const {
		return Vector3(x0y2, x1y2, x2y2);
	}

	static Matrix4x4 orthographic(double left, double right,
		double bottom, double top, double near, double far)
	{
		Matrix4x4 mat{};
		mat.x0y0 = 2.0 / (right - left);
		mat.x1y1 = -2.0 / (top - bottom);
		mat.x2y2 = -1.0 / (far - near);
		mat.x0y3 = -(right + left) / (right - left);
		mat.x1y3 = -(top + bottom) / (top - bottom);
		mat.x2y3 = -near / (far - near);
		mat.x3y3 = 1.0;
		return mat;
	}

	static Matrix4x4 perspective(double fovy,
		double aspect, double nearVal, double farVal)
	{
		double f, fn;
		Matrix4x4 mat{};
		f = 1.0 / tan(fovy * 0.5);
		fn = 1.0 / (nearVal - farVal);
		mat.x0y0 = f / aspect;
		mat.x1y1 = -f;
		mat.x2y2 = (nearVal + farVal) * fn;
		mat.x3y2 = -1.0;
		mat.x2y3 = 2.0 * nearVal * farVal * fn;
		return mat;
	}

	static Matrix4x4 projection_gl2vulkan(const Matrix4x4& mat) {
		Matrix4x4 res = mat;
		res.x1y1 *= -1.0F;
		return res;
	}

	static Vector4 row_vector(const Matrix4x4& matrix4x4, Mat4Row rowIndex) {
		return Vector4((double*)matrix4x4.linear + int(rowIndex) * 4);
	}

	static Vector4 col_vector(const Matrix4x4& matrix4x4, Mat4Col colIndex) {
#ifdef USE_SIMD
		vi32x4 cols = vi32x4_set(int(colIndex), int(colIndex), int(colIndex), int(colIndex));
		vi32x4 offs = vi32x4_set(0, 4, 8, 12);
		union { vi32x4 v; int i[4]; } res;
		res.v = vi32x4_add(cols, offs);
		return Vector4 {
			*(((double*)matrix4x4.linear) + res.i[0]),
				* (((double*)matrix4x4.linear) + res.i[1]),
				* (((double*)matrix4x4.linear) + res.i[2]),
				* (((double*)matrix4x4.linear) + res.i[3])
		};
#else
		return Vector4 {
			*(((double*)matrix4x4.linear) + int(colIndex)),
				* (((double*)matrix4x4.linear) + int(colIndex) + 4),
				* (((double*)matrix4x4.linear) + int(colIndex) + 8),
				* (((double*)matrix4x4.linear) + int(colIndex) + 12)
		};
#endif
	}

	Matrix4x4 operator*(const Matrix4x4& rhs) const {
		Matrix4x4 copy = *this;
		copy *= rhs;
		return copy;
	}

	void operator*=(const Matrix4x4& rhs) {
		Matrix4x4 copy = *this;
		Vector4 row, col;
		size_t matrixIndex = 0;
		for (size_t i = 0; i < 4; ++i) {
			row = row_vector(copy, (Mat4Row)i);
			for (size_t j = 0; j < 4; ++j) {
				col = col_vector(rhs, (Mat4Col)j);
				*(((double*)linear) + matrixIndex++) = Vector4::dot(row, col);
			}
		}
	}

	friend Vector4 operator*(const Vector4& lhs, const Matrix4x4& rhs) {
		Vector4 result;
		Vector4 row;
		for (size_t i = 0; i < 4; ++i) {
			row = Matrix4x4::row_vector(rhs, (Mat4Row)i);
			*((double*)&result + i) = Vector4::dot(row, lhs);
		}
		return result;
	}

	friend Vector4 operator*(const Matrix4x4& lhs, const Vector4& rhs) {
		Vector4 result;
		Vector4 row;
		for (size_t i = 0; i < 4; ++i) {
			row = Matrix4x4::col_vector(lhs, (Mat4Col)i);
			*((double*)&result + i) = Vector4::dot(row, rhs);
		}
		return result;
	}
};

#endif
