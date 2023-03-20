#ifndef LAMBDA_MATH_QUATERNION_H
#define LAMBDA_MATH_QUATERNION_H

#include <algorithm>
#include "vector3.h"
#include "vector4.h"
#include "geometry.h"
#include "matrix4x4.h"

namespace lambda {
	struct Quaternion {
		union {
			float linear[4];
			struct {
				float w;
				float x;
				float y;
				float z;
			};
		};

		Quaternion() : w(1.0), x(0.0), y(0.0), z(0.0) {}
		Quaternion(float w, float x, float y, float z) : w(w), x(x), y(y), z(z) {}
		Quaternion(float* wxyz) : w(wxyz[0]), x(wxyz[1]), y(wxyz[2]), z(wxyz[3]) {}
		Quaternion(Vector4& input) : w(input.w), x(input.x), y(input.y), z(input.z) {}

		Matrix4x4 to_matrix4x4() const {
			Matrix4x4 out{};
			float sqw = w * w;
			float sqx = x * x;
			float sqy = y * y;
			float sqz = z * z;
			// invs (inverse square length) is only required if quaternion is not already normalized
			float invs = 1.0F / (sqx + sqy + sqz + sqw);
			out.x0y0 = (sqx - sqy - sqz + sqw) * invs; // since sqw + sqx + sqy + sqz =1/invs*invs
			out.x1y1 = (-sqx + sqy - sqz + sqw) * invs;
			out.x2y2 = (-sqx - sqy + sqz + sqw) * invs;
			float tmp1 = x * y;
			float tmp2 = z * w;
			out.x1y0 = 2.0F * (tmp1 + tmp2) * invs;
			out.x0y1 = 2.0F * (tmp1 - tmp2) * invs;
			tmp1 = x * z;
			tmp2 = y * w;
			out.x2y0 = 2.0F * (tmp1 - tmp2) * invs;
			out.x0y2 = 2.0F * (tmp1 + tmp2) * invs;
			tmp1 = y * z;
			tmp2 = x * w;
			out.x2y1 = 2.0F * (tmp1 + tmp2) * invs;
			out.x1y2 = 2.0F * (tmp1 - tmp2) * invs;
			return out;
		}

		Vector3 to_euler() const {
			Vector3 out{};
			Matrix4x4 m = to_matrix4x4();
			out.y = rad2deg(asin(std::clamp(m.x0y2, -1.0F, 1.0F)));
			if (fabs(m.x0y2) < 0.9999999) {
				out.x = rad2deg(atan2(-m.x1y2, m.x2y2));
				out.z = rad2deg(atan2(-m.x0y1, m.x0y0));
			} else {
				out.x = rad2deg(atan2(m.x2y1, m.x1y1));
				out.z = 0.0;
			}
			return out;
		}

		void normalize() {
			Vector4::normalize_linear(linear);
		}

		void inverse() {
			float d = Vector4::dot_linear(linear, linear);
#if defined(USE_SIMD)
			vf32x4 l = vf32x4_set(linear[0], -linear[1], -linear[2], -linear[3]);
			vf32x4 r = vf32x4_set(d, d, d, d);
			vf32x4_st(linear, vf32x4_div(l, r));
#else
			w = w / d;
			x = -x / d;
			y = -y / d;
			z = -z / d;
#endif
		}

		void conjugate() {
#if defined(USE_SIMD)
			vf32x4 l = vf32x4_ld(linear);
			vf32x4 r = vf32x4_set(1.0F, -1.0F, -1.0F, -1.0F);
			vf32x4_st(linear, vf32x4_mul(l, r));
#else
			x = -x;
			y = -y;
			z = -z;
#endif
		}

		Vector3 multiply_vec3(const Vector3& rhs) {
			float v[12] = { x * 2.0F, y * 2.0F, z * 2.0F };
			v[3] = x * v[0];
			v[4] = y * v[1];
			v[5] = z * v[2];
			v[6] = x * v[1];
			v[7] = x * v[2];
			v[8] = y * v[2];
			v[9] = w * v[0];
			v[10] = w * v[1];
			v[11] = w * v[2];
			return Vector3(
				(1.0F - (v[4] + v[5])) * rhs.x + (v[6] - v[11]) * rhs.y + (v[7] + v[10]) * rhs.z,
				(v[6] + v[11]) * rhs.x + (1.0F - (v[3] + v[5])) * rhs.y + (v[8] - v[9]) * rhs.z,
				(v[7] - v[10]) * rhs.x + (v[8] + v[9]) * rhs.y + (1.0F - (v[3] + v[4])) * rhs.z);
		}

		std::string to_string() const {
			std::stringstream ss;
			ss << std::to_string(w) << ", " << std::to_string(x) << ", " << std::to_string(y) << ", " << std::to_string(z);
			return ss.str();
		}

		void from_string(const std::string& str) {
			char skip;
			std::stringstream ss;
			ss << str;
			ss >> w >> skip >> x >> skip >> y >> skip >> z;
		}

		void print() const {
			std::cout << "Quaternion<" << w << ", " << x << ", " << y << ", " << z << ">";
		}

		static Quaternion angle_axis(float angle, const Vector3& axis) {
			Vector3 cpy = axis * sin(angle * 0.5F);
			return Quaternion(cos(angle * 0.5F), cpy.x, cpy.y, cpy.z);
		}

		static Quaternion angle_between(const Vector3& lhs, const Vector3& rhs) {
			// It is important that the inputs are of equal length when
			// calculating the half-way vector.
			float k_cos_theta = Vector3::dot(lhs, rhs);
			float k = sqrt(pow(lhs.length(), 2.0F) * pow(rhs.length(), 2.0F));
			// TODO:  Approx here
			if (k_cos_theta / k == -1.0F) {
				// 180 degree rotation around any orthogonal vector
				const Vector3 o = lhs.orthogonal();
				Vector3 oNorm = o.normal();
				return Quaternion(0.0F, oNorm.x, oNorm.y, oNorm.z);
			}
			Vector3 c = Vector3::cross(lhs, rhs);
			Quaternion q = Quaternion(k_cos_theta + k, c.x, c.y, c.z);
			q.normalize();
			return q;
		}

		static Quaternion look_at(const Vector3& from, const Vector3& to) {
			const Vector3 diff = to - from;;
			Vector3 direction = diff.normal();
			const Vector3 back = Vector3::backward();
			float dot = Vector3::dot(back, direction);
			if (fabs(dot - (-1.0F)) < 0.000001F) {
				Vector3 u = Vector3::up();
				return angle_axis(rad2deg(M_PI), u);
			} else if (fabs(dot - (1.0F)) < 0.000001F)
				return Quaternion();
			float angle = -rad2deg(acos(dot));
			const Vector3 cross = Vector3::cross(back, direction);
			Vector3 nmlCross = cross.normal();
			return angle_axis(angle, nmlCross);
		}

		static Quaternion from_euler(const Vector3& from) {
			float x = deg2rad(from.x), y = deg2rad(from.y), z = deg2rad(from.z);
			float c1 = cos(x / 2.0F);
			float c2 = cos(y / 2.0F);
			float c3 = cos(z / 2.0F);
			float s1 = sin(x / 2.0F);
			float s2 = sin(y / 2.0F);
			float s3 = sin(z / 2.0F);
			return {
				c1 * c2 * c3 - s1 * s2 * s3,
				s1 * c2 * c3 + c1 * s2 * s3,
				c1 * s2 * c3 - s1 * c2 * s3,
				c1 * c2 * s3 + s1 * s2 * c3
			};
		}

		static Quaternion from_mat4(const Matrix4x4& mat) {
			const float* m = mat.linear;
			const float m00 = m[0], m10 = m[1], m20 = m[2],
				m01 = m[4], m11 = m[5], m21 = m[6],
				m02 = m[8], m12 = m[9], m22 = m[10];
			float t = m00 + m11 + m22;
			if (t > 0) {
				const float s = 0.5F / sqrt(t + 1.0F);
				return Quaternion(0.25F / s, (m12 - m21) * s, (m20 - m02) * s, (m01 - m10) * s);
			} else if (m00 > m11 && m00 > m22) {
				const float s = 2.0F * sqrt(1.0F + m00 - m11 - m22);
				return Quaternion((m12 - m21) / s, 0.25F * s, (m10 + m01) / s, (m20 + m02) / s);
			} else if (m11 > m22) {
				const float s = 2.0F * sqrt(1.0F + m11 - m00 - m22);
				return Quaternion((m20 - m02) / s, (m10 + m01) / s, 0.25F * s, (m21 + m12) / s);
			} else {
				const float s = 2.0F * sqrt(1.0F + m22 - m00 - m11);
				return Quaternion((m01 - m10) / s, (m20 + m02) / s, (m21 + m12) / s, 0.25F * s);
			}
		}

		static Quaternion quat_lerp(const Quaternion& from, const Quaternion& to, float factor) {
			Quaternion r;
			float t_ = 1.0F - factor;
			r.x = t_ * from.x + factor * to.x;
			r.y = t_ * from.y + factor * to.y;
			r.z = t_ * from.z + factor * to.z;
			r.w = t_ * from.w + factor * to.w;
			r.normalize();
			return r;
		}

		static Quaternion quat_slerp(const Quaternion& from, const Quaternion& to, float factor) {
			if (factor <= std::numeric_limits<float>::epsilon())
				return from;
			else if (factor >= 1.0F)
				return to;
			else {
				Quaternion r;
				const float x = from.x, y = from.y, z = from.z, w = from.w;
				float cosHalfTheta = w * to.w + x * to.x + y * to.y + z * to.z;
				if (cosHalfTheta < 0.0F) {
					r.w = -to.w;
					r.x = -to.x;
					r.y = -to.y;
					r.z = -to.z;
					cosHalfTheta = -cosHalfTheta;
				} else
					r = to;
				if (cosHalfTheta >= 1.0F) {
					r.w = w;
					r.x = x;
					r.y = y;
					r.z = z;
					return r;
				}
				const float sqrSinHalfTheta = 1.0F - cosHalfTheta * cosHalfTheta;
				if (sqrSinHalfTheta <= std::numeric_limits<float>::epsilon()) {
					const float s = 1.0F - factor;
					r.w = s * w + factor * r.w;
					r.x = s * x + factor * r.x;
					r.y = s * y + factor * r.y;
					r.z = s * z + factor * r.z;
					r.normalize();
					return r;
				}
				const float sinHalfTheta = sqrt(sqrSinHalfTheta);
				const float halfTheta = atan2(sinHalfTheta, cosHalfTheta);
				const float ratioA = sin((1.0F - factor) * halfTheta) / sinHalfTheta,
					ratioB = sin(factor * halfTheta) / sinHalfTheta;
				r.w = (w * ratioA + r.w * ratioB);
				r.x = (x * ratioA + r.x * ratioB);
				r.y = (y * ratioA + r.y * ratioB);
				r.z = (z * ratioA + r.z * ratioB);
				return r;
			}
		}

		Quaternion operator+(Quaternion& other) const {
#if defined(USE_SIMD)
			Quaternion ret = *this;
			vf32x4 l = vf32x4_ld(ret.linear);
			vf32x4 r = vf32x4_ld(other.linear);
			vf32x4_st(ret.linear, vf32x4_add(l, r));
			return ret;
#else
			return Quaternion(w + other.w, x + other.x, y + other.y, z + other.z);
#endif
		}

		Quaternion operator-(Quaternion& other) const {
#if defined(USE_SIMD)
			Quaternion ret = *this;
			vf32x4 l = vf32x4_ld(ret.linear);
			vf32x4 r = vf32x4_ld(other.linear);
			vf32x4_st(ret.linear, vf32x4_sub(l, r));
			return ret;
#else
			return Quaternion(w - other.w, x - other.x, y - other.y, z - other.z);
#endif
		}

		Quaternion operator*(float scalar) const {
#if defined(USE_SIMD)
			Quaternion ret = *this;
			vf32x4 l = vf32x4_ld(ret.linear);
			vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
			vf32x4_st(ret.linear, vf32x4_mul(l, r));
			return ret;
#else
			return Quaternion(w * scalar, x * scalar, y * scalar, z * scalar);
#endif
		}

		Quaternion operator*(Quaternion& other) const {
#if defined(USE_SIMD)
			Quaternion ret = *this;
			vf32x4 l = vf32x4_ld(ret.linear);
			vf32x4 r = vf32x4_ld(other.linear);
			vf32x4_st(ret.linear, vf32x4_mul(l, r));
			return ret;
#else
			return Quaternion(w * other.w, x * other.x, y * other.y, z * other.z);
#endif
		}

		Quaternion operator/(float scalar) const {
#if defined(USE_SIMD)
			Quaternion ret = *this;
			vf32x4 l = vf32x4_ld(ret.linear);
			vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
			vf32x4_st(ret.linear, vf32x4_div(l, r));
			return ret;
#else
			return Quaternion(w / scalar, x / scalar, y / scalar, z / scalar);
#endif
		}

		Quaternion operator/(Quaternion other) const {
#if defined(USE_SIMD)
			Quaternion ret = *this;
			vf32x4 l = vf32x4_ld(ret.linear);
			vf32x4 r = vf32x4_ld(other.linear);
			vf32x4_st(ret.linear, vf32x4_div(l, r));
			return ret;
#else
			return Quaternion(w / other.w, x / other.x, y / other.y, z / other.z);
#endif
		}

		void operator+=(Quaternion& other) {
#if defined(USE_SIMD)
			vf32x4 l = vf32x4_ld(linear);
			vf32x4 r = vf32x4_ld(other.linear);
			vf32x4_st(linear, vf32x4_add(l, r));
#else
			w += other.w;
			x += other.x;
			y += other.y;
			z += other.z;
#endif
		}

		void operator-=(Quaternion& other) {
#if defined(USE_SIMD)
			vf32x4 l = vf32x4_ld(linear);
			vf32x4 r = vf32x4_ld(other.linear);
			vf32x4_st(linear, vf32x4_sub(l, r));
#else
			w -= other.w;
			x -= other.x;
			y -= other.y;
			z -= other.z;
#endif
		}

		void operator*=(float scalar) {
#if defined(USE_SIMD)
			vf32x4 l = vf32x4_ld(linear);
			vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
			vf32x4_st(linear, vf32x4_mul(l, r));
#else
			w *= scalar;
			x *= scalar;
			y *= scalar;
			z *= scalar;
#endif
		}

		void operator*=(Quaternion& other) {
#if defined(USE_SIMD)
			vf32x4 l = vf32x4_ld(linear);
			vf32x4 r = vf32x4_ld(other.linear);
			vf32x4_st(linear, vf32x4_add(l, r));
#else
			w *= other.w;
			x *= other.x;
			y *= other.y;
			z *= other.z;
#endif
		}

		void operator/=(float scalar) {
#if defined(USE_SIMD)
			vf32x4 l = vf32x4_ld(linear);
			vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
			vf32x4_st(linear, vf32x4_div(l, r));
#else
			w /= scalar;
			x /= scalar;
			y /= scalar;
			z /= scalar;
#endif
		}

		void operator/=(Quaternion& other) {
#if defined(USE_SIMD)
			vf32x4 l = vf32x4_ld(linear);
			vf32x4 r = vf32x4_ld(other.linear);
			vf32x4_st(linear, vf32x4_div(l, r));
#else
			w /= other.w;
			x /= other.x;
			y /= other.y;
			z /= other.z;
#endif
		}

		bool operator==(Quaternion& other) const {
			float epsilon = std::numeric_limits<float>::epsilon();
			float d[4];
#if defined(USE_SIMD)
			vf32x4 l = vf32x4_ld(linear);
			vf32x4 r = vf32x4_ld(other.linear);
			vf32x4_st(d, vf32x4_sub(l, r));
#else
			d[0] = w - other.w;
			d[1] = x - other.x;
			d[2] = y - other.y;
			d[3] = z - other.z;
#endif
			// TODO:  Could probably use the abs of the sum rather than each
			return std::abs(d[0]) < epsilon
				&& std::abs(d[1]) < epsilon
				&& std::abs(d[2]) < epsilon
				&& std::abs(d[3]) < epsilon;
		}

		bool operator!=(Quaternion& other) const {
			float epsilon = std::numeric_limits<float>::epsilon();
			float d[4];
#if defined(USE_SIMD)
			vf32x4 l = vf32x4_ld(linear);
			vf32x4 r = vf32x4_ld(other.linear);
			vf32x4_st(d, vf32x4_sub(l, r));
#else
			d[0] = w - other.w;
			d[1] = x - other.x;
			d[2] = y - other.y;
			d[3] = z - other.z;
#endif
			// TODO:  Could probably use the abs of the sum rather than each
			return std::abs(d[0]) > epsilon
				&& std::abs(d[1]) > epsilon
				&& std::abs(d[2]) > epsilon
				&& std::abs(d[3]) > epsilon;
		}
	};
}

#endif
