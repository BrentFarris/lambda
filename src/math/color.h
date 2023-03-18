#ifndef LAMBDA_MATH_COLOR_H
#define LAMBDA_MATH_COLOR_H

#include <string>
#include <sstream>
#include <iostream>
#include "vector4.h"
#include "hardware_vector.h"

struct Color {
	union {
		float linear[4];
		struct {
			float r;
			float g;
			float b;
			float a;
		};
	};

	Color() : r(1.0F), g(1.0F), b(1.0F), a(1.0F) {}
	Color(float r, float g, float b, float a) : r(r), g(g), b(b), a(a) {}
	Color(float r, float g, float b) : r(r), g(g), b(b), a(1.0F) {}
	Color(float rgb) : r(rgb), g(rgb), b(rgb), a(1.0F) {}
	Color(float* rgba) : r(rgba[0]), g(rgba[1]), b(rgba[2]), a(rgba[3]) {}
	Color(const Vector4& vec) : r(vec.x), g(vec.y), b(vec.z), a(vec.w) {}
	Color(int r, int g, int b, int a) : r(r / 255.0F), g(g / 255.0F), b(b / 255.0F), a(a / 255.0F) {}
	Color(int rgb) : r(rgb / 255.0F), g(rgb / 255.0F), b(rgb / 255.0F), a(1.0F) {}

	static Color red() { return Color(1.0F, 0.0F, 0.0F); }
	static Color white() { return Color(1.0F); }
	static Color blue() { return Color(0.0F, 0.0F, 1.0F); }
	static Color black() { return Color(0.0F); }
	static Color green() { return Color(0.0F, 1.0F, 1.0F); }
	static Color yellow() { return Color(1.0F, 1.0F, 0.0F); }
	static Color orange() { return Color(1.0F, 0.647F, 0.0F); }
	static Color clear() { return Color(0.0F, 0.0F, 0.0F, 0.0F); }
	static Color gray() { return Color(0.5F); }
	// Begin generated colors
	static Color purple() { return Color(0.5F, 0.0F, 0.5F); }
	static Color brown() { return Color(0.647F, 0.165F, 0.165F); }
	static Color pink() { return Color(1.0F, 0.753F, 0.796F); }
	static Color cyan() { return Color(0.0F, 1.0F, 1.0F); }
	static Color magenta() { return Color(1.0F, 0.0F, 1.0F); }
	static Color teal() { return Color(0.0F, 0.5F, 0.5F); }
	static Color lime() { return Color(0.0F, 1.0F, 0.0F); }
	static Color maroon() { return Color(0.5F, 0.0F, 0.0F); }
	static Color olive() { return Color(0.5F, 0.5F, 0.0F); }
	static Color navy() { return Color(0.0F, 0.0F, 0.5F); }
	static Color silver() { return Color(0.753F, 0.753F, 0.753F); }
	static Color gold() { return Color(1.0F, 0.843F, 0.0F); }
	static Color sky() { return Color(0.529F, 0.808F, 0.922F); }
	static Color violet() { return Color(0.933F, 0.51F, 0.933F); }
	static Color indigo() { return Color(0.294F, 0.0F, 0.51F); }
	static Color turquoise() { return Color(0.251F, 0.878F, 0.816F); }
	static Color azure() { return Color(0.941F, 1.0F, 1.0F); }
	static Color chartreuse() { return Color(0.498F, 1.0F, 0.0F); }
	static Color coral() { return Color(1.0F, 0.498F, 0.314F); }
	static Color crimson() { return Color(0.863F, 0.078F, 0.235F); }
	static Color fuchsia() { return Color(1.0F, 0.0F, 1.0F); }
	static Color khaki() { return Color(0.941F, 0.902F, 0.549F); }
	static Color lavender() { return Color(0.902F, 0.902F, 0.98F); }
	static Color moccasin() { return Color(1.0F, 0.894F, 0.71F); }
	static Color salmon() { return Color(0.98F, 0.502F, 0.447F); }
	static Color sienna() { return Color(0.627F, 0.322F, 0.176F); }
	static Color tan() { return Color(0.824F, 0.706F, 0.549F); }
	static Color tomato() { return Color(1.0F, 0.388F, 0.278F); }
	static Color wheat() { return Color(0.961F, 0.871F, 0.702F); }
	static Color aqua() { return Color(0.0F, 1.0F, 1.0F); }
	static Color aquamarine() { return Color(0.498F, 1.0F, 0.831F); }
	static Color beige() { return Color(0.961F, 0.961F, 0.863F); }
	static Color bisque() { return Color(1.0F, 0.894F, 0.769F); }
	static Color blanchedalmond() { return Color(1.0F, 0.922F, 0.804F); }
	static Color blueviolet() { return Color(0.541F, 0.169F, 0.886F); }
	static Color burlywood() { return Color(0.871F, 0.722F, 0.529F); }
	static Color cadetblue() { return Color(0.373F, 0.62F, 0.627F); }
	static Color chocolate() { return Color(0.824F, 0.412F, 0.118F); }
	static Color cornflowerblue() { return Color(0.392F, 0.584F, 0.929F); }
	static Color cornsilk() { return Color(1.0F, 0.973F, 0.863F); }
	static Color darkblue() { return Color(0.0F, 0.0F, 0.545F); }
	static Color darkcyan() { return Color(0.0F, 0.545F, 0.545F); }
	static Color darkgoldenrod() { return Color(0.722F, 0.525F, 0.043F); }
	static Color darkgray() { return Color(0.663F, 0.663F, 0.663F); }
	static Color darkgreen() { return Color(0.0F, 0.392F, 0.0F); }
	static Color darkkhaki() { return Color(0.741F, 0.718F, 0.42F); }
	static Color darkmagenta() { return Color(0.545F, 0.0F, 0.545F); }
	static Color darkolivegreen() { return Color(0.333F, 0.42F, 0.184F); }
	static Color darkorange() { return Color(1.0F, 0.549F, 0.0F); }
	static Color darkorchid() { return Color(0.6F, 0.196F, 0.8F); }
	static Color darkred() { return Color(0.545F, 0.0F, 0.0F); }
	static Color darksalmon() { return Color(0.914F, 0.588F, 0.478F); }
	static Color darkseagreen() { return Color(0.561F, 0.737F, 0.561F); }
	static Color darkslateblue() { return Color(0.282F, 0.239F, 0.545F); }
	static Color darkslategray() { return Color(0.184F, 0.31F, 0.31F); }
	static Color darkturquoise() { return Color(0.0F, 0.808F, 0.82F); }
	static Color darkviolet() { return Color(0.58F, 0.0F, 0.827F); }
	static Color deeppink() { return Color(1.0F, 0.078F, 0.576F); }
	static Color deepskyblue() { return Color(0.0F, 0.749F, 1.0F); }
	static Color dimgray() { return Color(0.412F, 0.412F, 0.412F); }
	static Color dodgerblue() { return Color(0.118F, 0.565F, 1.0F); }
	static Color firebrick() { return Color(0.698F, 0.133F, 0.133F); }
	static Color floralwhite() { return Color(1.0F, 0.98F, 0.941F); }
	static Color forestgreen() { return Color(0.133F, 0.545F, 0.133F); }
	static Color gainsboro() { return Color(0.863F, 0.863F, 0.863F); }
	static Color ghostwhite() { return Color(0.973F, 0.973F, 1.0F); }
	static Color goldenrod() { return Color(0.855F, 0.647F, 0.125F); }
	static Color greenyellow() { return Color(0.678F, 1.0F, 0.184F); }
	static Color honeydew() { return Color(0.941F, 1.0F, 0.941F); }
	static Color hotpink() { return Color(1.0F, 0.412F, 0.706F); }
	static Color indianred() { return Color(0.804F, 0.361F, 0.361F); }
	static Color ivory() { return Color(1.0F, 1.0F, 0.941F); }
	static Color lavenderblush() { return Color(1.0F, 0.941F, 0.961F); }
	static Color lawngreen() { return Color(0.486F, 0.988F, 0.0F); }
	static Color lemonchiffon() { return Color(1.0F, 0.98F, 0.804F); }
	static Color lightblue() { return Color(0.678F, 0.847F, 0.902F); }
	static Color lightcoral() { return Color(0.941F, 0.502F, 0.502F); }
	static Color lightcyan() { return Color(0.878F, 1.0F, 1.0F); }
	static Color lightgoldenrodyellow() { return Color(0.98F, 0.98F, 0.824F); }
	static Color lightgreen() { return Color(0.565F, 0.933F, 0.565F); }
	static Color lightgrey() { return Color(0.827F, 0.827F, 0.827F); }
	static Color lightpink() { return Color(1.0F, 0.714F, 0.757F); }
	static Color lightsalmon() { return Color(1.0F, 0.627F, 0.478F); }
	static Color lightseagreen() { return Color(0.125F, 0.698F, 0.667F); }
	static Color lightskyblue() { return Color(0.529F, 0.808F, 0.98F); }
	static Color lightslategray() { return Color(0.467F, 0.533F, 0.6F); }
	static Color lightsteelblue() { return Color(0.69F, 0.769F, 0.871F); }
	static Color lightyellow() { return Color(1.0F, 1.0F, 0.878F); }
	static Color limegreen() { return Color(0.196F, 0.804F, 0.196F); }
	static Color linen() { return Color(0.98F, 0.941F, 0.902F); }
	static Color mediumaquamarine() { return Color(0.4F, 0.804F, 0.667F); }
	static Color mediumblue() { return Color(0.0F, 0.0F, 0.804F); }
	static Color mediumorchid() { return Color(0.729F, 0.333F, 0.827F); }
	static Color mediumpurple() { return Color(0.576F, 0.439F, 0.859F); }
	static Color mediumseagreen() { return Color(0.235F, 0.702F, 0.443F); }
	static Color mediumslateblue() { return Color(0.482F, 0.408F, 0.933F); }
	static Color mediumspringgreen() { return Color(0.0F, 0.98F, 0.604F); }
	static Color mediumturquoise() { return Color(0.282F, 0.82F, 0.8F); }
	static Color mediumvioletred() { return Color(0.78F, 0.082F, 0.522F); }
	static Color midnightblue() { return Color(0.098F, 0.098F, 0.439F); }
	static Color mintcream() { return Color(0.961F, 1.0F, 0.98F); }
	static Color mistyrose() { return Color(1.0F, 0.894F, 0.882F); }
	static Color navajowhite() { return Color(1.0F, 0.871F, 0.678F); }
	static Color oldlace() { return Color(0.992F, 0.961F, 0.902F); }
	static Color olivedrab() { return Color(0.42F, 0.557F, 0.137F); }
	static Color orangered() { return Color(1.0F, 0.271F, 0.0F); }
	static Color orchid() { return Color(0.855F, 0.439F, 0.839F); }
	static Color palegoldenrod() { return Color(0.933F, 0.91F, 0.667F); }
	static Color palegreen() { return Color(0.596F, 0.984F, 0.596F); }
	static Color paleturquoise() { return Color(0.686F, 0.933F, 0.933F); }
	static Color palevioletred() { return Color(0.859F, 0.439F, 0.576F); }
	static Color papayawhip() { return Color(1.0F, 0.937F, 0.835F); }
	static Color peachpuff() { return Color(1.0F, 0.855F, 0.725F); }
	static Color peru() { return Color(0.804F, 0.522F, 0.247F); }
	static Color plum() { return Color(0.867F, 0.627F, 0.867F); }
	static Color powderblue() { return Color(0.69F, 0.878F, 0.902F); }
	static Color rosybrown() { return Color(0.737F, 0.561F, 0.561F); }
	static Color royalblue() { return Color(0.255F, 0.412F, 0.882F); }
	static Color saddlebrown() { return Color(0.545F, 0.271F, 0.075F); }
	static Color sandybrown() { return Color(0.957F, 0.643F, 0.376F); }
	static Color seagreen() { return Color(0.18F, 0.545F, 0.341F); }
	static Color seashell() { return Color(1.0F, 0.961F, 0.933F); }
	static Color skyblue() { return Color(0.529F, 0.808F, 0.922F); }
	static Color slateblue() { return Color(0.416F, 0.353F, 0.804F); }
	static Color slategray() { return Color(0.439F, 0.502F, 0.565F); }
	static Color slategrey() { return Color(0.439F, 0.502F, 0.565F); }
	static Color snow() { return Color(1.0F, 0.98F, 0.98F); }
	static Color springgreen() { return Color(0.0F, 1.0F, 0.498F); }
	static Color steelblue() { return Color(0.275F, 0.51F, 0.706F); }
	static Color thistle() { return Color(0.847F, 0.749F, 0.847F); }
	static Color whitesmoke() { return Color(0.961F); }
	static Color yellowgreen() { return Color(0.604F, 0.804F, 0.196F); }
	// End generated colors

	std::string to_string() {
		std::stringstream ss;
		ss << std::to_string(r) << ", " << std::to_string(g) << ", " << std::to_string(b) << ", " << std::to_string(a);
		return ss.str();
	}

	void from_string(const std::string& str) {
		char skip;
		std::stringstream ss;
		ss << str;
		ss >> r >> skip >> g >> skip >> b >> skip >> a;
	}

	void print() {
		std::cout << "Color<" << r << ", " << g << ", " << b << ", " << a << ">";
	}

	static Color min(const Color& a, const Color& b) {
		return Color(std::min(a.r, b.r), std::min(a.g, b.g), std::min(a.b, b.b), std::min(a.a, b.a));
	}

	static Color max(const Color& a, const Color& b) {
		return Color(std::max(a.r, b.r), std::max(a.g, b.g), std::max(a.b, b.b), std::max(a.a, b.a));
	}

	static Color lerp(const Color& from, const Color& to, float t) {
#if defined(USE_SIMD)
		Color ret = to;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_ld(from.linear);
		vf32x4 vt = vf32x4_set(t, t, t, t);
		vf32x4_st(ret.linear, vf32x4_add(l, vf32x4_mul(vf32x4_sub(l, r), vt)));
		return ret;
#else
		return Color(from.r + (to.r - from.r) * t,
			from.g + (to.g - from.g) * t,
			from.b + (to.b - from.b) * t,
			from.a + (to.a - from.a) * t);
#endif
	}

	Color operator+(const Color& other) const {
#if defined(USE_SIMD)
		Color ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(ret.linear, vf32x4_add(l, r));
		return ret;
#else
		return Color(r + other.r, g + other.g, b + other.b, a + other.a);
#endif
	}

	Color operator-(const Color& other) const {
#if defined(USE_SIMD)
		Color ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(ret.linear, vf32x4_sub(l, r));
		return ret;
#else
		return Color(r - other.r, g - other.g, b - other.b, a - other.a);
#endif
	}

	Color operator*(float scalar) const {
#if defined(USE_SIMD)
		Color ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
		vf32x4_st(ret.linear, vf32x4_mul(l, r));
		return ret;
#else
		return Color(r * scalar, g * scalar, b * scalar, a * scalar);
#endif
	}

	Color operator*(const Color& other) const {
#if defined(USE_SIMD)
		Color ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(ret.linear, vf32x4_mul(l, r));
		return ret;
#else
		return Color(r * other.r, g * other.g, b * other.b, a * other.a);
#endif
	}

	Color operator/(float scalar) const {
#if defined(USE_SIMD)
		Color ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
		vf32x4_st(ret.linear, vf32x4_div(l, r));
		return ret;
#else
		return Color(r / scalar, g / scalar, b / scalar, a / scalar);
#endif
	}

	Color operator/(const Color& other) const {
#if defined(USE_SIMD)
		Color ret = *this;
		vf32x4 l = vf32x4_ld(ret.linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(ret.linear, vf32x4_div(l, r));
		return ret;
#else
		return Color(r / other.r, g / other.g, b / other.b, a / other.a);
#endif
	}

	void operator+=(const Color& other) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(linear, vf32x4_add(l, r));
#else
		r += other.r;
		g += other.g;
		b += other.b;
		a += other.a;
#endif
	}

	void operator-=(const Color& other) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(linear, vf32x4_sub(l, r));
#else
		r -= other.r;
		g -= other.g;
		b -= other.b;
		a -= other.a;
#endif
	}

	void operator*=(float scalar) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
		vf32x4_st(linear, vf32x4_mul(l, r));
#else
		r *= scalar;
		g *= scalar;
		b *= scalar;
		a *= scalar;
#endif
	}

	void operator*=(const Color& other) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(linear, vf32x4_add(l, r));
#else
		r *= other.r;
		g *= other.g;
		b *= other.b;
		a *= other.a;
#endif
	}

	void operator/=(float scalar) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_set(scalar, scalar, scalar, scalar);
		vf32x4_st(linear, vf32x4_div(l, r));
#else
		r /= scalar;
		g /= scalar;
		b /= scalar;
		a /= scalar;
#endif
	}

	void operator/=(const Color& other) {
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(linear, vf32x4_div(l, r));
#else
		r /= other.r;
		g /= other.g;
		b /= other.b;
		a /= other.a;
#endif
	}

	bool operator==(const Color& other) const {
		float epsilon = std::numeric_limits<float>::epsilon();
		float d[4];
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(d, vf32x4_sub(l, r));
#else
		d[0] = r - other.r;
		d[1] = g - other.g;
		d[2] = b - other.b;
		d[3] = a - other.a;
#endif
		// TODO:  Could probably use the abs of the sum rather than each
		return std::abs(d[0]) < epsilon
			&& std::abs(d[1]) < epsilon
			&& std::abs(d[2]) < epsilon
			&& std::abs(d[3]) < epsilon;
	}

	bool operator!=(const Color& other) const {
		float epsilon = std::numeric_limits<float>::epsilon();
		float d[4];
#if defined(USE_SIMD)
		vf32x4 l = vf32x4_ld(linear);
		vf32x4 r = vf32x4_ld(other.linear);
		vf32x4_st(d, vf32x4_sub(l, r));
#else
		d[0] = r - other.r;
		d[1] = g - other.g;
		d[2] = b - other.b;
		d[3] = a - other.a;
#endif
		// TODO:  Could probably use the abs of the sum rather than each
		return std::abs(d[0]) > epsilon
			&& std::abs(d[1]) > epsilon
			&& std::abs(d[2]) > epsilon
			&& std::abs(d[3]) > epsilon;
	}
};

#endif
