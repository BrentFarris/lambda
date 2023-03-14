#ifdef ANDROID
	#if defined(__arm__) || defined(__aarch64__)
		#define SIMD_NEON
		#include <arm_neon.h>
	#endif
#elif (defined(__x86_64__) || defined(__amd64) || defined(_WIN64))
	#define SIMD
	#include <immintrin.h>
#endif

#ifdef SIMD
typedef __m128 vf32x4;
typedef __m256 vf32x8;
typedef __m128i vi32x4;
typedef __m256d vf64x4;

#define vf32x4_set(f0, f1, f2, f3) _mm_setr_ps((f0), (f1), (f2), (f3))
#define vf32x4_ld(f32Arr4) _mm_load_ps(f32Arr4)
#define vf32x4_st(f32Arr4, f32x4) _mm_store_ps((f32Arr4), (f32x4))
#define vf32x4_add(L_f32x4, R_f32x4) _mm_add_ps((L_f32x4), (R_f32x4))
#define vf32x4_sub(L_f32x4, R_f32x4) _mm_sub_ps((L_f32x4), (R_f32x4))
#define vf32x4_mul(L_f32x4, R_f32x4) _mm_mul_ps((L_f32x4), (R_f32x4))
#define vf32x4_div(L_f32x4, R_f32x4) _mm_div_ps((L_f32x4), (R_f32x4))
#ifdef __linux__
#define vf32x4_sum(f32x4) (f32x4)[0] + (f32x4)[1] + (f32x4)[2] + (f32x4)[3]
#else
#define vf32x4_sum(f32x4) (f32x4).m128_f32[0] + (f32x4).m128_f32[1] + (f32x4).m128_f32[2] + (f32x4).m128_f32[3]
#endif

#define vf64x4_set(f0, f1, f2, f3) _mm256_setr_pd((f0), (f1), (f2), (f3))
#define vf64x4_ld(f64Arr4) _mm256_load_pd(f64Arr4)
#define vf64x4_st(f64Arr4, f64x4) _mm256_store_pd((f64Arr4), (f64x4))
#define vf64x4_add(L_f64x4, R_f64x4) _mm256_add_pd((L_f64x4), (R_f64x4))
#define vf64x4_sub(L_f64x4, R_f64x4) _mm256_sub_pd((L_f64x4), (R_f64x4))
#define vf64x4_mul(L_f64x4, R_f64x4) _mm256_mul_pd((L_f64x4), (R_f64x4))
#define vf64x4_div(L_f64x4, R_f64x4) _mm256_div_pd((L_f64x4), (R_f64x4))
#ifdef __linux__
#define vf64x4_sum(f64x4) (f64x4)[0] + (f64x4)[1] + (f64x4)[2] + (f64x4)[3]
#else
#define vf64x4_sum(f64x4) (f64x4).m256d_f64[0] + (f64x4).m256d_f64[1] + (f64x4).m256d_f64[2] + (f64x4).m256d_f64[3]
#endif

#define vf32x8_set(f0, f1, f2, f3, f4, f5, f6, f7) _mm256_setr_ps((f0), (f1), (f2), (f3), (f4), (f5), (f6), (f7))
#define vf32x8_ld(f32Arr8) _mm256_load_ps(f32Arr8)
#define vf32x8_st(f32Arr8, f32x8) _mm256_store_ps((f32Arr8), (f32x8))
#define vf32x8_add(outName, L_f32x8, R_f32x8) vf32x8 (outName) = _mm256_add_ps((L_f32x8), (R_f32x8))
#define vf32x8_sub(outName, L_f32x8, R_f32x8) vf32x8 (outName) = _mm256_sub_ps((L_f32x8), (R_f32x8))
#define vf32x8_mul(outName, L_f32x8, R_f32x8) vf32x8 (outName) = _mm256_mul_ps((L_f32x8), (R_f32x8))
#define vf32x8_div(outName, L_f32x8, R_f32x8) vf32x8 (outName) = _mm256_div_ps((L_f32x8), (R_f32x8))

#define vf64x8_st(f64Arr8, f64x8) _mm256_store_pd((f64Arr8), (f64x8))
#define vf64x4n_add(outName, L_f64x8, R_f64x8) vf64x4 (outName) = _mm256_add_pd((L_f64x8), (R_f64x8))
#define vf64x4n_sub(outName, L_f64x8, R_f64x8) vf64x4 (outName) = _mm256_sub_pd((L_f64x8), (R_f64x8))
#define vf64x4n_mul(outName, L_f64x8, R_f64x8) vf64x4 (outName) = _mm256_mul_pd((L_f64x8), (R_f64x8))
#define vf64x4n_div(outName, L_f64x8, R_f64x8) vf64x4 (outName) = _mm256_div_pd((L_f64x8), (R_f64x8))

#define vi32x4_set(i0, i1, i2, i3) _mm_setr_epi32((i0), (i1), (i2), (i3))
#define vi32x4_ld(i32Arr4) _mm_load_epi32(i32Arr4)
#define vi32x4_st(i32Arr4, i32x4) _mm_store_epi32((i32Arr4), (i32x4))
#define vi32x4_add(L_i32x4, R_i32x4) _mm_add_epi32((L_i32x4), (R_i32x4))
#define vi32x4_sub(L_i32x4, R_i32x4) _mm_sub_epi32((L_i32x4), (R_i32x4))
#define vi32x4_mul(L_i32x4, R_i32x4) _mm_mul_epi32((L_i32x4), (R_i32x4))
#define vi32x4_div(L_i32x4, R_i32x4) _mm_div_epi32((L_i32x4), (R_i32x4))
#elif defined(SIMD_NEON)
typedef float32x4_t vf32x4;
typedef float32x4x2_t vf32x8;
typedef int32x4_t vi32x4;

// TODO:  NEON version of __m128d and __m256d

#define vf32x4_set(f0, f1, f2, f3) { (f0), (f1), (f2), (f3) }
#define vf32x4_ld(floatArr4) vld1q_f32(floatArr4)
#define vf32x4_st(floatArr4, f32x4) vst1q_f32((floatArr4), (f32x4))
#define vf32x4_add(L_f32x4, R_f32x4) vaddq_f32((L_f32x4), (R_f32x4))
#define vf32x4_sub(L_f32x4, R_f32x4) vsubq_f32((L_f32x4), (R_f32x4))
#define vf32x4_mul(L_f32x4, R_f32x4) vmulq_f32((L_f32x4), (R_f32x4))
#define vf32x4_sum(f32x4) f32x4[0] + f32x4[1] + f32x4[2] + f32x4[3];

#define vf32x8_set(f0, f1, f2, f3, f4, f5, f6, f7) { { [0] = (f0), (f1), (f2), (f3), [1] = (f4), (f5), (f6), (f7) } }
#define vf32x8_ld(f32Arr8) vld2q_f32(f32Arr8)
#define vf32x8_st(f32Arr8, f32x8) do { vst1q_f32((f32Arr8), (f32x8.val[0])); vst1q_f32((f32Arr8) + 4, (f32x8.val[1])); } while(0)
#define vf32x8_add(outName, L_f32x8, R_f32x8) vf32x8 (outName); do {	\
	(outName).val[0] = vaddq_f32((L_f32x8.val[0]), (R_f32x8.val[0]));	\
	(outName).val[1] = vaddq_f32((L_f32x8.val[1]), (R_f32x8.val[1]));	\
} while(false)
#define vf32x8_sub(outName, L_f32x8, R_f32x8)  vf32x8 (outName); do {	\
	(outName).val[0] = vsubq_f32((L_f32x8.val[0]), (R_f32x8.val[0]));	\
	(outName).val[1] = vsubq_f32((L_f32x8.val[1]), (R_f32x8.val[1]));	\
} while(false)
#define vf32x8_mul(outName, L_f32x8, R_f32x8)  vf32x8 (outName); do {	\
	(outName).val[0] = vmulq_f32((L_f32x8.val[0]), (R_f32x8.val[0]));	\
	(outName).val[1] = vmulq_f32((L_f32x8.val[1]), (R_f32x8.val[1]));	\
} while(false)

#define vi32x4_set(i0, i1, i2, i3) { (i0), (i1), (i2), (i3) }
#define vi32x4_ld(i32Arr4) vld1q_s32(i32Arr4)
#define vi32x4_st(i32Arr4, i32x4) vst1q_s32((i32Arr4), (i32x4))
#define vi32x4_add(L_i32x4, R_i32x4) vaddq_s32((L_i32x4), (R_i32x4))
#define vi32x4_sub(L_i32x4, R_i32x4) vsubq_s32((L_i32x4), (R_i32x4))
#define vi32x4_mul(L_i32x4, R_i32x4) vmulq_s32((L_i32x4), (R_i32x4))
#endif

#if defined(SIMD) || defined(SIMD_NEON)
#define USE_SIMD
#endif