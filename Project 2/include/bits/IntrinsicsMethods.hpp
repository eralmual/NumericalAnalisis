//
// Created by erick on 07/10/18.
//

#ifndef PROYECTO2_INTRINSICSMETHODS_H
#define PROYECTO2_INTRINSICSMETHODS_H

#include "Intrinsics.hpp"

//-------------------------------------------- AVX only --------------------------------------------

/**
 * Load method for alligned registers
 * @tparam T        Datatype
 * @tparam regType  Register datatype
 * @return Pointer to a register loaded with the pointer data
 */
template<typename T, class regType>
regType mm_loadRegister(const T *);
#ifdef  __AVX__

    template<>
    inline __m256d __attribute__((__always_inline__))
    mm_loadRegister<double>(const double *a) {
        return _mm256_load_pd(a);
    }

    template<>
    inline __m256 __attribute__((__always_inline__))
    mm_loadRegister<float>(const float *a) {
        return _mm256_load_ps(a);
    }

#endif


/**
 * Load method for unalligned registers
 * @tparam T        Datatype
 * @tparam regType  Register datatype
 * @return Pointer to a register loaded with the pointer data
 */
template<typename T, class regType>
regType mm_loadRegisteru(const T *);
#ifdef  __AVX__

    template<>
    inline __m256d __attribute__((__always_inline__))
    mm_loadRegisteru<double>(const double *a) {
        return _mm256_loadu_pd(a);
    }

    template<>
    inline __m256 __attribute__((__always_inline__))
    mm_loadRegisteru<float>(const float *a) {
        return _mm256_loadu_ps(a);
    }

#endif


/**
 * Implementation of substraction
 * @tparam T        Datatype
 * @tparam regType  Register datatype
 * @param a         Register to be substracted
 * @param b         Register with the values to substract
 * @return          Register with the result of the operation
 */
template<typename T, class regType>
regType mm_sub(regType, regType);

#ifdef __AVX__

    template<>
    inline __m256d __attribute__((__always_inline__))
    mm_sub<double>(__m256d a, __m256d b) {
        return _mm256_sub_pd(a, b);
    }

    template<>
    inline __m256 __attribute__((__always_inline__))
    mm_sub<float>(__m256 a, __m256 b) {
        return _mm256_sub_ps(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<uint64_t>(__m256i a, __m256i b) {
        return _mm256_sub_epi64(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<int64_t>(__m256i a, __m256i b) {
        return _mm256_sub_epi64(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<uint32_t>(__m256i a, __m256i b) {
        return _mm256_sub_epi32(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<int32_t>(__m256i a, __m256i b) {
        return _mm256_sub_epi32(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<uint16_t>(__m256i a, __m256i b) {
        return _mm256_sub_epi16(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<int16_t>(__m256i a, __m256i b) {
        return _mm256_sub_epi16(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<uint8_t>(__m256i a, __m256i b) {
        return _mm256_sub_epi8(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_sub<int8_t>(__m256i a, __m256i b) {
        return _mm256_sub_epi8(a, b);
    }
#endif






/**
 * Implementation of addition
 * @tparam T        Datatype
 * @tparam regType  Register datatype
 * @param a         Register to be added
 * @param b         Register with the values to add
 * @return          Register with the result of the operation
 */
template<typename T, class regType>
regType mm_add(regType, regType);


#ifdef __AVX__

    template<>
    inline __m256d __attribute__((__always_inline__))
    mm_add<double>(__m256d a, __m256d b) {
        return _mm256_add_pd(a, b);
    }

    template<>
    inline __m256 __attribute__((__always_inline__))
    mm_add<float>(__m256 a, __m256 b) {
        return _mm256_add_ps(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint64_t>(__m256i a, __m256i b) {
        return _mm256_add_epi64(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int64_t>(__m256i a, __m256i b) {
        return _mm256_add_epi64(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint32_t>(__m256i a, __m256i b) {
        return _mm256_add_epi32(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int32_t>(__m256i a, __m256i b) {
        return _mm256_add_epi32(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint16_t>(__m256i a, __m256i b) {
        return _mm256_add_epi16(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int16_t>(__m256i a, __m256i b) {
        return _mm256_add_epi16(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<uint8_t>(__m256i a, __m256i b) {
        return _mm256_add_epi8(a, b);
    }

    template<>
    inline __m256i __attribute__((__always_inline__))
    mm_add<int8_t>(__m256i a, __m256i b) {
        return _mm256_add_epi8(a, b);
    }

#endif


/**
 * Implementation of division
 * @tparam T        Datatype
 * @tparam regType  Register datatype
 * @param a         Register to be divided
 * @param b         Register with the values to divide
 * @return          Register with the result of the operation
 */
template<typename T, class regType>
regType mm_div(regType, regType);

#ifdef  __AVX__

    template<>
    inline __m256d __attribute__((__always_inline__))
    mm_div<double>(__m256d a, __m256d b) {
        return _mm256_div_pd(a, b);
    }

    template<>
    inline __m256 __attribute__((__always_inline__))
    mm_div<float>(__m256 a, __m256 b) {
        return _mm256_div_ps(a, b);
    }

#endif



/**
 * Implementation of multiplication
 * @tparam T        Datatype
 * @tparam regType  Register datatype
 * @param a         Register to be added
 * @param b         Register with the values to add
 * @return          Register with the result of the operation
 */
template<typename T, class regType>
regType mm_mult(regType, regType);


#ifdef __AVX__

template<>
inline __m256d __attribute__((__always_inline__))
mm_mult<double>(__m256d a, __m256d b) {
    return _mm256_mul_pd(a, b);
}

template<>
inline __m256 __attribute__((__always_inline__))
mm_mult<float>(__m256 a, __m256 b) {
    return _mm256_mul_ps(a, b);
}

template<>
inline __m256i __attribute__((__always_inline__))
mm_mult<uint64_t>(__m256i a, __m256i b) {
    return _mm256_mullo_epi64(a, b);
}

template<>
inline __m256i __attribute__((__always_inline__))
mm_mult<int64_t>(__m256i a, __m256i b) {
    return _mm256_mullo_epi64(a, b);
}

template<>
inline __m256i __attribute__((__always_inline__))
mm_mult<uint32_t>(__m256i a, __m256i b) {
    return _mm256_mul_epi32(a, b);
}

template<>
inline __m256i __attribute__((__always_inline__))
mm_mult<int32_t>(__m256i a, __m256i b) {
    return _mm256_mul_epi32(a, b);
}

template<>
inline __m256i __attribute__((__always_inline__))
mm_mult<uint16_t>(__m256i a, __m256i b) {
    return _mm256_mullo_epi16(a, b);
}

template<>
inline __m256i __attribute__((__always_inline__))
mm_mult<int16_t>(__m256i a, __m256i b) {
    return _mm256_mullo_epi16(a, b);
}

#endif







#endif //PROYECTO2_INTRINSICSMETHODS_H
