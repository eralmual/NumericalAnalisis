/*
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Pablo Alvarado
 * @Date:   28.12.2017
 */

#ifndef ANPI_MATRIX_ARITHMETIC_HPP
#define ANPI_MATRIX_ARITHMETIC_HPP

#include "Intrinsics.hpp"
#include <type_traits>
#include "Matrix.hpp"
#include "Exception.hpp"
#include "IntrinsicsMethods.hpp"
#include <functional>
#include <iostream>

namespace anpi {
    namespace fallback {
        /*
         * Sum
         */

        // Fallback implementation

        // In-copy implementation c=a+b
        template<typename T, class Alloc>
        inline void add(const Matrix <T, Alloc> &a,
                        const Matrix <T, Alloc> &b,
                        Matrix <T, Alloc> &c) {

            assert((a.rows() == b.rows()) &&
                   (a.cols() == b.cols()));

            const size_t tentries = a.rows() * a.dcols();
            c.allocate(a.rows(), a.cols());

            T *here = c.data();
            T *const end = here + tentries;
            const T *aptr = a.data();
            const T *bptr = b.data();

            for (; here != end;) {
                *here++ = *aptr++ + *bptr++;
            }
        }

        // In-place implementation a = a+b
        template<typename T, class Alloc>
        inline void add(Matrix <T, Alloc> &a,
                        const Matrix <T, Alloc> &b) {

            assert((a.rows() == b.rows()) &&
                   (a.cols() == b.cols()));

            const size_t tentries = a.rows() * a.dcols();

            T *here = a.data();
            T *const end = here + tentries;

            const T *bptr = b.data();

            for (; here != end;) {
                *here++ += *bptr++;
            }
        }


        /*
         * Subtraction
         */

        // Fall back implementations

        // In-copy implementation c=a-b
        template<typename T, class Alloc>
        inline void subtract(const Matrix <T, Alloc> &a,
                             const Matrix <T, Alloc> &b,
                             Matrix <T, Alloc> &c) {

            assert((a.rows() == b.rows()) &&
                   (a.cols() == b.cols()));

            const size_t tentries = a.rows() * a.dcols();
            c.allocate(a.rows(), a.cols());

            T *here = c.data();
            T *const end = here + tentries;
            const T *aptr = a.data();
            const T *bptr = b.data();

            for (; here != end;) {
                *here++ = *aptr++ - *bptr++;
            }
        }

        // In-place implementation a = a-b
        template<typename T, class Alloc>
        inline void subtract(Matrix <T, Alloc> &a,
                             const Matrix <T, Alloc> &b) {

            assert((a.rows() == b.rows()) &&
                   (a.cols() == b.cols()));

            const size_t tentries = a.rows() * a.dcols();

            T *here = a.data();
            T *const end = here + tentries;

            const T *bptr = b.data();

            for (; here != end;) {
                *here++ -= *bptr++;
            }
        }


        // In-place implementation a = a*b
        template<typename T, class Alloc>
        inline void multiplicate(Matrix <T, Alloc> &a,
                                 const Matrix <T, Alloc> &b) {

            //Se verifica que las matrices sean compatibles para
            //realizar la multiplicacion
            if (a.cols() != b.rows()) {
                throw anpi::Exception("Las dimensiones de las matrices no son compatibles");
            } else {
                //Se crea la matriz resultante y se inicializa en 0
                Matrix<T> mult;
                mult.allocate(a.rows(), b.cols());
                mult.fill(T(0));
                // Se hace la multiplicacion y se almacena en la matriz resultante
                for (unsigned int i = 0; i < a.rows(); ++i) {

                    for (unsigned int j = 0; j < b.cols(); ++j) {

                        for (unsigned int k = 0; k < a.cols(); ++k) {

                            mult(i, j) += a(i, k) * b(k, j);
                        }
                    }
                }
                a = mult;
            }
        }


        // In-copy implementation c=a*b
        template<typename T, class Alloc>
        inline void multiplicate(const Matrix <T, Alloc> &a,
                                 const Matrix <T, Alloc> &b,
                                 Matrix <T, Alloc> &mult) {

            //Se verifica que las matrices sean compatibles para
            //realizar la multiplicacion
            if (a.cols() != b.rows()) {
                throw anpi::Exception("Las dimensiones de las matrices no son compatibles");
            } else {
                //Se crea la matriz resultante y se inicializa en 0
                mult.fill(T(0));
                // Se hace la multiplicacion y se almacena en la matriz resultante
                for (unsigned int i = 0; i < a.rows(); ++i) {

                    for (unsigned int j = 0; j < b.cols(); ++j) {

                        for (unsigned int k = 0; k < a.cols(); ++k) {

                            mult(i, j) += a(i, k) * b(k, j);
                        }
                    }
                }
            }
        }

        // In-copy implementation C = A/b
        template<typename T, typename U, class Alloc>
        inline void divide(const Matrix <T, Alloc> &a,
                           Matrix <T, Alloc> &div,
                           const U factor) {

            if (factor == U(0)) {
                throw anpi::Exception("Se intentó dividir todos los elementos de la matriz por cero");
            } else {

                T fct = T(factor);
                //Se itera dividiendo entre cada elemnto de la matriz
                for (size_t i = 0; i < a.rows(); ++i) {

                    for (size_t j = 0; j < a.cols(); ++j) {
                        div(i, j) = a(i, j) / fct;
                    }

                }

            }
        }

    } // namespace fallback


    namespace simd {


        // On-copy implementation c=a+b
        template<typename T, class Alloc, typename regType>
        inline void addSIMD(const Matrix <T, Alloc> &a,
                            const Matrix <T, Alloc> &b,
                            Matrix <T, Alloc> &c) {

            // This method is instantiated with unaligned allocators.  We
            // allow the instantiation although externally this is never
            // called unaligned
            static_assert(!extract_alignment<Alloc>::aligned ||
                          (extract_alignment<Alloc>::value >= sizeof(regType)),
                          "Insufficient alignment for the registers used");

            const size_t tentries = a.rows() * a.dcols();
            c.allocate(a.rows(), a.cols());

            regType *here = reinterpret_cast<regType *>(c.data());
            const size_t blocks = (tentries * sizeof(T) + (sizeof(regType) - 1)) /
                                  sizeof(regType);
            regType *const end = here + blocks;
            const regType *aptr = reinterpret_cast<const regType *>(a.data());
            const regType *bptr = reinterpret_cast<const regType *>(b.data());

            for (; here != end;) {
                *here++ = mm_add<T>(*aptr++, *bptr++);
            }

        }

        // On-copy implementation c=a+b for SIMD-capable types
        template<typename T,
                class Alloc,
                typename std::enable_if<is_simd_type<T>::value, int>::type= 0>
        inline void add(const Matrix <T, Alloc> &a,
                        const Matrix <T, Alloc> &b,
                        Matrix <T, Alloc> &c) {

            assert((a.rows() == b.rows()) &&
                   (a.cols() == b.cols()));


            if (is_aligned_alloc<Alloc>::value) {
#ifdef __AVX512F__
                addSIMD<T,Alloc,typename avx512_traits<T>::reg_type>(a,b,c);
#elif  __AVX__
                addSIMD<T, Alloc, typename avx_traits<T>::reg_type>(a, b, c);
#elif  __SSE2__
                addSIMD<T,Alloc,typename sse2_traits<T>::reg_type>(a,b,c);
#else
        ::anpi::fallback::add(a,b,c);
#endif
            } else { // allocator seems to be unaligned
                ::anpi::fallback::add(a, b, c);
            }
        }


        // On-copy implementation c=a-b
        template<typename T, class Alloc, typename regType>
        inline void subSIMD(const Matrix <T, Alloc> &a,
                            const Matrix <T, Alloc> &b,
                            Matrix <T, Alloc> &c) {

            // This method is instantiated with unaligned allocators.  We
            // allow the instantiation although externally this is never
            // called unaligned
            static_assert(!extract_alignment<Alloc>::aligned ||
                          (extract_alignment<Alloc>::value >= sizeof(regType)),
                          "Insufficient alignment for the registers used");

            const size_t tentries = a.rows() * a.dcols();
            c.allocate(a.rows(), a.cols());

            regType *here = reinterpret_cast<regType *>(c.data());
            const size_t blocks = (tentries * sizeof(T) + (sizeof(regType) - 1)) /
                                  sizeof(regType);
            regType *const end = here + blocks;
            const regType *aptr = reinterpret_cast<const regType *>(a.data());
            const regType *bptr = reinterpret_cast<const regType *>(b.data());


            for (; here != end;) {
                *here++ = mm_sub<T>(*aptr++, *bptr++);
            }

        }


        // On-copy implementation c=a-b for SIMD-capable types
        template<typename T,
                class Alloc,
                typename std::enable_if<is_simd_type<T>::value, int>::type= 0>
        inline void subtract(const Matrix <T, Alloc> &a,
                             const Matrix <T, Alloc> &b,
                             Matrix <T, Alloc> &c) {

            assert((a.rows() == b.rows()) &&
                   (a.cols() == b.cols()));


            if (is_aligned_alloc<Alloc>::value) {
#ifdef  __AVX__
                subSIMD<T, Alloc, typename avx_traits<T>::reg_type>(a, b, c);

#endif
            } else { // allocator seems to be unaligned
                ::anpi::fallback::subtract(a, b, c);
            }
        }

        // On-copy implementation C = A/b
        template<typename T, class Alloc, typename regType>
        inline void divSIMD(const Matrix <T, Alloc> &a,
                            const Matrix <T, Alloc> &factor,
                            Matrix <T, Alloc> &c) {

            // This method is instantiated with unaligned allocators.  We
            // allow the instantiation although externally this is never
            // called unaligned
            static_assert(!extract_alignment<Alloc>::aligned ||
                          (extract_alignment<Alloc>::value >= sizeof(regType)),
                          "Insufficient alignment for the registers used");

            const size_t tentries = a.rows() * a.dcols();
            c.allocate(a.rows(), a.cols());

            regType *here = reinterpret_cast<regType *>(c.data());
            const size_t blocks = (tentries * sizeof(T) + (sizeof(regType) - 1)) /
                                  sizeof(regType);
            regType *const end = here + blocks;
            const regType *aptr = reinterpret_cast<const regType *>(a.data());
            const regType *bptr = reinterpret_cast<const regType *>(factor.data());

            for (; here != end;) {
                *here++ = mm_div<T>(*aptr++, *bptr++);
            }

        }


        // On-copy implementation c=a-b for SIMD-capable types
        template<typename T,
                class Alloc,
                typename std::enable_if<is_simd_type<T>::value, int>::type= 0>
        inline void divide(const Matrix <T, Alloc> &a,
                           const Matrix <T, Alloc> &b,//factor
                           Matrix <T, Alloc> &c) {

            assert((a.rows() == b.rows()) &&
                   (a.cols() == b.cols()));


            if (is_aligned_alloc<Alloc>::value) {
#ifdef __AVX__
                divSIMD<T, Alloc, typename avx_traits<T>::reg_type>(a, b, c);
#endif
            } else { // allocator seems to be unaligned
                ::anpi::fallback::divide(a, c, b(0, 0));
            }
        }

        // Non-SIMD types such as complex
        template<typename T,
                class Alloc,
                typename std::enable_if<!is_simd_type<T>::value, int>::type = 0>
        inline void add(const Matrix <T, Alloc> &a,
                        const Matrix <T, Alloc> &b,
                        Matrix <T, Alloc> &c) {

            ::anpi::fallback::add(a, b, c);
        }

        // In-place implementation a = a+b
        template<typename T, class Alloc>
        inline void add(Matrix <T, Alloc> &a,
                        const Matrix <T, Alloc> &b) {

            add(a, b, a);
        }


        /*
         * Subtraction
         */

        // Fall back implementations

        // In-copy implementation c=a-b
        template<typename T,
                class Alloc,
                typename std::enable_if<!is_simd_type<T>::value, int>::type = 0>
        inline void subtract(const Matrix <T, Alloc> &a,
                             const Matrix <T, Alloc> &b,
                             Matrix <T, Alloc> &c) {
            ::anpi::fallback::subtract(a, b, c);
        }

        // In-place implementation a = a-b
        template<typename T, class Alloc>
        inline void subtract(Matrix <T, Alloc> &a,
                             const Matrix <T, Alloc> &b) {

            subtract(a, b, a);
        }



        /*
         * Division
         */

        // Fall back implementations

        // In-copy implementation c=a-b
        template<typename T,
                class Alloc,
                typename std::enable_if<!is_simd_type<T>::value, int>::type = 0>
        inline void divide(const Matrix <T, Alloc> &a,
                           const Matrix <T, Alloc> &factor,
                           Matrix <T, Alloc> &c) {
            ::anpi::fallback::divide(a, c, factor(0, 0));
        }


        // In-place implementation  A = A/b
        template<typename T, class Alloc>
        inline void divide(Matrix <T, Alloc> &a,
                           const Matrix <T, Alloc> &factor) {

            divide(a, factor, a);
        }


    } // namespace simd


    // The arithmetic implementation (aimpl) namespace
    // dispatches to the corresponding methods
#ifdef ANPI_ENABLE_SIMD
    namespace aimpl = simd;
#else
    namespace aimpl=fallback;
#endif

} // namespace anpi

#endif
