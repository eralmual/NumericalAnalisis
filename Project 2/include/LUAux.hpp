/*
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Erick Muñoz
 * @Date:   06.10.2018
 */

#include <cmath>
#include <limits>
#include <algorithm>


#include "Matrix.hpp"


#ifndef PROYECTO2_LUAUX_H
#define PROYECTO2_LUAUX_H


namespace anpi {


    namespace fallback {

        /**
        * this method swaps an entire row of the matrix with another
        * @tparam T     Typename
        * @param A      The matrix where ee store the rows
        * @param row1   The row we'll pivot
        * @param row2   The other row we'll pivot
        */
        template<typename T, class Alloc>
        void swapRows(Matrix <T, Alloc> &A, size_t row1, size_t row2) {
            for (size_t i = 0; i < A.cols(); ++i) {
                T temp = A[row1][i];
                A[row1][i] = A[row2][i];
                A[row2][i] = temp;
            }
        }


        template<typename T, class Alloc>
        inline void luDoolittle(Matrix <T, Alloc> &LU,
                                std::vector<size_t> &permut) {

            size_t rows = LU.rows();

            for (size_t k = 0; k < rows - 1; ++k) {
                pivot(LU, k, permut);
                for (size_t i = k + 1; i < rows; ++i) {
                    if (std::abs(LU[k][k]) < std::numeric_limits<T>::epsilon()) {
                        throw anpi::Exception("error, division by 0 in pivot -> singular matrix cannot perform LU");
                    }

                    // Calc the elimination factor
                    const T eliminationFactor = LU[i][k] / LU[k][k];

                    // Store the factor in the L matrix
                    LU[i][k] = eliminationFactor;

                    // Substract the factor from the row
                    for (size_t j = k + 1; j < rows; ++j) {
                        LU[i][j] -= eliminationFactor * LU[k][j];
                    }
                }
            }
        }

    } //namespace fallback


    namespace simd {
        template<typename T, class Alloc, typename regType>
        inline void swapRowsSIMD(Matrix <T, Alloc> &A,
                                 size_t row1,
                                 size_t row2) {

            // This method is instantiated with unaligned allocators.  We
            // allow the instantiation although externally this is never
            // called unaligned
            static_assert(!extract_alignment<Alloc>::aligned ||
                          (extract_alignment<Alloc>::value >= sizeof(regType)),
                          "Insufficient alignment for the registers used");



            //We calc the lenght of the rows in memory
            const size_t tentries = A.dcols(); //Applied to only one row
            const size_t blocks = (tentries * sizeof(T) + (sizeof(regType) - 1)) /
                                  sizeof(regType);

            // We calc how many numbers we can load per block so later we can
            // move the pointers to avoid loading the same numbers multiple times
            const size_t sizeOfBlocks = sizeof(regType) / sizeof(T);

            //We get the pointers to the rows we are going to switch
            // Pointer to row1
            // We need a deep copy of this pointer because we are going to overwrite its contents
            T a[tentries];
            T *row1ptr = A.data() + row1 * A.dcols();
            std::copy(row1ptr, row1ptr + A.dcols(), a);
            row1ptr = a;

            //Pointer to row2
            T *row2ptr = A.data() + row2 * A.dcols();

            //We load the rows we need into registers
            regType *row1Reg = reinterpret_cast<regType *>(A.data() + row1 * A.dcols());
            regType *row2Reg = reinterpret_cast<regType *>(A.data() + row2 * A.dcols());

            //Calc the end of the row
            regType *const endRow1 = row1Reg + blocks;


            for (; row1Reg != endRow1;) {
                // row1 = row2
                *row1Reg++ = mm_loadRegister<T, regType>(row2ptr);

                // row2 = row1deepCopy
                *row2Reg++ = mm_loadRegisteru<T, regType>(row1ptr);

                // Move the registers to the next block
                row1ptr += sizeOfBlocks;
                row2ptr += sizeOfBlocks;
            }
        }

        // On-copy implementation c=a-b for SIMD-capable types
        template<typename T,
                class Alloc,
                typename std::enable_if<is_simd_type<T>::value, int>::type = 0>
        inline void swapRows(Matrix <T, Alloc> &A,
                             size_t row1,
                             size_t row2) {


            if (is_aligned_alloc<Alloc>::value && (std::is_same<T, float>::value || std::is_same<T, double>::value)) {
#ifdef  __AVX__
                swapRowsSIMD<T, Alloc, typename avx_traits<T>::reg_type>(A, row1, row2);
#else
                ::anpi::fallback::swapRows(A, row1, row2);
#endif
            } else { // allocator seems to be unaligned
                ::anpi::fallback::swapRows(A, row1, row2);
            }
        }

        // Non-SIMD types such as complex
        template<typename T,
                class Alloc,
                typename std::enable_if<!is_simd_type<T>::value, int>::type = 0>
        inline void swapRows(Matrix <T, Alloc> &A,
                             size_t row1,
                             size_t row2) {

            ::anpi::fallback::swapRows(A, row1, row2);
        }


        /**
         * Applies gaussian elimination to a LU matrix and stores the L matrix in situ
         *
         * @tparam T        Typename
         * @tparam Alloc    Allocator
         * @tparam regType  Registertype
         * @param LU        The matrix that we wiĺl work on
         * @param permut    The permutation vector for the matrix
         */
        template<typename T, class Alloc, typename regType>
        inline void luDoolittleSIMD(Matrix <T, Alloc> &LU,
                                    std::vector<size_t> &permut) {

            size_t rows = LU.rows();// Remember that in LU rows = cols

            //We get the total number of elements in the matrix
            const size_t tentries = LU.rows() * LU.dcols();
            //We calc the lenght of the matrix in in terms of blocks
            const size_t blocks = (tentries * sizeof(T) + (sizeof(regType) - 1)) /
                                  sizeof(regType);

            // We create block-counting related variables
            const size_t blocksByRow = blocks / rows;
            const size_t columnsPerBlock = LU.dcols() / blocksByRow;


            // The pointer to the begining of the register that contains the matrix
            regType *luReg = reinterpret_cast<regType *>(LU.data());


            // We create auxiliar objects
            const size_t sizeOfDatatype = sizeof(T);

            // Lets create the vector where we will store the elimination factors
            std::vector<T, Alloc> eliminationVector(sizeof(regType) / sizeOfDatatype);
            regType *eliminationRegister = reinterpret_cast<regType *>(eliminationVector.data());

            // These are the registers used to go through the matrix
            regType *luRegRowBeforePtr = nullptr;
            regType *luRegPtr = nullptr;
            regType *luEndRegPtr = nullptr;



            // We perform gauss elimination and store the L matrix in situ
            for (size_t i = 0; i < rows - 1; ++i) {

                // We pivot rows if necessary
                pivot(LU, i, permut);

                //printf("%zu",i);

                // Gets the number of blocks we need to move in order to
                // avoid changing the value of the L matrix
                // the aux is needed because the ceil function does not accept
                // numbers like a/b
                const float aux = float(i + 1) / columnsPerBlock;
                const size_t blockOffset = size_t(ceil(aux));

                for (size_t j = i + 1; j < rows; ++j) {
                    if (std::abs(LU[i][i]) < std::numeric_limits<T>::epsilon()) {
                        throw anpi::Exception("error, division by 0 in pivot -> singular matrix cannot perform LU");
                    }

                    // We calc the elimination factor
                    const T eliminationFactor = LU[j][i] / LU(i, i);

                    // The factor is stored in L
                    LU[j][i] = eliminationFactor;


                    // If we almost finish the row there is no need for SIMD operations
                    if (blockOffset < blocksByRow) {

                        // The elimination vector is filled with the factor
                        std::fill(eliminationVector.begin(), eliminationVector.end(), eliminationFactor);

                        // The register is positioned on the required row and is
                        // moved to the nearest block to avoid operating L
                        luRegPtr = luReg + j * blocksByRow + blockOffset;

                        // We position the register in the row before so we can do the multiplication
                        luRegRowBeforePtr = luReg + i * blocksByRow + blockOffset;

                        //The end of the row is defined
                        luEndRegPtr = luReg + (j + 1) * blocksByRow;

                        for (; luRegPtr != luEndRegPtr;) {
                            // First we multiply the row for the elimination factor
                            // then the result is subtracted from the row
                            // and finally its stored in the register
                            *luRegPtr++ = mm_sub<T>(*luRegPtr,
                                                    mm_mult<T, regType>(*luRegRowBeforePtr++, *eliminationRegister));
                        }
                    }

                    // Now we eliminate the elements that could not be eliminated
                    // using SIMD in a secuential manner
                    for (size_t k = i + 1; k < columnsPerBlock * blockOffset; ++k) {
                        LU[j][k] -= eliminationFactor * LU[i][k];
                    }
                }
            }

        }


        // On-copy implementation of gaussian elimination for SIMD-capable types
        template<typename T,
                class Alloc,
                typename std::enable_if<is_simd_type<T>::value, int>::type = 0>
        inline void gaussElimRegType(Matrix <T, Alloc> &LU,
                                     std::vector<size_t> &permut) {


            if (is_aligned_alloc<Alloc>::value && (std::is_same<T, float>::value || std::is_same<T, double>::value)) {
#ifdef __AVX__
                luDoolittleSIMD<T, Alloc, typename avx_traits<T>::reg_type>(LU, permut);
#else
                anpi::fallback::luDoolittle(LU, permut);
#endif
            } else {
                anpi::fallback::luDoolittle(LU, permut);
            }
        }


    }  //namespace simd





    template<typename T,
            class Alloc>
    inline void gaussElimination(Matrix <T, Alloc> &LU,
                                 std::vector<size_t> &permut) {

        if (is_simd_type<T>::value) {
            anpi::simd::gaussElimRegType(LU, permut);
        } else {
            anpi::fallback::luDoolittle(LU, permut);
        }

    }


// The arithmetic implementation (luimpl) namespace
// dispatches to the corresponding methods
#ifdef ANPI_ENABLE_SIMD
    namespace luimpl = simd;
#else
    namespace luimpl = fallback;
#endif

} // namespace anpi


#endif //PROYECTO2_LUAUX_H
