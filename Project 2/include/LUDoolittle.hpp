#include <limits>
#include <algorithm>

#include "Exception.hpp"
#include "Matrix.hpp"
#include "LUAux.hpp"


#ifndef ANPI_LU_DOOLITTLE_HPP
#define ANPI_LU_DOOLITTLE_HPP

namespace anpi {

    /**
     * pivot the element of the diagonal
     * @tparam T
     * @param A matrix
     * @param row current row
     * @param permut permutations vector
     */
    template<typename T>
    void pivot(Matrix<T> &A, size_t row, std::vector<size_t> &permut) {
        //we start assuming the first element is the greatest
        T max = A[row][row];
        size_t maxI = row;
        for (size_t p = row + 1; p < A.rows(); ++p) {
            if (std::abs(A[p][row]) > std::abs(max)) {
                //save the new max
                maxI = p;
                max = A[p][row];
            }
        }
        if (maxI != row) {
            std::swap(permut[row], permut[maxI]);
            anpi::luimpl::swapRows(A, row, maxI);

        }
    }

    /**
     * Auxiliary method used to debug LU decomposition.
     *
     * It separates a packed LU matrix into the lower triangular matrix
     * L and the upper triangular matrix U, such that the diagonal of L
     * is composed by 1's.
     */
    template<typename T>
    void unpackDoolittle(const Matrix<T> &LU,
                         Matrix<T> &L,
                         Matrix<T> &U) {

        L = LU;
        U = LU;
        size_t rows = LU.rows();
        for (size_t i = 0; i < rows; ++i) {//iterate the rows
            for (size_t j = 0; j < rows; ++j) {//iterate the cols
                if (i == j) {
                    L[i][j] = 1;// diagonal of L has 1's
                }
                if (j > i) {
                    L[i][j] = 0;
                }
                if (j < i) {
                    U[i][j] = 0;
                }
            }
        }
    }

    /**
     * Decompose the matrix A into a lower triangular matrix L and an
     * upper triangular matrix U.  The matrices L and U are packed into
     * a single matrix LU.
     *
     * The L matrix will have in the Doolittle's LU decomposition a
     * diagonal of 1's
     *
     * @param[in] A a square matrix
     * @param[out] LU matrix encoding the L and U matrices
     * @param[out] permut permutation vector, holding the indices of the
     *             original matrix falling into the corresponding element.
     *             For example if permut[5]==3 holds, then the fifth row
     *             of the LU decomposition in fact is dealing with the third
     *             row of the original matrix.
     *
     * @throws anpi::Exception if matrix cannot be decomposed, or input
     *         matrix is not square.
     */
    template<typename T>
    void luDoolittle(const Matrix<T> &A,
                     Matrix<T> &LU,
                     std::vector<size_t> &permut) {
        if (A.cols() != A.rows()) {
            throw anpi::Exception("cannot solve rectangular matrix");
        }
        LU = A;
        size_t rows = A.rows();
        permut.resize(rows);// resize de vector
        for (size_t i = 0; i < rows; ++i) {
            //the permutations vector starts with the items ordered from 0 to n
            permut[i] = i;
        }

        //we perform gauss elimination while creating the L and U matrices
        anpi::gaussElimination(LU, permut);
    }

    template<typename T>
    void luFallbackTest(const Matrix<T> &A,
                        Matrix<T> &LU,
                        std::vector<size_t> &permut) {
        if (A.cols() != A.rows()) {
            throw anpi::Exception("cannot solve rectangular matrix");
        }
        LU = A;
        size_t rows = A.rows();
        permut.resize(rows);// resize de vector
        for (size_t i = 0; i < rows; ++i) {
            //the permutations vector starts with the items ordered from 0 to n
            permut[i] = i;
        }

        //we perform gauss elimination while creating the L and U matrices
        anpi::fallback::luDoolittle(LU, permut);
    }

    template<typename T>
    void luSIMDTest(const Matrix<T> &A,
                    Matrix<T> &LU,
                    std::vector<size_t> &permut) {
        if (A.cols() != A.rows()) {
            throw anpi::Exception("cannot solve rectangular matrix");
        }
        LU = A;
        size_t rows = A.rows();
        permut.resize(rows);// resize de vector
        for (size_t i = 0; i < rows; ++i) {
            //the permutations vector starts with the items ordered from 0 to n
            permut[i] = i;
        }

        //we perform gauss elimination while creating the L and U matrices
        anpi::simd::gaussElimRegType(LU, permut);
    }

} // namespace anpi

#endif
