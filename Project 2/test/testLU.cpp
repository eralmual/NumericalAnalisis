/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <boost/test/unit_test.hpp>


#include "LUDoolittle.hpp"

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>

#include <functional>

#include <cmath>

namespace anpi {
    namespace test {

        /// Test the given closed root finder
        template<typename T>
        void luTest(const std::function<void(const Matrix<T> &,
                                             Matrix<T> &,
                                             std::vector<size_t> &)> &decomp,
                    const std::function<void(const Matrix<T> &,
                                             Matrix<T> &,
                                             Matrix<T> &)> &unpack) {

            // The result
            Matrix<T> LU;

            // Test if a non-square matrix is successfully detected
            {
                Matrix<T> A = {{1, 7,  6,  4},
                               {2, 17, 27, 17}};
                std::vector<size_t> p;
                try {
                    decomp(A, LU, p);
                    BOOST_CHECK_MESSAGE(false, "Rectangular matrix not properly catched");
                }
                catch (anpi::Exception &exc) {
                    BOOST_CHECK_MESSAGE(true, "Rectangular matrix properly detected");
                }
            }

            // Test pivoting
            {
                anpi::Matrix<T> A = {{-1, -2, 1, 2},
                                     {2,  0,  1, 2},
                                     {-1, -1, 0, 1},
                                     {1,  1,  1, 1}};
                std::vector<size_t> p;
                decomp(A, LU, p);

                std::vector<size_t> gp = {1, 0, 3, 2};
                BOOST_CHECK(gp == p);
            }

            // Test decomposition
            {
                // same matrix as before, but already permuted to force a
                // clean decomposition
                anpi::Matrix<T> A = {{2,  0,  1, 2},
                                     {-1, -2, 1, 2},
                                     {1,  1,  1, 1},
                                     {-1, -1, 0, 1}};
                std::vector<size_t> p;
                decomp(A, LU, p);
                Matrix<T> L, U;
                unpack(LU, L, U);
                Matrix<T> Ar = L * U;

                const T eps = std::numeric_limits<T>::epsilon();

                BOOST_CHECK(Ar.rows() == A.rows());
                BOOST_CHECK(Ar.cols() == A.cols());

                for (size_t i = 0; i < Ar.rows(); ++i) {
                    for (size_t j = 0; j < Ar.cols(); ++j) {
                        BOOST_CHECK(std::abs(Ar(i, j) - A(i, j)) < eps);
                    }
                }
            }

            // Test bigger decomposition
            {
                // same matrix as before, but already permuted to force a
                // clean decomposition
                anpi::Matrix<T> A = {{2,  3,  1,  9,   0,  -8, -7, 3,   4,  5,  -3},
                                     {-1, -2, 8,  -8,  -1, 5,  5,  9,   -6, -4, 2},
                                     {1,  1,  2,  7,   2,  2,  3,  -1,  8,  6,  4},
                                     {-1, -1, 3,  5,   3,  7,  -1, 7,   3,  7,  -9},
                                     {2,  4,  -6, -6,  -4, 4,  5,  4,   -3, -3, 8},
                                     {-5, -1, 14, 4,   5,  -1, 9,  -6,  2,  2,  -5},
                                     {8,  6,  -9, 3,   6,  0,  8,  2,   10, 1,  -4},
                                     {-4, 8,  -7, 2,   -7, 9,  -5, 8,   6,  4,  7},
                                     {-3, -4, 7,  1,   8,  6,  2,  -13, 11, 7,  -10},
                                     {7,  3,  0,  0,   9,  -3, 4,  0,   8,  8,  3},
                                     {4,  -5, 3,  -11, 10, 8,  -6, 9,   -7, -9, 1},};


                std::cout << "A =" << std::endl << "{";
                for (int i = 0; i < A.rows(); ++i) {
                    std::cout << "{";
                    for (int j = 0; j < A.cols(); ++j) {
                        std::cout << A(i, j) << ",\t";
                    }
                    std::cout << "}" << std::endl;
                }
                std::cout << "}" << std::endl;
                std::cout << std::endl;


                std::vector<size_t> p;
                decomp(A, LU, p);
                Matrix<T> L, U;
                unpack(LU, L, U);
                Matrix<T> Ar = L * U;

                std::cout << "A =" << std::endl << "{";
                for (int i = 0; i < A.rows(); ++i) {
                    std::cout << "{";
                    for (int j = 0; j < A.cols(); ++j) {
                        std::cout << A(p[i], j) << ",\t";
                    }
                    std::cout << "}" << std::endl;
                }
                std::cout << "}" << std::endl;
                std::cout << std::endl;


                std::cout << "Ar =" << std::endl << "{";
                for (int i = 0; i < Ar.rows(); ++i) {
                    std::cout << "{";
                    for (int j = 0; j < Ar.cols(); ++j) {
                        std::cout << Ar(i, j) << ",\t";
                    }
                    std::cout << "}" << std::endl;
                }
                std::cout << "}" << std::endl;
                std::cout << std::endl;

                std::cout << "L =" << std::endl << "{";
                for (int i = 0; i < L.rows(); ++i) {
                    std::cout << "{";
                    for (int j = 0; j < L.cols(); ++j) {
                        std::cout << L(i, j) << ",\t";
                    }
                    std::cout << "}" << std::endl;
                }
                std::cout << "}" << std::endl;
                std::cout << std::endl;

                std::cout << "U =" << std::endl << "{";
                for (int i = 0; i < U.rows(); ++i) {
                    std::cout << "{";
                    for (int j = 0; j < U.cols(); ++j) {
                        std::cout << U(i, j) << ",\t";
                    }
                    std::cout << "}" << std::endl;
                }
                std::cout << "}" << std::endl;
                std::cout << std::endl;

                const T eps = std::numeric_limits<T>::epsilon();

                BOOST_CHECK(Ar.rows() == A.rows());
                BOOST_CHECK(Ar.cols() == A.cols());

                for (size_t i = 0; i < Ar.rows(); ++i) {
                    for (size_t j = 0; j < Ar.cols(); ++j) {
                        std::cout << Ar(i, j) - A(p[i], j) << std::endl;
                        BOOST_CHECK(std::abs(Ar(i, j) - A(p[i], j)) < eps * 100);
                    }
                }
            }
        }

    } // test
}  // anpi

BOOST_AUTO_TEST_SUITE(LU)

    BOOST_AUTO_TEST_CASE(Doolittle) {
        anpi::test::luTest<float>(anpi::luDoolittle<float>,
                                  anpi::unpackDoolittle<float>);
        anpi::test::luTest<double>(anpi::luDoolittle<double>,
                                   anpi::unpackDoolittle<double>);
    }


BOOST_AUTO_TEST_SUITE_END()
