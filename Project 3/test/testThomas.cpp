/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Erick Muñoz Alvarado
 * @Date  : 4.11.2018
 * @File  : testThomas.cpp
 */
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <Thomas.hpp>
#include <vector>


/**
 * Unit tests for the matrix class
 */

#include "Matrix.hpp"
#include "Allocator.hpp"



BOOST_AUTO_TEST_SUITE(Thomas)

    template<typename T>
    void testThomasDecomposition() {
        std::vector<T> l = {0, -1, -1, -1};
        std::vector<T> m = {2.04, 2.04, 2.04, 2.04};
        std::vector<T> u = {-1, -1, -1, 0};
        std::vector<T> r = {40.08, 0.8, 0.8, 200.8};
        std::vector<T> answer = {65.4256, 93.3883, 124.286, 159.356};
        std::vector<T> x = anpi::thomasDecomposition(l,
                                        m,
                                        u,
                                        r);
        T epsilon = std::numeric_limits<T>::epsilon();
        if(std::is_same<T, float>::value){
            epsilon *= 10000;
        }
        if(std::is_same<T, double>::value){
            epsilon = 0.001;
        }

        for (size_t j = 0; j < x.size(); ++j) {
            BOOST_CHECK(abs(x[j] - answer[j]) <= epsilon);
        }


    }


/**
 * Instantiate and test the methods of the Matrix class
 */
    BOOST_AUTO_TEST_CASE(Decomposition) {
        testThomasDecomposition<double>();
        testThomasDecomposition<float>();
    }

BOOST_AUTO_TEST_SUITE_END()