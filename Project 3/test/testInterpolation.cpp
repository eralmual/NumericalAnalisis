/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Erick Muñoz Alvarado
 * @Date  : 4.11.2018
 * @File  : testInterpolation.cpp
 */
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <Thomas.hpp>
#include <vector>
#include <Interpolation.hpp>


/**
 * Unit tests for the matrix class
 */

#include "Matrix.hpp"
#include "Allocator.hpp"



BOOST_AUTO_TEST_SUITE(Interpolation)

    template<typename T>
    void testCubicSplines() {
        std::vector<double> x = {0, 12, 24, 49};
        std::vector<double> y = {10, 40, 50, 100};
        std::vector<double> r = {10, 12.5, 15, 17.5, 20, 22.5, 25, 27.5, 30,
                                 32.5, 35, 37.5, 40, 40.8333, 41.6667, 42.5,
                                 43.3333, 44.1667, 45, 45.8333, 46.6667, 47.5,
                                 48.3333, 49.1667, 50, 52, 54, 56, 58, 60, 62,
                                 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86,
                                 88, 90, 92, 94, 96, 98, 100};
        std::vector<double> z = anpi::cubicSplinesInterpolation(x, y, 50);


        T epsilon = std::numeric_limits<T>::epsilon();
        if(std::is_same<T, float>::value){
            epsilon *= 10000;
        }
        if(std::is_same<T, double>::value){
            epsilon = 0.001;
        }

        for (size_t j = 0; j < z.size(); ++j) {
            BOOST_CHECK(abs(z[j] - r[j]) <= epsilon);
        }


    }


/**
 * Instantiate and test the methods of the Matrix class
 */
    BOOST_AUTO_TEST_CASE(cubicSplines) {
        testCubicSplines<double>();
        testCubicSplines<float>();
    }

BOOST_AUTO_TEST_SUITE_END()