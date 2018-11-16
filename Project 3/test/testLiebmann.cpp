/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 */

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <chrono>


/**
 * Unit tests for the matrix class
 */

#include "Matrix.hpp"
#include "Liebmann.hpp"
#include "Allocator.hpp"


BOOST_AUTO_TEST_SUITE(LiebmannImplementation)


    BOOST_AUTO_TEST_CASE(Liebmann) {
        anpi::Matrix<double> borders = anpi::Matrix<double>(4, 1000);
        borders.fillRow(100, 0);
        borders.fillRow(50, 1);
        borders.fillRow(250, 2);
        borders.fillRow(33, 3);


        std::vector<bool> bordersIsolation = {false, false, false, false};

        std::cout << "Starting Liebmann test..." << std::endl;

        auto start = std::chrono::steady_clock::now();

        anpi::Matrix<double> b = liebmann(borders, 1000, 1000, bordersIsolation);

        auto end = std::chrono::steady_clock::now();
        //b.print();
        std::cout << "Finished Liebmann test in " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;

        //25K x 25K = 67717.6 ms; No OpenM
        //25K x 25K = 48123.4 ms; OpenMP
    }


BOOST_AUTO_TEST_SUITE_END()