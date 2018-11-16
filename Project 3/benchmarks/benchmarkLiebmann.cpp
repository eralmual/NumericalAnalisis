/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Erick Muñoz Alvarado
 * @Date  : 4.11.2018
 */

#include <boost/test/unit_test.hpp>


#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <Liebmann.hpp>

/**
 * Unit tests for the matrix class
 */
#include "benchmarkFramework.hpp"
#include "Matrix.hpp"
#include "Allocator.hpp"

BOOST_AUTO_TEST_SUITE(LiebmannBench)

/// Benchmark for addition operations
    template<typename T>
    class benchLiebmann {
    protected:
        /// Maximum allowed size for the square matrices
        const size_t _verticalLenght, _horizontalLength;


        /// State of the benchmarked evaluation
        anpi::Matrix<T> _borderConditions = anpi::Matrix<T>(25000, 25000);
        std::vector<bool> _isIsolated = {false, false, false, false};
        bool _isUsingOpenMP;
    public:
        /// Construct
        benchLiebmann(const size_t horizontalLength, const size_t verticalLenght, const bool isUsingOpenMP)
                : _horizontalLength(horizontalLength), _verticalLenght(verticalLenght) {

            _isUsingOpenMP = isUsingOpenMP;
            _borderConditions.fillRow(250, 0);
            _borderConditions.fillRow(125, 1);
            _borderConditions.fillRow(300, 2);
            _borderConditions.fillRow(50, 3);
        }

        /// Prepare the evaluation of given size
        void prepare(const size_t size) {
            assert (size <= 25000);
        }
    };

/// Provide the evaluation method for in-place addition
    template<typename T>
    class parallelLiebmannFloat : public benchLiebmann<T> {
    public:
        /// Constructor
        parallelLiebmannFloat(const size_t horizontalLength,
                              const size_t verticalLenght) : benchLiebmann<T>(horizontalLength, verticalLenght, true) {}

        // Evaluate in-place
        inline void eval() {
            anpi::liebmann<float>(this->_borderConditions,
                                  this->_verticalLenght,
                                  this->_horizontalLength,
                                  this->_isIsolated,
                                  1,
                                  this->_isUsingOpenMP);
        }
    };

/// Provide the evaluation method for on-copy addition
    template<typename T>
    class parallelLiebmannDouble : public benchLiebmann<T> {
    public:
        /// Constructor
        parallelLiebmannDouble(const size_t horizontalLength,
                               const size_t verticalLenght) : benchLiebmann<T>(horizontalLength, verticalLenght,
                                                                               true) {}

        // Evaluate in-place
        inline void eval() {
            anpi::liebmann<double>(this->_borderConditions,
                                   this->_verticalLenght,
                                   this->_horizontalLength,
                                   this->_isIsolated,
                                   1,
                                   this->_isUsingOpenMP);
        }
    };


    /// Provide the evaluation method for in-place addition
    template<typename T>
    class notParallelLiebmannDouble : public benchLiebmann<T> {
    public:
        /// Constructor
        notParallelLiebmannDouble(const size_t horizontalLength,
                                  const size_t verticalLenght) : benchLiebmann<T>(horizontalLength, verticalLenght,
                                                                                  false) {}

        // Evaluate in-place
        inline void eval() {
            anpi::liebmann<float>(this->_borderConditions,
                                  this->_verticalLenght,
                                  this->_horizontalLength,
                                  this->_isIsolated,
                                  1,
                                  this->_isUsingOpenMP);
        }
    };

/// Provide the evaluation method for on-copy addition
    template<typename T>
    class notParallelLiebmannFloat : public benchLiebmann<T> {
    public:
        /// Constructor
        notParallelLiebmannFloat(const size_t horizontalLength,
                                 const size_t verticalLenght) : benchLiebmann<T>(horizontalLength, verticalLenght,
                                                                                 false) {}

        // Evaluate in-place
        inline void eval() {
            anpi::liebmann<double>(this->_borderConditions,
                                   this->_verticalLenght,
                                   this->_horizontalLength,
                                   this->_isIsolated,
                                   1,
                                   this->_isUsingOpenMP);
        }
    };

/**
 * Instantiate and test the methods of the Matrix class
 */
    BOOST_AUTO_TEST_CASE(DoolittleBenchmrk) {

        std::vector<size_t> sizes = {256, 512, 1024, 2048,
                                     4096, 8192, 16384, 25000};

        const size_t n = sizes.back();
        const size_t repetitions = 10;
        std::vector<anpi::benchmark::measurement> times;

        {
            parallelLiebmannFloat<float> baoc(n, n);

            // Measure on-copy add
            ANPI_BENCHMARK(sizes, repetitions, times, baoc);

            ::anpi::benchmark::write("Liebmann_float_OPM.txt", times);
            ::anpi::benchmark::plotRange(times, "FloatOpenMP", "r");
        }

        {
            parallelLiebmannDouble<double> baoc(n,n);

            // Measure on-copy add
            ANPI_BENCHMARK(sizes, repetitions, times, baoc);

            ::anpi::benchmark::write("Liebmann_double_OPM.txt", times);
            ::anpi::benchmark::plotRange(times, "DoubleOpenMP", "g");
        }
        {
            notParallelLiebmannDouble<float> baoc(n,n);

            // Measure on-copy add
            ANPI_BENCHMARK(sizes, repetitions, times, baoc);

            ::anpi::benchmark::write("Liebmann_float_serial.txt", times);
            ::anpi::benchmark::plotRange(times, "FloatSerial", "b");
        }

        {
            notParallelLiebmannFloat<double> baoc(n,n);

            // Measure on-copy add
            ANPI_BENCHMARK(sizes, repetitions, times, baoc);

            ::anpi::benchmark::write("Liebmann_double_serial.txt", times);
            ::anpi::benchmark::plotRange(times, "DoubleSerial", "m");
        }


        ::anpi::benchmark::show();
    }

BOOST_AUTO_TEST_SUITE_END()