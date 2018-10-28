/**
 * Copyright (C) 2017
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @author Erick Munoz Alvarado
 * @date   01.10.2018
 */


#include <boost/test/unit_test.hpp>


#include <iostream>
#include <exception>
#include <cstdlib>
#include <complex>
#include <LUDoolittle.hpp>

/**
 * Unit tests for the matrix class
 */
#include "benchmarkFramework.hpp"
#include "Matrix.hpp"
#include "Allocator.hpp"

BOOST_AUTO_TEST_SUITE(Matrix)

/// Benchmark for addition operations
    template<typename T>
    class benchLU {
    protected:
        /// Maximum allowed size for the square matrices
        const size_t _maxSize;

        /// A large matrix holding
        anpi::Matrix<T> _data;

        /// State of the benchmarked evaluation
        anpi::Matrix<T> _LU;
        anpi::Matrix<T> _L;
        anpi::Matrix<T> _U;
        std::vector<size_t> _p;
    public:
        /// Construct
        benchLU(const size_t maxSize)
                : _maxSize(maxSize), _data(maxSize, maxSize, anpi::DoNotInitialize) {

            size_t idx = 0;
            for (size_t r = 0; r < _maxSize; ++r) {
                for (size_t c = 0; c < _maxSize; ++c) {
                    //_data(r,c)=idx++;
                    _data(r, c) = r - c + 1;
                }
            }
        }

        /// Prepare the evaluation of given size
        void prepare(const size_t size) {
            assert (size <= this->_maxSize);
        }
    };

/// Provide the evaluation method for in-place addition
    template<typename T>
    class decomposeLUDoolittleFloatFallback : public benchLU<T> {
    public:
        /// Constructor
        decomposeLUDoolittleFloatFallback(const size_t n) : benchLU<T>(n) {}

        // Evaluate add in-place
        inline void eval() {
            anpi::luFallbackTest<float>(this->_data, this->_LU, this->_p);
        }
    };

/// Provide the evaluation method for on-copy addition
    template<typename T>
    class decomposeLUDoolittleDoubleFallback : public benchLU<T> {
    public:
        /// Constructor
        decomposeLUDoolittleDoubleFallback(const size_t n) : benchLU<T>(n) {}

        // Evaluate add in-place
        inline void eval() {
            anpi::luFallbackTest<double>(this->_data, this->_LU, this->_p);
        }
    };


    /// Provide the evaluation method for in-place addition
    template<typename T>
    class decomposeLUDoolittleFloatSIMD : public benchLU<T> {
    public:
        /// Constructor
        decomposeLUDoolittleFloatSIMD(const size_t n) : benchLU<T>(n) {}

        // Evaluate add in-place
        inline void eval() {
            anpi::luSIMDTest<float>(this->_data, this->_LU, this->_p);
        }
    };

/// Provide the evaluation method for on-copy addition
    template<typename T>
    class decomposeLUDoolittleDoubleSIMD : public benchLU<T> {
    public:
        /// Constructor
        decomposeLUDoolittleDoubleSIMD(const size_t n) : benchLU<T>(n) {}

        // Evaluate add in-place
        inline void eval() {
            anpi::luSIMDTest<double>(this->_data, this->_LU, this->_p);
        }
    };

/**
 * Instantiate and test the methods of the Matrix class
 */
    BOOST_AUTO_TEST_CASE(DoolittleBenchmrk) {

        std::vector<size_t> sizes = {24, 32, 48, 64,
                                     96, 128, 192, 256,
                                     384};/*, 512, 768,1024};,
                                       1536,2048,3072,4096};*/

        const size_t n = sizes.back();
        const size_t repetitions = 50;
        std::vector<anpi::benchmark::measurement> times;

        {
            decomposeLUDoolittleFloatFallback<float> baoc(n);

            // Measure on-copy add
            ANPI_BENCHMARK(sizes, repetitions, times, baoc);

            ::anpi::benchmark::write("doolittle_float_fb.txt", times);
            ::anpi::benchmark::plotRange(times, "FloatFallback", "r");
        }

        {
            decomposeLUDoolittleDoubleFallback<double> baoc(n);

            // Measure on-copy add
            ANPI_BENCHMARK(sizes, repetitions, times, baoc);

            ::anpi::benchmark::write("doolittle_double_fb.txt", times);
            ::anpi::benchmark::plotRange(times, "DoubleFallback", "g");
        }
        {
            decomposeLUDoolittleFloatSIMD<float> baoc(n);

            // Measure on-copy add
            ANPI_BENCHMARK(sizes, repetitions, times, baoc);

            ::anpi::benchmark::write("doolittle_float_simd.txt", times);
            ::anpi::benchmark::plotRange(times, "FloatSIMD", "b");
        }

        {
            decomposeLUDoolittleDoubleSIMD<double> baoc(n);

            // Measure on-copy add
            ANPI_BENCHMARK(sizes, repetitions, times, baoc);

            ::anpi::benchmark::write("doolittle_double_simd.txt", times);
            ::anpi::benchmark::plotRange(times, "DoubleSIMD", "m");
        }


        ::anpi::benchmark::show();
    }

BOOST_AUTO_TEST_SUITE_END()