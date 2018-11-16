/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Erick Muñoz Alvarado
 * @Date  : 4.11.2018
 * @File  : Thomas.hpp
 */

#ifndef PROYECTO3_THOMAS_HPP
#define PROYECTO3_THOMAS_HPP


#include <vector>

namespace anpi {

    /**
     * Solves the equation system of a tridiagonal matrix using Thomas method
     *
     * @tparam T    :   Data type template
     * @param lower :   Lower diagonal
     * @param mid   :   Middle diagonal
     * @param upper :   Upper diagonal
     * @param r     :   Result vector
     * @return      :   Vector that when multiplied by a
     *                  tridiagonal matrix composed by
     *                  upper, middle and lower is equal to r
     */
    template<typename T>
    std::vector<T> thomasDecomposition(std::vector<T> lower,
                                       std::vector<T> mid,
                                       std::vector<T> upper,
                                       std::vector<T> r) {

        size_t n = mid.size();
        std::vector<T> x(n);

        // Decomposition
        for (size_t i = 1; i < n; ++i) {

            lower[i] /= mid[i - 1];
            mid[i] -= lower[i]*upper[i - 1];

        }

        // Forward substitution
        for (size_t j = 1; j < n; ++j) {

            r[j] -= lower[j] * r[j - 1];

        }

        // Back substitution
        x[n - 1] = r[n - 1] / mid[n - 1];
        for (long k = n - 2; 0 <= k; --k) {

            x[k] = ((r[k] - upper[k]* x[k + 1]) / mid[k]);

        }

        return x;

    }
}

#endif //PROYECTO3_THOMAS_HPP
