//
// Created by allan on 19/09/18.
//

#ifndef TAREA04_SUBSTITUTION_H
#define TAREA04_SUBSTITUTION_H

#include "Matrix.hpp"
#include <limits>

namespace anpi {

    template<typename T>
    void backSUB(const anpi::Matrix<T> &U, std::vector<T> &x, const std::vector<T> &result) {
        int n = static_cast<int>(U.rows());
        x.resize(static_cast<unsigned long>(n));
        x[n - 1] = result[n - 1] / U[n - 1][n - 1];
        T sum;
        for (int i = n - 2; i >= 0; --i) {
            sum = T(0);
            for (int j = i + 1; j < n; ++j) {
                sum += U[i][j] * x[j];
            }
            T value = T(1) / U[i][i] * (result[i] - sum);
            x[i] = value;
        }
    }

    template<typename T>
    void forwardSUB(const anpi::Matrix<T> &L, std::vector<T> &x, const std::vector<T> &result) {
        size_t n = L.rows();
        x.resize(n);
        x[0] = result[0] / L[0][0];
        T sum;
        for (size_t i = 1; i < n; ++i) {
            sum = T(0);
            for (size_t j = 0; j < i; ++j) {
                sum += L[i][j] * x[j];
            }
            T value = T(1) / L[i][i] * (result[i] - sum);
            x[i] =value;
        }
    }
}
#endif //TAREA04_SUBSTITUTION_H
