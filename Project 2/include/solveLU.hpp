//
// Created by allan on 19/09/18.
//

#ifndef TAREA04_SOLVELU_H
#define TAREA04_SOLVELU_H

#include "Matrix.hpp"
#include "LUDoolittle.hpp"
#include "Substitution.hpp"

namespace anpi {
    template<typename T>
    inline void lu(const anpi::Matrix<T> &A,
                   anpi::Matrix<T> &LU,
                   std::vector<size_t> &p) {

        anpi::luDoolittle(A, LU, p);
    }

    template<typename T>
    bool solveLU(const anpi::Matrix<T> &A,
                 std::vector<T> &x,
                 const std::vector<T> &b) {

        anpi::Matrix<T> LU;
        std::vector<size_t> p;
        anpi::lu(A, LU, p);


        std::vector<T> bpermuted(b.size());
        for (size_t i = 0; i < b.size(); ++i) {
            bpermuted[i] = b[p[i]];
        }
        anpi::Matrix<T> L, U;
        anpi::unpackDoolittle(LU, L, U);

        std::vector<T> y;
        anpi::forwardSUB(L, y, bpermuted);
        anpi::backSUB(U, x, y);

        return true;

    }

    template<typename T>
    void invert(const anpi::Matrix<T> &A,
                anpi::Matrix<T> &Ai) {

        Ai.allocate(A.rows(), A.cols());
        anpi::Matrix<T> I;
        I.allocate(A.rows(), A.cols());
        I.fill(T(0));
        std::vector<T> x;
        x.resize(A.cols());
        for (int i = 0; i < A.rows(); i++) {
            I[i][i] = 1;
        }
        for (int i = 0; i < A.rows(); i++) {
            anpi::solveLU(A, x, I.column(i));
            for (int j = 0; j < A.rows(); j++) {
                Ai[j][i] = x[j];
            }
        }
    }
}

#endif //TAREA04_SOLVELU_H
