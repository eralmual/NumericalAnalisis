/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Erick Muñoz Alvarado
 * @Date  : 7.11.2018
 * @File  : Spline.hpp
 */

#ifndef PROYECTO3_SPLINE_HPP
#define PROYECTO3_SPLINE_HPP


#include <vector>
#include "Matrix.hpp"
#include <Thomas.hpp>
namespace anpi {


    /**
     * Contains the necessary information to interpolate any value in a given range
     * into a given function
     *
     * @tparam T : Data type
     */
    template<typename T>
    class Spline {
    public:
        Spline() {}

        /// Constructor, recives two vectors one with the x's and other with f(x) to interpolate
        Spline<T>(const std::vector<T> x, const std::vector<T> y);

        /**
         * Gets the interpolated value of a given x
         *
         * @param inputX : Value to interpolate
         * @return       : Interpolated value
         */
        inline T interpolate(const T inputX) {
            size_t i = 1;

            // If the value to interpolate is out of bounds a exeption is raised
            if (inputX > this->x[this->x.size() - 1]) throw anpi::Exception("Value to interpolate is out of bounds");

            // Find the interpolation range
            while (inputX > this->x[i]) {
                ++i;
            }

            // Interpolate using the cubic splines formula
            T result = coefficients(0, i - 1) * pow(inputX - this->x[i], 3)
                       + coefficients(1, i - 1) * pow(inputX - this->x[i - 1], 3)
                       + coefficients(2, i - 1) * (inputX - this->x[i])
                       + coefficients(3, i - 1) * (inputX - this->x[i - 1]);
            return result;
        }

    private:
        // Atributtes
        /// All values of x
        std::vector<T> x;
        /// All values of y
        std::vector<T> y;
        /// All values of the second derivatives
        std::vector<T> ddfx;
        /// All the coefficients used to interpolate
        anpi::Matrix<T> coefficients;

        // Methods
        /**
         * Calculates the derivatives of the function using Thomas algorithm to solve the
         * tridiagonal equation system.
         *
         * @param x : Values of x
         * @param y : Values of y
         */
        inline void calculateDerivatives(const std::vector<T> &x, const std::vector<T> &y) {
            // We have to create a tridiagonal matrix in order to
            // solve the splines, so to avoid wasting memory, lets
            // just create three vectors
            std::vector<T> upper, lower, mid, r;
            size_t size = y.size();

            // Initial conditions are set
            lower.push_back(0);
            upper.push_back(x[2] - x[1]);
            mid.push_back(2 * (x[2] - x[0]));
            r.push_back(6 * ((y[2] - y[1]) / (x[2] - x[1]) -
                             (y[1] - y[0]) / (x[1] - x[0])));

            // Here we fill the upper and middle diagonal but since the
            // middle one has one more element we will have to add it
            // later. Also upper and lower have the same values so, we just
            // make the operations on upper and then copy it to lower
            for (size_t n = 1; n < size - 3; ++n) {

                // Formula to get the upper row
                upper.push_back(x[n + 1] - x[n]);

                // The lower row will always be equal to the upper last row
                lower.push_back(upper[n - 2]);

                // Formula to get the middle diagonal
                mid.push_back(2 * (x[n + 2] - x[n]));

                // Formula to get the r value
                r.push_back(6 * ((y[n + 2] - y[n + 1]) / (x[n + 2] - x[n + 1]) -
                                 (y[n + 1] - y[n]) / (x[n + 1] - x[n])));
            }

            // Final conditions are set
            lower.push_back(x[size - 1] - x[size - 2]);
            mid.push_back(2 * (x[size - 1] - x[size - 3]));
            upper.push_back(0);
            r.push_back(6 * ((y[size - 1] - y[size - 2]) / (x[size - 1] - x[size - 2]) -
                             (y[size - 2] - y[size - 3]) / (x[size - 2] - x[size - 3])));

            // We get the solution to the equation system for the second derivatives
            // using Thomas method
            this->ddfx = thomasDecomposition(lower, mid, upper, r);

            // Lets insert the border derivatives, they are asumed zero for
            // Thomas algorithm
            this->ddfx .push_back(T(0));
            this->ddfx .insert(ddfx.begin(), T(0));
        }

        /**
         * Calculates the coefficients used to interpolate values in teh given range
         */
        inline void calculateCoefficients() {

            // Now lets calculate the coefficient matrix to avoid calculating them in
            // each iteration
            for (size_t i = 1; i < x.size(); ++i) {
                this->coefficients(0, i - 1) = ddfx[i - 1] / (6 * (x[i - 1] - x[i]));
                this->coefficients(1, i - 1) = ddfx[i] / (6 * (x[i] - x[i - 1]));

                this->coefficients(2, i - 1) = y[i - 1] / (x[i - 1] - x[i]) - ddfx[i - 1] * (x[i - 1] - x[i]) / 6;
                this->coefficients(3, i - 1) = y[i] / (x[i] - x[i - 1]) - ddfx[i] * (x[i] - x[i - 1]) / 6;
            }

        }

    };

    /**
     * Constructor of the spline class
     * @tparam T : Data type
     * @param x  : Array of x values
     * @param y  : Array of y values
     */
    template<typename T>
    Spline<T>::Spline(const std::vector<T> x, const std::vector<T> y) {
        // Check the zise of the input data
        if (x.size() != y.size()) {
            throw anpi::Exception("X and Y have a different size");
        }

        // Lets save the input vectors
        this->x = x;
        this->y = y;

        // Calculates the derivatives needed to create the coefficients matrix
        calculateDerivatives(this->x, this->x);

        // Finally lets calculate the coefficients and store them
        this->coefficients = anpi::Matrix<T>(4, x.size() - 1);
        calculateCoefficients();
    }


} // namespace anpi


#endif //PROYECTO3_SPLINE_HPP
