/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Erick Muñoz Alvarado
 * @Date  : 7.11.2018
 * @File  : Interpolation.hpp
 */

#ifndef PROYECTO3_INTERPOLATION_HPP
#define PROYECTO3_INTERPOLATION_HPP

#include <vector>
#include "Spline.hpp"

namespace anpi{

    /**
     * Uses linear interpolation to get the values between two given points
     * in a function
     *
     * @tparam T    : Data type
     * @param x1    : Initial value of x
     * @param x2    : Final value of x
     * @param y1    : Function value at initial x
     * @param y2    : Function value at final x
     * @param size  : Amount of intermediate numbers to interpolate
     * @return      : A vector with the interpolated values
     */
    template<typename T>
    std::vector<T> linearInterpolation(T x1, T x2, T y1, T y2, size_t size){

        // We create our result vector
        std::vector<T> y;

        // Define the 'm' and 'b' for the linear ecuation
        T m = (y2 - y1) / (x2 - x1);
        T b = y1 - m*x1;
        T step = (x2 - x1) / (size - 1);

        // Finally we make the calculations for each 'x'
        // and store it in 'y'
        for (size_t i = 0; i < size; ++i) {
            y.push_back(m*(x1 + T(i)*step) + b);
        }

        return y;
    }


    /**
     * Uses cubic splines to interpolate the values in a given array of points
     * in a function
     *
     * @tparam T    :   Data type
     * @param x     :   Vector of x's in a function
     * @param y     :   Vector containing the f(x) value of the given x's
     * @param size  :   Number of values we want to interpolate
     * @return      :   A vector with the interpolated values
     */
    template<typename T>
    std::vector<T> cubicSplinesInterpolation(std::vector<T> x, std::vector<T> y, size_t size) {

        assert(x.size() == y.size());

        // Lets initialize the spline
        anpi::Spline<T> spline(x, y);

        // Create the vector with the interpolated values
        std::vector<T> interpolation(size);
        interpolation[0] = y[0];
        interpolation[size - 1] = y[y.size() - 1];

        // Auxiliary variables
        T step = ceil(size / y.size());
        size_t yPosition = step;
        size_t yIndex = 1;

        for (size_t i = 1; i < size - 1; ++i) {

            if(i != yPosition){

                interpolation[i] = spline.interpolate(i);
                yPosition += step;

            }
            else{
                interpolation[i] = y[yIndex++];
            }

        }

        return interpolation;
    }



}

#endif //PROYECTO3_INTERPOLATION_HPP
