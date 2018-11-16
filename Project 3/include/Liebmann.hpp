/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, ITCR, Costa Rica
 *
 * This file is part of the numerical analysis lecture CE3102 at TEC
 *
 * @Author: Erick Muñoz Alvarado
 * @Date  : 4.11.2018
 */


#ifndef ANPI_LIEBMANN_H
#define ANPI_LIEBMANN_H


#include <cstdlib>
#include <functional>
#include <cmath>

#include "Matrix.hpp"
#include "Exception.hpp"

namespace anpi {

    /**
     *  Used to get the average value in a chunk on the border pixels, this
     *  chunk is alligned with the size of the pixel block used on the iteration
     *  and changes are saved on the lastIteration matrix since this matrix
     *  is the one used to calculate the new values in the next iteration.
     *
     *
     * @tparam T            :   Data type
     * @tparam Alloc        :   Allocator used for row allignment in the matrix values
     * @param lastIteration :   Result of the operations on the previous iteration on the matrix
     * @param isIsolated    :   Vector describing if a border is isolated. Convention used: {top; bot; left; right}
     * @param rowIndexes    :   Array of positions in which the matrix rows are grouped
     * @param columnIndexes :   Array of positions in which the matrix columns are grouped
     * @param sizes         :   Size of the arrays.
     */
    template<typename T, class Alloc>
    void fixFrontierConditions(anpi::Matrix<T, Alloc> &lastIteration,
                               const std::vector<bool> isIsolated,
                               T *rowIndexes,
                               T *columnIndexes,
                               size_t sizes) {

        // We set the indexes to get the mean of the frontier conditions
        size_t iStart = 1;
        size_t iEnd = *rowIndexes++;
        size_t jStart = 1;
        size_t jEnd = *columnIndexes++;


        for (size_t i = 0; i < sizes; ++i) {


            // Fix chunk of rows
            lastIteration.fillRow(lastIteration.averageRow(2 * isIsolated.at(0),
                                                           jStart,
                                                           jEnd), // Mean of the frontier conditions in range
                                  0,
                                  jStart,                       // Determine the range
                                  jEnd);


            lastIteration.fillRow(lastIteration.averageRow(lastIteration.rows() - (1 + 2 * isIsolated.at(1)),
                                                           jStart,
                                                           jEnd), // Mean of the frontier conditions in range
                                  lastIteration.rows() - 1,
                                  jStart,                       // Determine the range
                                  jEnd);


            // Fix chunk of columns
            lastIteration.fillColumn(lastIteration.averageColumn(2 * isIsolated.at(2),
                                                                 iStart,
                                                                 iEnd), // Mean of the frontier conditions in range
                                     0,
                                     iStart,                       // Determine the range
                                     iEnd);

            lastIteration.fillColumn(lastIteration.averageColumn(lastIteration.cols() - (1 + 2 * isIsolated.at(3)),
                                                                 iStart,
                                                                 iEnd), // Mean of the frontier conditions in range
                                     lastIteration.cols() - 1,
                                     iStart,                       // Determine the range
                                     iEnd);


            // Here we we get boundaries for the next iteration ready
            jStart = jEnd;
            iStart = iEnd;
            if (i + 2 == sizes) {
                jEnd = lastIteration.cols() - 1;
                iEnd = lastIteration.rows() - 1;
            } else {
                jEnd = *columnIndexes++;
                iEnd = *rowIndexes++;
            }

        }

    }


    /**
     *  Makes the calculation for the new value of a chunk of pixels
     *
     * @tparam T                : Data type
     * @tparam Alloc            : Allocator used for row allignment in the matrix values
     * @param operationMatrix   : Matrix where we write the result of the operations
     * @param lastIteration     : Operation matrix after the operations on the last iteration and the borders fixed
     * @param iStart            : First row of the chunk
     * @param iEnd              : Limit row of the chunk
     * @param jStart            : First column of the chunk
     * @param jEnd              : Limit column of the chunk
     * @param lambda            : Relaxation coefficient
     */
    template<typename T, class Alloc>
    void operateOnChunk(anpi::Matrix<T, Alloc> &operationMatrix,
                        const anpi::Matrix<T, Alloc> &lastIteration,
                        const size_t iStart,
                        const size_t iEnd,
                        const size_t jStart,
                        const size_t jEnd,
                        const T lambda) {


        // Calculate the values using the neighborhood pixels of the chunk
        T newValue = (lastIteration(iEnd, jStart) + lastIteration(iStart - 1, jStart) + lastIteration(iStart, jEnd) +
                    lastIteration(iStart, jStart - 1)) / 4;


        for (size_t i = iStart; i < iEnd; ++i) {

            for (size_t j = jStart; j < jEnd; ++j) {

                // Write the value in the operation matrix chunk
                operationMatrix(i, j) = lambda * newValue + (1 - lambda) * lastIteration(i, j);

            }
        }

    }


    /**
     * Gets the indexes to divide the matrix rows by half on each iteration,first row is first
     * iteration and so on, the las iteration is the maximum row division possible before dividing
     * the matrix rows in groups of one
     *
     * @tparam T                : Data type
     * @tparam Alloc            : Allocator used for row allignment in the matrix values
     * @param operationMatrix   : Matrix we will be dividing
     * @return                  : Matrix containing the indexes we need to divide in each iteration
     */
    template<typename T, class Alloc>
    anpi::Matrix<T, Alloc> getRowIndexMatrix(const anpi::Matrix<T, Alloc> &operationMatrix) {
        // This is used to determine how many rows and columns the index matrix needs
        size_t n = ceil(log2(operationMatrix.rows())) - 1;
        anpi::Matrix<T, Alloc> rowIndex = anpi::Matrix<T, Alloc>(n, pow(2, n) + 1, T(1));

        // The first row is set
        rowIndex(0, 1) = operationMatrix.rows() / 2;
        rowIndex(0, 2) = operationMatrix.rows() - 1;

        // we now calc the index cuts needed to divide the rows in each iteration
        // the index starts at 1 because we need the initial conditions set befores in order
        // for the algorithm to work properly
        for (size_t i = 1; i < n; ++i) {

            // This index is an auxiliary index to copy elements in the row before
            // without losing their location
            T previousRowIndex = 1;

            for (size_t j = 1; j < pow(2, i) + 1; ++j) {

                // Completes the first half of the matrix

                // When
                if (j % 2 != 0) {
                    rowIndex(i, j) = ceil(rowIndex(i - 1, previousRowIndex - 1) +
                                          (rowIndex(i - 1, previousRowIndex) - rowIndex(i - 1, previousRowIndex - 1)) /
                                          2);


                } else {
                    rowIndex(i, j) = rowIndex(i - 1, previousRowIndex);
                    ++previousRowIndex;
                }

                // Completes the other half of the matrix
                rowIndex(i, pow(2, i + 1) - j) = operationMatrix.rows() - rowIndex(i, j);
                //std::cout << "RI: " << previousRowIndex << "   j: " << j << "   i: " << i << std::endl;
            }

            // We always start from the position 1 and en in the pos rows - 1
            // to avoid operating the border conditions
            rowIndex(i, pow(2, i + 1)) = operationMatrix.rows() - 1;
        }

        return rowIndex;
    }


    /**
     * Auxiliary function to Liebmann(*), it operates on chunks of data until maximum division is achived,
     * then it starts to iterate each pixel individually and finished when the substraction of last iteration
     * and the current is lower or equal to a threshold
     *
     * @tparam T                : Data type
     * @tparam Alloc            : Allocator used for row allignment in the matrix values
     * @param operationMatrix   : Matrix in the which we will write our calculations
     * @param isIsolated        : Vector describing if a border is isolated. Convention used: {top; bot; left; right}
     * @param lambda            : Relaxation coefficient
     * @param isUsingOpenMP     : Flag to activate openMP
     */
    template<typename T, class Alloc>
    void liebmannAux(anpi::Matrix<T, Alloc> &operationMatrix,
                     const std::vector<bool> isIsolated,
                     const T lambda,
                     const bool isUsingOpenMP = true) {

        // We create the auxiliary variables
        anpi::Matrix<T, Alloc> lastIteration, rowIndex, columnIndex;
        T factor = 0;
        size_t limit;
        bool isTransposed = false;

        // Lets check which dimension is the smallest
        // we always assume that rows < cols, in case this is not
        // true, we transpose the matrix and work normally
        // a flag is set to transpose the matrix again in the end
        if (operationMatrix.rows() > operationMatrix.cols()) {
            operationMatrix = operationMatrix.copyTransposed();
            isTransposed = true;
        }


        // Create the last iteration copy to read from it and
        // gets the index necessary to reduce the rows and the columns to
        // individual pixels
        lastIteration = operationMatrix;
        rowIndex = getRowIndexMatrix(operationMatrix);
        columnIndex = getRowIndexMatrix(operationMatrix.copyTransposed());


        // We iterate over the rows of rowIndex, they contain the indexes
        // to separate the matrix in row chunks, at the same time we also
        // iterate over columnIndex who hat the same function but for the
        // columns
        for (size_t i = 0; i < rowIndex.rows(); ++i) {

            // We get the mean of the frontier conditions for them to align
            // to the chunk that will need their information
            fixFrontierConditions(lastIteration,
                                  isIsolated,
                                  rowIndex[i],
                                  columnIndex[i],
                                  size_t(pow(2, i) - 1));

            // we need to save the limit in a variable in order to make
            // OpenMP work

            limit = size_t(pow(2, i + 1));



            // lets activate OpenMP when the matrix has more than 100 rows and
            // the user says they want to use it
#ifdef ANPI_ENABLE_OpenMP
#pragma omp parallel for num_threads(4) if (isUsingOpenMP)
#endif
            for (size_t j = 0; j <= limit; ++j) {

                for (size_t k = 0; k <= limit; ++k) {

                    // We calculate the value of a chunk delimited by the indexes
                    operateOnChunk(operationMatrix,
                                   lastIteration,
                                   rowIndex(i, j),
                                   rowIndex(i, j + 1),
                                   columnIndex(i, k),
                                   columnIndex(i, k + 1),
                                   lambda);

                }

            }

            // We make a copy of our matrix for the calculations of the next iteration
            lastIteration = operationMatrix;

        }


        // This is in case we can divide the matrix even more before running pixel by pixel
        // that means operationMatrix.cols > operationMatrix.rows
        if (operationMatrix.rows() != operationMatrix.cols()) {
            size_t rowIndexRow = rowIndex.rows() - 1;
            factor = 15; //this is to increase speed
            size_t rowLimit = rowIndex.cols();

            for (size_t i = rowIndexRow + 1; i < columnIndex.rows(); ++i) {
                // We get the mean of the frontier conditions for them to align
                // to the chunk that will need their information
                fixFrontierConditions(lastIteration,
                                      isIsolated,
                                      rowIndex[rowIndexRow],
                                      columnIndex[i],
                                      size_t(pow(2, i) - 1));

                // we need to save the limit in a variable in order to make
                // OpenMP work
                limit = size_t(pow(2, i + 1));

                // lets activate OpenMP when the matrix has more than 100 rows and
                // the user says they want to use it
#ifdef ANPI_ENABLE_OpenMP
#pragma omp parallel for num_threads(4) if (isUsingOpenMP)
#endif
                for (size_t j = 0; j < rowLimit; ++j) {

                    for (size_t k = 0; k <= limit; ++k) {

                        // We calculate the value of a chunk delimited by the indexes
                        operateOnChunk(operationMatrix,
                                       lastIteration,
                                       rowIndex(rowIndexRow, j),
                                       rowIndex(rowIndexRow, j + 1),
                                       columnIndex(i, k),
                                       columnIndex(i, k + 1),
                                       lambda);

                    }

                }

                // We make a copy of our matrix for the calculations of the next iteration
                if (i + 1 < columnIndex.rows()) {
                    lastIteration = operationMatrix;
                }


            }
        }

        T ip1, im1, jp1, jm1;
        //operationMatrix.print('O');
        // Finally we make and individual pixel run until convergence is achieved
        // Reduce the factor to increase precision
        while (!operationMatrix.hasConverged(lastIteration, factor)) {
            T newValue;
            lastIteration = operationMatrix;

            for (size_t i = 1; i < operationMatrix.rows() - 1; ++i) {

                for (size_t j = 1; j < operationMatrix.cols() - 1; ++j) {

                    // Top and bot isolation
                    if(isIsolated.at(0)){
                        ip1 = lastIteration(i - 1, j);
                    }
                    else{
                        ip1 = lastIteration(i + 1, j);
                    }
                    if(isIsolated.at(1)){
                        im1 = lastIteration(i + 1, j);
                    }
                    else{
                        im1 = lastIteration(i - 1, j);
                    }

                    // left and right
                    if(isIsolated.at(2)){
                        jp1 = lastIteration(i, j + 1);
                    }
                    else{
                        jp1 = lastIteration(i, j - 1);
                    }
                    if(isIsolated.at(3)){
                        jm1 = lastIteration(i, j - 1);
                    }
                    else{
                        jm1 = lastIteration(i, j + 1);
                    }

                    newValue = (ip1 + im1 + jp1 + jm1) / 4;

                    operationMatrix(i, j) = lambda * newValue + (1 - lambda) * lastIteration(i, j);

                }
            }
        }


        //  If the matrix was transposed we return it as it was
        if (isTransposed) {
            operationMatrix = operationMatrix.copyTransposed();
        }

    }


    /**
     * Master function it takes the frontier conditions and the isolation vector to get the average of
     * vales and set it as the initial value for all the matrix
     * @tparam T                : Data type
     * @tparam Alloc            : Allocator used for row allignment in the matrix values
     * @param frontierConditions: Sets the state of the borders of the plate. Convention used: {top; bot; left; right}
     * @param verticalLength    : Number of rows in the operation matrix
     * @param horizontalLength  : Number of columns in the operation matrix
     * @param isIsolated        : Vector describing if a border is isolated. Convention used: {top; bot; left; right}
     * @param lambda            : Relaxation coefficient
     * @param isUsingOpenMP     : Flag signaling the use of OpenMP
     * @return                  : A matrix with the heat distribution given the border conditions
     */
    template<typename T, class Alloc>
    anpi::Matrix<T, Alloc> liebmann(anpi::Matrix<T, Alloc> &frontierConditions,
                                    const size_t verticalLength,
                                    const size_t horizontalLength,
                                    const std::vector<bool> isIsolated,
                                    T lambda = 1,
                                    const bool isUsingOpenMP = true) {

        // The +2 is added to insert the frontierConditions into the matrix
        const size_t cols = horizontalLength + 2;
        const size_t rows = verticalLength + 2;

        // Lets initialize the result matrix
        T averageFrontierCondition = 0;
        T numBorders = 0;

        // Take into account the frontier conditions at the start
        // Top
        if(!isIsolated[0]){
            averageFrontierCondition  += frontierConditions.averageRow(0, 0, horizontalLength);
            ++numBorders;
        }
        else if(!isIsolated[1]){
            averageFrontierCondition  += frontierConditions.averageRow(1, 0, horizontalLength);
            ++numBorders;
        }

        // Bottom
        if(!isIsolated[1]){
            averageFrontierCondition  += frontierConditions.averageRow(1, 0, horizontalLength);
            ++numBorders;
        }
        else if(!isIsolated[0]){
            averageFrontierCondition  += frontierConditions.averageRow(0, 0, horizontalLength);
            ++numBorders;
        }

        // Left
        if(!isIsolated[2]){
            averageFrontierCondition  += frontierConditions.averageRow(2, 0, verticalLength);
            ++numBorders;
        }
        else if(!isIsolated[3]){
            averageFrontierCondition  += frontierConditions.averageRow(3, 0, verticalLength);
            ++numBorders;
        }

        // Right
        if(!isIsolated[3]){
            averageFrontierCondition  += frontierConditions.averageRow(3, 0, verticalLength);
            ++numBorders;
        }
        else if(!isIsolated[2]){
            averageFrontierCondition  += frontierConditions.averageRow(2, 0, verticalLength);
            ++numBorders;
        }

        // All the plate is isolated
        else if(isIsolated[0] && isIsolated[1] && isIsolated[2] && isIsolated[3]){
            throw anpi::Exception("El sistema esta completamente aislado");
        }

        // cal the average
        averageFrontierCondition /= numBorders;

        anpi::Matrix<T, Alloc> operationMatrix = anpi::Matrix<T, Alloc>(rows, cols, averageFrontierCondition);

        // We set the top and bottom conditions
        // the first and last columns are ignored
        // because they are frontier conditions
        for (size_t j = 1; j < cols - 1; ++j) {
            operationMatrix(0, j) = frontierConditions(0, j - 1);           // Top condition is copied
            operationMatrix(rows - 1, j) = frontierConditions(1, j - 1);    // Bottom condition is copied
        }

        // Lets now continue with the left and right
        for (size_t i = 1; i < rows - 1; ++i) {
            operationMatrix(i, 0) = frontierConditions(2, i - 1);           // Left condition is copied
            operationMatrix(i, cols - 1) = frontierConditions(3, i - 1);    // Right condition is copied

        }

        liebmannAux(operationMatrix, isIsolated, lambda, isUsingOpenMP);

        return operationMatrix;

    }


} //namespace anpi



#endif //ANPI_LIEBMANN_H
