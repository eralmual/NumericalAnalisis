#include <utility>

//
// Created by allan on 22/09/18.
//

#ifndef PROYECTO2_MAP_H
#define PROYECTO2_MAP_H

#include "Matrix.hpp"
#include "Exception.hpp"
#include <AnpiConfig.hpp>
#include <limits>

#include <string>

#include <opencv2/core.hpp>    // For cv::Mat
#include <opencv2/highgui.hpp> // For cv::imread/imshow
#include <iostream>
#include "solveLU.hpp"
#include "PlotPy.hpp"

namespace anpi {
    float epsFloat = std::numeric_limits<float>::epsilon();
    struct indexPair {
        std::size_t row1;
        std::size_t col1;
        std::size_t row2;
        std::size_t col2;
    };

    struct node {
        size_t i, j;
    };

    class ResistorGrid {
    private:
        const float highR = 1000000.f;
        const float lowR = 1.f;
        std::string filename;
        size_t blockSize = std::numeric_limits<size_t>::quiet_NaN();
        size_t numberOfResistors = std::numeric_limits<size_t>::quiet_NaN();
        node &initialNode, &finalNode;
        std::vector<float> b_{};


    public:
        Matrix<float> A_;
        std::vector<float> currents;

        ResistorGrid(std::string path, node &pInitialNode, node &pFinalNode);


        /**
         * Construct the grid from the given @param filename
         * @return true if successful or false otherwise
         */
        bool build(const std::string &filename);

        /**
         * Compute the internal data to navigate between the given @param nodes
         * @return
         */
        bool navigate();

        /**
         * Converts a pair of nodes to a linear index
         * @param row1
         * @param col1
         * @param row2
         * @param col2
         * @return
         */
        std::size_t nodesToIndex(std::size_t row1,
                                 std::size_t col1,
                                 std::size_t row2,
                                 std::size_t col2);

        /**
         * Converts an index to a pair of nodes containing the resistor wanted
         * @param idx
         * @return the indexes of the starting and ending node
         */
        indexPair indexToNodes(std::size_t idx);

        size_t nodesToIndex(indexPair &indexPair1);

        void populateA();

        void solveCurrents();

        void checkValues();

        bool build();

        Matrix<float> rawMap_{};

        //------------------------------------------------------------------------------------------------------------------------------
        //Pasos en las columnas(x)
        std::vector<float> stepsCol;
        //pasos en las filas(y)
        std::vector<float> stepsRow;
        //public:
        Matrix<float> componente_x;
        Matrix<float> componente_y;
        float magnitud_mayor;
        Matrix<float> matriz_magnitudes;
        Matrix<float> componente_x_normalizada;
        Matrix<float> componente_y_normalizada;
        float alpha = 0.2;
        float eps = 0.00001;
        bool nav2 = false;

        float arreglarX(const indexPair &nodes, const int y, const int x);

        float arreglarY(const indexPair &nodes, const int y, const int x);

        void normalizar();

        void componenteX();

        void componenteY();

        void plotNav();

        float BilinearInterpolation(float q11, float q12, float q21, float q22, float x1, float x2, float y1, float y2,
                                    float x, float y);

        indexPair compararIndices(const indexPair &nodes, const std::vector<float> &x);

        std::size_t nodesToIndex2(std::size_t row1,
                                  std::size_t col1,
                                  std::size_t row2,
                                  std::size_t col2);

        float magnitud_max();
    };
}

namespace anpi {
    ResistorGrid::ResistorGrid(std::string path, node &pInitialNode, node &pFinalNode)
            : filename(std::move(path)), initialNode{pInitialNode}, finalNode{pFinalNode} {
        // se calcula el numero total de resistencias
        if ((pInitialNode.i == 0 && pInitialNode.j == 0) ||
            (pFinalNode.i == 0 && pFinalNode.j == 0)) {
            throw Exception("error: must not choose pixel 0,0");

        }
        componente_x.allocate(rawMap_.rows(), rawMap_.cols());
        componente_y.allocate(rawMap_.rows(), rawMap_.cols());
        componente_x_normalizada.allocate(rawMap_.rows(), rawMap_.cols());
        componente_y_normalizada.allocate(rawMap_.rows(), rawMap_.cols());
    }

    size_t ResistorGrid::nodesToIndex(const std::size_t row1, const std::size_t col1, const std::size_t row2,
                                      const std::size_t col2) {

        //we validate that the coordinates are valid
        if ((row1 == row2 && col1 == col2) ||
            (row2 != row1 && col1 != col2) ||
            (row2 < row1 || col2 < col1) ||
            (row2 > row1 + 1 || col2 > col1 + 1) ||
            (row2 >= rawMap_.rows()) ||
            (col2 >= rawMap_.cols())) {
            throw anpi::Exception("check indexes");
        }
        size_t index = blockSize * row1;
        if (row1 == row2) {
            index += col2;
        } else {
            index += rawMap_.cols() + col2;
        }

        return index - 1;
    }

    size_t ResistorGrid::nodesToIndex(struct indexPair &indexPair1) {
        return ResistorGrid::nodesToIndex(indexPair1.row1, indexPair1.col1, indexPair1.row2, indexPair1.col2);
    }

    indexPair ResistorGrid::indexToNodes(const std::size_t idx) {
        if (idx >= numberOfResistors) {
            throw anpi::Exception("index greater than number of resistors");
        }
        size_t realQuant = idx + 1;//indexes begin in 0, we want the actual number of element

        size_t block = realQuant / blockSize;//calculate the number of block we are in
        size_t remainder = realQuant % blockSize;

        if (remainder < rawMap_.cols()) {
            //we are in a row
            return indexPair{block, remainder - 1, block, remainder};
        } else {
            //we are in a column
            size_t colNumber = remainder - rawMap_.cols();

            return indexPair{block, colNumber, block + 1, colNumber};
        }
    }

    bool ResistorGrid::build(const std::string &filename) {
        // Build the name of the image in the data path
        std::string mapPath = std::string(ANPI_DATA_PATH) + "/" + filename;

        // Read the image using the OpenCV
        cv::Mat_<float> map;

        cv::imread(mapPath.c_str(),
                   CV_LOAD_IMAGE_GRAYSCALE).convertTo(map, CV_32FC1);
        map /= 255.0f; // normalize image range to 0 .. 255


        // Convert the OpenCV matrix into an anpi matrix
        // We have to use the std::allocator to avoid an exact stride
        anpi::Matrix<float, std::allocator<float> > amapTmp(static_cast<const size_t>(map.rows),
                                                            static_cast<const size_t>(map.cols),
                                                            map.ptr<float>());
        // And transform it to a SIMD-enabled matrix
        anpi::Matrix<float> amap(amapTmp);
        rawMap_ = amap;

        // se calcula el numero total de resistencias
        numberOfResistors = 2 * rawMap_.cols() * rawMap_.rows() - rawMap_.cols() - rawMap_.rows();
        blockSize = 2 * rawMap_.cols() - 1;

        return true;
    }

    void ResistorGrid::populateA() {
        b_.resize(numberOfResistors, 0.f);
        A_.allocate(numberOfResistors, numberOfResistors);
        //fill de A matrix with 0's
        A_.fill(0.f);

        //we add the first nxm equations to the A matrix (kirchoff current law
        size_t leftIdx = 0, rightIdx = 0, upIdx = 0, dwnIdx = 0;
        indexPair
                leftNode{0, 0, 0, 0},
                rightNode{0, 0, 0, 0},
                upNode{0, 0, 0, 0},
                dwnNode{0, 0, 0, 0};

        size_t current = 0;//keeps the current current row index of the A matrix
        for (size_t i = 0; i < rawMap_.rows(); ++i) {
            for (size_t j = 0; j < rawMap_.cols(); ++j) {
                if (i > 0) {//if not the first row we determine the upper resistor
                    upNode.row1 = i - 1;
                    upNode.row2 = i;
                    upNode.col1 = upNode.col2 = j;
                    upIdx = nodesToIndex(upNode);
                    A_[current][upIdx] = 1.f;//positive current entering from the top

                }
                if (i != rawMap_.rows() - 1) {//if not the last row we determine the resistor below the node
                    dwnNode.row1 = i;
                    dwnNode.row2 = i + 1;
                    dwnNode.col1 = dwnNode.col2 = j;
                    dwnIdx = nodesToIndex(dwnNode);
                    A_[current][dwnIdx] = -1.f;//negative current flowing out of the node to the bottom
                }
                if (j > 0) {//if not the first column
                    leftNode.row1 = leftNode.row2 = i;
                    leftNode.col1 = j - 1;
                    leftNode.col2 = j;
                    leftIdx = nodesToIndex(leftNode);
                    A_[current][leftIdx] = 1.f;//positive current entering from the left

                }
                if (j != rawMap_.cols() - 1) {//if not the last column
                    rightNode.row1 = rightNode.row2 = i;
                    rightNode.col1 = j;
                    rightNode.col2 = j + 1;
                    rightIdx = nodesToIndex(rightNode);
                    A_[current][rightIdx] = -1.f;
                }
                if (i == initialNode.i && j == initialNode.j) {
                    b_[current] = -1.f;
                }
                if (i == finalNode.i && j == finalNode.j) {
                    b_[current] = 1.f;
                }
                ++current;//we increment the index of the row of the A matrix
            }
        }

        size_t r3_index, r4_index;
        float dwnValue, rightValue, r3Value, r4Value;
        indexPair
                r3_node{0, 0, 0, 0},
                r4_node{0, 0, 0, 0};
        size_t idx = 0;
        //we continue with the kirchhoff voltage law
        for (size_t i = 0; i < rawMap_.rows() - 1; ++i) {
            for (size_t j = 0; j < rawMap_.cols() - 1; ++j) {
                dwnNode.row1 = i;
                dwnNode.row2 = i + 1;
                dwnNode.col1 = dwnNode.col2 = j;
                dwnIdx = nodesToIndex(dwnNode);
                dwnValue = (std::abs(rawMap_[i][j] * rawMap_[i + 1][j] - 1.f) < epsFloat) ? lowR
                                                                                          : highR;//1Mohm or 1 ohm depending on pixel value
                rightNode.row1 = rightNode.row2 = i;
                rightNode.col1 = j;
                rightNode.col2 = j + 1;
                rightIdx = nodesToIndex(rightNode);
                rightValue = (std::abs(rawMap_[i][j] * rawMap_[i][j + 1] - 1.f) < epsFloat) ? lowR
                                                                                            : highR;//1Mohm or 1 ohm depending on pixel value

                r3_node.col1 = r3_node.col2 = rightNode.col2;
                r3_node.row1 = i;
                r3_node.row2 = i + 1;
                r3_index = nodesToIndex(r3_node);
                r3Value = (std::abs(rawMap_[i][j + 1] * rawMap_[i + 1][j + 1] - 1.f) < epsFloat) ? lowR
                                                                                                 : highR;//1Mohm or 1 ohm depending on pixel value

                r4_node.row1 = r4_node.row2 = i + 1;
                r4_node.col1 = j;
                r4_node.col2 = j + 1;
                r4_index = nodesToIndex(r4_node);
                r4Value = (std::abs(rawMap_[i + 1][j] * rawMap_[i + 1][j + 1] - 1.f) < epsFloat) ? lowR
                                                                                                 : highR;//1Mohm or 1 ohm depending on pixel value

                if (i == 0 && j == 0) {
                    //replace the first Ax=b equation with a net equation
                    for (size_t col = 0; col < A_.cols(); ++col) {
                        A_[0][col] = 0.f;//we restore the row with 0's
                    }
                    idx = 0;
                } else {
                    idx = current;
                    ++current;//we increment the index of the row of the A matrix
                }
                //net equation
                A_[idx][dwnIdx] = -dwnValue;
                A_[idx][rightIdx] = rightValue;
                A_[idx][r3_index] = r3Value;
                A_[idx][r4_index] = -r4Value;
            }
        }
    }

    void ResistorGrid::solveCurrents() {
        currents.resize(numberOfResistors, 0.f);
        int count = 0;
        for (auto b:b_) {
            if (std::abs(b) >= epsFloat) count++;
        }
        if (count != 2)throw Exception("check currents at start and end");

        std::cout << "solving currents......." << std::endl;
        solveLU(A_, currents, b_);
        std::cout << "currents solved." << std::endl;
    }

    /**
     * this method checks that there's no column o row with all
     */
    void ResistorGrid::checkValues() {
        for (size_t i = 0; i < A_.rows(); ++i) {
            bool flag = false;
            for (size_t j = 0; j < A_.cols(); ++j) {
                flag = flag || (std::abs(A_[i][j]) >= epsFloat);
            }
            if (!flag)throw anpi::Exception("null row");
        }
        for (size_t j = 0; j < A_.cols(); ++j) {
            bool flag = false;
            std::vector<float> column = A_.column(j);
            for (auto item : column) {
                flag = flag || (std::abs(item) >= epsFloat);
            }
            if (!flag)throw anpi::Exception("null col");
        }

    }

    bool ResistorGrid::build() {
        std::cout << "building rawMap from image ......." << std::endl;
        ResistorGrid::build(filename);
        std::cout << "rawMap built." << std::endl;
        std::cout << "inserting equations in A.........." << std::endl;
        populateA();
        std::cout << "matrix A built." << std::endl;
        std::cout << "checking A for errors......" << std::endl;
        checkValues();
        std::cout << "successful build." << std::endl;
        return true;
    }
    //-------------------------------------------------------------------------------------------------------------------------------------------------


    size_t ResistorGrid::nodesToIndex2(const std::size_t row1, const std::size_t col1, const std::size_t row2,
                                       const std::size_t col2) {


        //we validate that the coordinates are valid
        if ((row1 == row2 && col1 == col2) ||
            (row2 != row1 && col1 != col2) ||
            //(row2 < row1 || col2 < col1) ||
            (row2 > row1 + 1 || col2 > col1 + 1) ||
            (row2 >= rawMap_.rows()) ||
            (col2 >= rawMap_.cols())) {
            plotNav();
            throw anpi::Exception("check indexes");
        }

        size_t index = blockSize * row1;
        if (row1 == row2 && col1 <= col2) {
            index += col2;
        } else if (row1 == row2 && col1 > col2) {
            index += col2 + 1;
        } else if (row1 > row2) {
            index = blockSize * row2;
            index += rawMap_.cols() + col2;
        } else {
            index += rawMap_.cols() + col2;
        }

        return index - 1;
    }


    bool ResistorGrid::navigate() {
        anpi::indexPair trayectoria;
        trayectoria.row1 = initialNode.i;
        trayectoria.col1 = initialNode.j;
        trayectoria.row2 = finalNode.i;
        trayectoria.col2 = finalNode.j;
        nav2 = true;
        if (trayectoria.col1 == trayectoria.col2 &&
            trayectoria.row1 == trayectoria.row2) {
            //true si ya se esta en la posicion que se queria.
            return true;
        }
        if (alpha > 1 || alpha <= 0) {
            throw anpi::Exception("alpha incorrecto");
        }
        componente_x.allocate(rawMap_.rows(), rawMap_.cols());
        componente_y.allocate(rawMap_.rows(), rawMap_.cols());
        componente_x_normalizada.allocate(rawMap_.rows(), rawMap_.cols());
        componente_y_normalizada.allocate(rawMap_.rows(), rawMap_.cols());
        matriz_magnitudes.allocate(rawMap_.rows(), rawMap_.cols());
        //genero matriz de vectores en x
        componenteX();
        //genero matriz de vectores en y
        componenteY();
        //busco la mayor magnitud de todos los vectores de la matriz y lo divido las matrices x y y entre esa magnitud
        normalizar();
        //posicion inicial en x y y
        float y = trayectoria.row1;
        float x = trayectoria.col1;
        //componentes normalizadas con la magnitud mas grande
        float d_x = componente_x[trayectoria.row1][trayectoria.col1];
        float d_y = componente_y[trayectoria.row1][trayectoria.col1];

        int limit = 0;
        while (!(floor(x) == trayectoria.col2 && floor(-y + rawMap_.rows()) == rawMap_.rows() - trayectoria.row2 - 1) &&
               limit < (rawMap_.rows() * rawMap_.cols()) * 5) {
            limit++;
            //reviso si ya almaceno las posiciones para posterior plot pero antes reviso si hay algun nan
            if (isnan(x) || isnan(y) || isnan(d_y) || isnan(d_y)) {
                break;
            } else {
                stepsCol.push_back(x);
                stepsRow.push_back(-y + rawMap_.rows());
                //std::cout << "Step: "<<limit-1 << " x= " << floor(x) << " y= " << -y+rawMap_.rows() << std::endl;
            }
            y += d_y * alpha;
            x += d_x * alpha;
            //para mejorar la velocidad, si el paso es muy pequeno entonces mejor se usa la magnitud del mismo
            if (abs(y - floor(y)) < eps && abs(x - floor(x)) < 0.1) {
                d_y = componente_y[(int) floor(y)][(int) floor(x)];
                d_x = componente_x[(int) floor(y)][(int) floor(x)];
            } else if (abs(d_y) < eps && abs(d_x) < eps) {
                y += arreglarY(trayectoria, y, x);
                x += arreglarX(trayectoria, y, x);
            } else {
                //finalmente si el paso no es tan pequen
                d_x = BilinearInterpolation(componente_x_normalizada[(int) (y)][(int) (x)],
                                            componente_x_normalizada[(int) (y) + 1][(int) (x)],
                                            componente_x_normalizada[(int) (y)][(int) (x) + 1],
                                            componente_x_normalizada[(int) (y) + 1][(int) (x) + 1],
                                            x, x + 1, y, y + 1, x, y);
                d_y = BilinearInterpolation(componente_y_normalizada[(int) y][(int) x],
                                            componente_y_normalizada[(int) y + 1][(int) x],
                                            componente_y_normalizada[(int) y][(int) x + 1],
                                            componente_y_normalizada[(int) y + 1][(int) x + 1],
                                            x, x + 1, y, y + 1, x, y);
            }
        }
        stepsCol.push_back(trayectoria.col2);
        stepsRow.push_back(rawMap_.rows() - trayectoria.row2);
        return true;
    }

    float ResistorGrid::arreglarX(const indexPair &nodes, const int y, const int x) {
        if (rawMap_[y][x + 1] == 0) {
            return -magnitud_mayor;
        }
        if (rawMap_[y][x - 1] == 0) {
            return magnitud_mayor;
        } else {
            return 0;
        }
    }

    float ResistorGrid::arreglarY(const indexPair &nodes, const int y, const int x) {
        if (rawMap_[y + 1][x] == 0) {
            return magnitud_mayor;
        } else if (rawMap_[y - 1][x] == 0) {
            return -magnitud_mayor;
        } else {
            return 0;
        }
    }


    void ResistorGrid::normalizar() {
        float magnitud = magnitud_max();
        componente_x_normalizada = componente_x / magnitud;
        componente_y_normalizada = componente_y / magnitud;
    }

    float ResistorGrid::magnitud_max() {
        float temp = 0;
        for (int i = 0; i < rawMap_.rows(); i++) {
            for (int j = 0; j < rawMap_.cols(); j++) {
                matriz_magnitudes[i][j] = sqrt(
                        componente_x[i][j] * componente_x[i][j] + componente_y[i][j] * componente_y[i][j]);
                if (temp < matriz_magnitudes[i][j]) {
                    temp = matriz_magnitudes[i][j];
                }
            }
        }
        //std::cout <<"Magnitud max:"<< temp << std::endl;
        return temp;
    }

    void ResistorGrid::componenteX() {
        for (int i = 0; i < rawMap_.rows(); i++) {//rows
            for (int j = 0; j < rawMap_.cols(); j++) {//cols
                if (j == 0) {
                    componente_x[i][j] = currents[nodesToIndex2(i, j, i, j + 1)];
                } else if (j == rawMap_.cols() - 1) {
                    componente_x[i][j] = currents[nodesToIndex2(i, j, i, j - 1)];
                } else {
                    componente_x[i][j] =
                            currents[nodesToIndex2(i, j, i, j - 1)] + currents[nodesToIndex2(i, j, i, j + 1)];
                }
            }
        }
    }

    void ResistorGrid::componenteY() {
        for (int i = 0; i < rawMap_.rows(); i++) {
            for (int j = 0; j < rawMap_.cols(); j++) {
                if (i == 0) {
                    componente_y[i][j] = currents[nodesToIndex2(i, j, i + 1, j)];
                } else if (i == rawMap_.rows() - 1) {
                    componente_y[i][j] = currents[nodesToIndex2(i, j, i - 1, j)];
                } else {
                    componente_y[i][j] =
                            currents[nodesToIndex2(i, j, i - 1, j)] + currents[nodesToIndex2(i, j, i + 1, j)];
                }
            }
        }
    }


    void ResistorGrid::plotNav() {
        std::cout << "Iniciando plot" << std::endl;
        anpi::Plot2d<float> plotter;
        plotter.initialize(0);
        plotter.plot2(stepsCol, stepsRow, rawMap_.cols(), rawMap_.rows(), std::string(ANPI_DATA_PATH) + "/" + filename);
        if (nav2) {
            plotter.plotQuiver(componente_x_normalizada, componente_y_normalizada);
        }
        plotter.show();
    }

    float ResistorGrid::BilinearInterpolation(float q11, float q12, float q21, float q22, float x1, float x2, float y1,
                                              float y2, float x, float y) {
        float x2x1, y2y1, x2x, y2y, yy1, xx1;
        x2x1 = x2 - x1;
        y2y1 = y2 - y1;
        x2x = x2 - x;
        y2y = y2 - y;
        yy1 = y - y1;
        xx1 = x - x1;
        return 1 / (x2x1 * y2y1) * (
                q11 * x2x * y2y +
                q21 * xx1 * y2y +
                q12 * x2x * yy1 +
                q22 * xx1 * yy1
        );
    }

}
namespace anpi {
    template<typename T>
    void printMatrix(Matrix<T> &M) {
        size_t rows = M.rows();
        size_t cols = M.cols();
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                std::cout << M[i][j] << ",";
            }
            std::cout << "\n";
        }
    }
}


#endif //PROYECTO2_MAP_H
