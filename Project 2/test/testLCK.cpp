/**
 * Copyright (C) 2017 
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: Pablo Alvarado
 * @Date  : 10.02.2018
 */

#include <boost/test/unit_test.hpp>
#include "LUDoolittle.hpp"
#include "ResistorGrid.hpp"

namespace anpi {
    namespace test {
        void LCK(const std::string &filename, node &ni, node &nf, float tolerance) {

            anpi::ResistorGrid rg(filename, ni, nf);
            rg.build();
            rg.solveCurrents();


            indexPair
                    leftNode{0, 0, 0, 0},
                    rightNode{0, 0, 0, 0},
                    upNode{0, 0, 0, 0},
                    dwnNode{0, 0, 0, 0};

            size_t leftIdx = 0, rightIdx = 0, upIdx = 0, dwnIdx = 0;

            float upI, dnI, lI, rI;
            float saliendo, entrando, suma;

            for (size_t i = 0; i < rg.rawMap_.rows(); ++i) {
                for (size_t j = 0; j < rg.rawMap_.cols(); ++j) {
                    //we set the values to 0;
                    upI = 0.f;
                    dnI = 0.f;
                    lI = 0.f;
                    rI = 0.f;
                    saliendo = 0.f;
                    entrando = 0.f;
                    suma = 0.f;

                    if (i > 0) {//if not the first row we determine the upper resistor
                        upNode.row1 = i - 1;
                        upNode.row2 = i;
                        upNode.col1 = upNode.col2 = j;
                        upIdx = rg.nodesToIndex(upNode);
                        upI = rg.currents[upIdx];
                        if (upI < 0.f) {
                            saliendo += std::abs(upI);
                        } else {
                            entrando += upI;
                        }


                    }
                    if (i != rg.rawMap_.rows() - 1) {//if not the last row we determine the resistor below the node
                        dwnNode.row1 = i;
                        dwnNode.row2 = i + 1;
                        dwnNode.col1 = dwnNode.col2 = j;
                        dwnIdx = rg.nodesToIndex(dwnNode);
                        dnI = rg.currents[dwnIdx];
                        if (dnI < 0.) {
                            entrando += std::abs(dnI);
                        } else {
                            saliendo += dnI;
                        }
                    }
                    if (j > 0) {//if not the first column
                        leftNode.row1 = leftNode.row2 = i;
                        leftNode.col1 = j - 1;
                        leftNode.col2 = j;
                        leftIdx = rg.nodesToIndex(leftNode);
                        lI = rg.currents[leftIdx];
                        if (lI < 0.) {
                            saliendo += std::abs(lI);
                        } else {
                            entrando += lI;
                        }

                    }
                    if (j != rg.rawMap_.cols() - 1) {//if not the last column
                        rightNode.row1 = rightNode.row2 = i;
                        rightNode.col1 = j;
                        rightNode.col2 = j + 1;
                        rightIdx = rg.nodesToIndex(rightNode);
                        rI = rg.currents[rightIdx];
                        if (rI < 0.f) {
                            entrando += std::abs(rI);
                        } else {
                            saliendo += rI;
                        }
                    }
                    suma = std::abs(entrando - saliendo);
                    if ((i == ni.i && j == ni.j) || (i == nf.i && j == nf.j)) {
                        BOOST_CHECK_MESSAGE(std::abs(suma - 1) <= tolerance,
                                            "sum of currents not 1 in initial or final node");

                    } else {
                        BOOST_CHECK_MESSAGE(suma <= tolerance, "sum of currents not 0");
                    }

                }
            }

            rg.navigate();
            rg.plotNav();
        }
    }
}

BOOST_AUTO_TEST_SUITE(ResistorGrid)

    BOOST_AUTO_TEST_CASE(LCK) {
        std::cout << "LCK test: checking LCK in resistor grid....\n";
        anpi::node ni{5, 0}, nf{10, 28};
        anpi::test::LCK("mapa25x29.png", ni, nf, 1.0e-5f);

    }


BOOST_AUTO_TEST_SUITE_END()
