/**
 * Copyright (C) 2018
 * Área Académica de Ingeniería en Computadoras, TEC, Costa Rica
 *
 * This file is part of the CE3102 Numerical Analysis lecture at TEC
 *
 * @Author: 
 * @Date  : 24.02.2018
 */

#include <cstdlib>
#include <iostream>
#include <Matrix.hpp>
#include <ResistorGrid.hpp>


int main(int argc, char *argv[]) {

    if(argc!=2){
        throw anpi::Exception("run with ./proyecto2 MapName.png");
    }
    anpi::node ni{5, 0}, nf{10, 28};
    anpi::ResistorGrid rg(argv[1], ni, nf);

    rg.build();
    rg.solveCurrents();

    rg.navigate();
    rg.plotNav();

    return EXIT_SUCCESS;
}
  
