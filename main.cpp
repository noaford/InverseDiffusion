//
//  main.cpp
//  InverseDiffusion
//
//  Created by Noah Ford on 12/7/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//

#include <iostream>

#include "Grid.hpp"

int main(int argc, const char * argv[]) {
    //Initialize Grid
    int gridsizex = 100;
    int gridsizey = 100;
    double dx = .01;
    double dy = .01;
    int numdoublelayers = 20;
    int numintlayers = 4;
    
    Grid* theGrid = new Grid(gridsizex,gridsizey,dx,dy,numdoublelayers,numintlayers);
    
    //theGrid -> ReadParameters();
    
    //Do forward diffusion problem to get data
    double sourceX = .75;
    double sourceY = .75;
    double sourceR = .2;
    theGrid->InitializeData(sourceX, sourceY, sourceR, kTemp1);
    
    
    //Initialize Level Set
    double initializeX = .5;
    double initializeY = .5;
    double initializeR = .3;
    theGrid->InitializeLevelSet(initializeX,initializeY,initializeR);
    
    double inversionspeedmultiplier = 100; //Relaxation Parameter
    int maxiter = 100; //Number of iterations in inversion
    theGrid->DoInversion(maxiter, inversionspeedmultiplier);
    
    //Output forward diffusion and level set to file
    theGrid->OutputToFile();
    
    //Clean up the Grid
    delete theGrid;

    std::cout << "Finished With Simulation" << std::endl;
    return 0;
    
    
}
