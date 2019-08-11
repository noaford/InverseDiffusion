//
//  Grid.hpp
//  InverseDiffusion
//
//  Created by Noah Ford on 12/9/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//

#ifndef Grid_hpp
#define Grid_hpp

#include "uniformmesh2d.h"
#include <stdio.h>


// INDEXING SHORTCUTS
#define sind(IS,I,J) ((IS)*maxi+(I))*maxj+(J)
#define cind(I,J) (I)*maxj+(J)
enum {kSurface=0, // levelset slice
    kSpeed,
    kForward,
    kBackward,
    kmu_s, //Scattering info
    kmu_a, //absorption info
    kD, //Diffusion info
    kindex, //indices of refraction ?
    //kindex_delta, //boundary location
    kincident_radiation, //Incident Radiation ?
    kData, //Boundary Data from foward model
    kSource, //Location of discs
    kDeriv,
    kNorm,
    kdx,
    kdy,
    kdxx,
    kdyy,
    kdxy,
    kcurv,
    kTemp1,
    kTemp2
};

enum {kIsData = 0,
    kitemp1,
    kitemp2,
    kitemp3
};

class Grid : public levelset::UniformMesh2D{
    double SourceStrength = 1.;
    int nX;
    int nY;
    double dx;
    double dy;
    int numdoublelayers;
    int numintlayers;
    
public:
    Grid(const int nX, const int nY, const double dx, const double dy, const int numdoublelayers, const int numintlayers);
    
    ~Grid(void);
    
    void ReadParameters(void);
    
    void InitializeData(double centerx, double centery, double radius, int const temp);
    
    void ForwardDiffusion(int const surface, int const output, int const ktemp);
    
    void BackwardDiffusion(int const source, int const output);
    
    void GetSpeed(double inversionspeedmultiplier);
    
    void InitializeLevelSet(double centerx, double centery, double radius);
    
    void InverseStep(double inversionspeedmultiplier);
        
    void DoInversion(int maxiter, double inversionspeedmultiplier);
    
    void OutputToFile(void);
    
    void GetCurvature(void);
    
    
    
    
};
#endif /* Grid_hpp */
