//
//  Grid.cpp
//  InverseDiffusion
//
//  Created by Noah Ford on 12/9/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//

#include "Grid.hpp"
#include "um2boundary.h"
#include "um2linear.h"
#include "um2periodic.h"
#include "um2xperiodic.h"
#include "um2yperiodic.h"

#define sind(IS,I,J) ((IS)*maxi+(I))*maxj+(J)
#define cind(I,J) (I)*maxj+(J)

using namespace levelset;

Grid::Grid(const int nX, const int nY, const double dx, const double dy, const int numdoublelayers, const int numintlayers): nX(nX),nY(nY),dx(dx),dy(dy),numdoublelayers(numdoublelayers),numintlayers(numintlayers),
UniformMesh2D(StepSize,nX,nY,numdoublelayers,numintlayers,dx,dy,0.0,0.0,*(new UM2_LinearBdry(2,2,2,2) /* inputs are for ghost point width*/))
/* These inputs are X grid size, Y grid size, number of double layers, number of integer layers, dx, dy, x reference coordinate, y reference coordinat, boundary condition)*/ {

    
}

Grid::~Grid(){
    
}

void Grid::ReadParameters(){
    // File handles for parameters
    FILE *mu_s_fid;
    FILE *mu_a_fid;
    FILE *index_fid;
    FILE *index_delta_fid;
    FILE *disc_position_fid;
    FILE *incident_radiation_fid;
    
    double* mu_s_in = new double[(maxi)*(maxj)];
    double* mu_a_in = new double[(maxi)*(maxj)];
    double* index_in = new double[(maxi)*(maxj)];
    int* Is_Data_in = new int[(maxi)*(maxj)];
    int* Source_in = new int[(maxi)*(maxj)];
    double* incident_radiation_in = new double[(maxi)*(maxj)];
    
    //Read in Parameters
    mu_s_fid=fopen("mus_N250_cs3_hsp75.txt","r");
    mu_a_fid=fopen("mua_N250_ca1_hap25.txt","r");
    index_fid=fopen("index1p3.txt","r");
    index_delta_fid=fopen("index_change_bool.txt","r");
    disc_position_fid=fopen("disc_mesh.txt","r");
    incident_radiation_fid=fopen("incident_radiation.txt","r");
    
    fread(mu_s_in, sizeof(double), (maxi)*(maxj), mu_s_fid);
    fread(mu_a_in, sizeof(double), (maxi)*(maxj),mu_a_fid);
    fread(index_in, sizeof(double), (maxi)*(maxj),index_fid);
    fread(Is_Data_in, sizeof(int), (maxi)*(maxj),index_delta_fid);
    fread(Source_in, sizeof(int), (maxi)*(maxj),disc_position_fid);
    fread(incident_radiation_in, sizeof(double), (maxi)*(maxj),incident_radiation_fid);
    
    
    fclose(mu_s_fid);
    fclose(mu_a_fid);
    fclose(index_fid);
    fclose(index_delta_fid);
    fclose(disc_position_fid);
    
    double g = 0.; //REPLACE WITH REAL VALUE
    for(int i =0; i< maxi; i++){
        for( int j=0; j<maxj; j++){
            data_(i,j,kmu_s) = mu_s_in[cind(i,j)];
            data_(i,j,kmu_a) = mu_a_in[cind(i,j)];
            data_(i,j,kindex) = index_in[cind(i,j)];
            idata_(i,j,kIsData) = Is_Data_in[cind(i,j)];
            data_(i,j,kSource) = -(double)Source_in[cind(i,j)] +.5;
            data_(i,j,kincident_radiation) = incident_radiation_in[cind(i,j)];
            
            
            data_(i,j,kD) = (1.)/3./((1.-g)*data_(i,j,kmu_s) + data_(i,j,kmu_a));
        }
    }
    SourceStrength = 1.;
    
    delete[] mu_s_in;
    delete[] mu_a_in;
    delete[] index_in;
    delete[] Is_Data_in;
    delete[] Source_in;
    delete[] incident_radiation_in;
    
}


void Grid::InitializeData(double centerx, double centery, double radius, int const temp){
 
    //double centerx = 3.*dx*maxi/4.;
    //double centery = 3.*dy*maxj/4.;
   
    for(int i=0; i<maxi; i++){
        for(int j=0;j<maxj;j++){
            data_(i,j,temp) = pow(i*dx - centerx,2) + pow(j*dy-centery,2) - pow(min(centery,centerx),2)*pow(radius,2);
            data_(i,j,kD) = 1.;
            data_(i,j,kmu_a) = 0.;
            idata_(i,j,kIsData) = 0;
        }
    }
    
    for(int i=1; i<maxi-1; i++){
        idata_(i,1,kIsData) = 1;
        idata_(i, maxj-2, kIsData) = 1;
    }
    for(int j=2;j<maxj-2;j++){
        idata_(1,j,kIsData) = 1;
        idata_(maxi-2, j, kIsData) = 1;
    }
    
    
    Reinitialize(temp, kSource, kitemp1, kitemp2);
    double minedgevalue = 1.;
    for(int i=1; i<maxi-1; i++){
        minedgevalue = std::min(minedgevalue,data_(i,1,kSource));
        minedgevalue = std::min(minedgevalue,data_(i,maxj-2,kSource));
    }
    for(int j=2;j<maxj-2;j++){
        minedgevalue = std::min(minedgevalue,data_(1,j,kSource));
        minedgevalue = std::min(minedgevalue, data_(maxi-2,j,kSource));
    }
    if(minedgevalue < 0.){
        std::cout << "Source overlaps the domain boundary. Make sure source is fully within domain." << std::endl;
        exit(1);
    }
    
    ForwardDiffusion(kSource,kData,temp);
    
}

void Grid::InitializeLevelSet(double centerx, double centery, double radius){
    
    
    //double centerx = dx*maxi/2.;
    //double centery = dy*maxj/2.;
    
    for(int i=0; i<maxi;i++)
        for(int j=0; j<maxj;j++)
            data_(i,j,kTemp1) = pow(i*dx - centerx,2) + pow(j*dy-centery,2) - pow(min(centery,centerx),2)*pow(radius,2);
    
    Reinitialize(kTemp1, kSurface, kitemp1, kitemp2);
    double minedgevalue = 1.;
    for(int i=1; i<maxi-1; i++){
        minedgevalue = std::min(minedgevalue,data_(i,1,kSurface));
        minedgevalue = std::min(minedgevalue,data_(i,maxj-2,kSurface));
    }
    for(int j=2;j<maxj-2;j++){
        minedgevalue = std::min(minedgevalue,data_(1,j,kSurface));
        minedgevalue = std::min(minedgevalue, data_(maxi-2,j,kSurface));
    }
    if(minedgevalue < 0.){
        std::cout << "Initialized Surface overlaps the domain boundary. Make sure this surface is fully within domain." << std::endl;
        exit(1);
    }
    
    
}

void Grid::ForwardDiffusion(int const surface, int const output, int const ktemp){
    double omega = 1.5;
    
    //double coeff = D*(2/dx/dx + 2/dy/dy);
    //double coeffx = D/dx/dx;
    //double coeffy = D/dy/dy;
    
    for(int i=0; i<maxi; i++){
        for(int j=0;j<maxj;j++){
            data_(i,j,ktemp) = (data_(i,j,surface) < 0.) ? SourceStrength : 0.;
        }
    }
    
    int maxiters = 10000;
    int iters;
    for(iters = 0; iters < maxiters; iters++){
        double maxchange = 0.;
        for(int i=1; i<maxi-1; i++){
            for(int j=1;j<maxj-1;j++){
                double coeff = (data_(i-1,j,kD)/2. + data_(i,j,kD) + data_(i+1,j,kD)/2.)/dx/dx
                         + (data_(i,j-1,kD)/2. + data_(i,j,kD) + data_(i,j+1,kD)/2.)/dy/dy + data_(i,j,kmu_a);
                double coeffxl = (data_(i-1,j,kD) + data_(i,j,kD))/2./dx/dx;
                double coeffxr = (data_(i+1,j,kD) + data_(i,j,kD))/2./dx/dx;
                double coeffyd = (data_(i,j-1,kD) + data_(i,j,kD))/2./dy/dy;
                double coeffyu = (data_(i,j+1,kD) + data_(i,j,kD))/2./dy/dy;
                
                double change = (data_(i,j,ktemp) + data_(i-1,j,output)*coeffxl + data_(i+1,j,output)*coeffxr + data_(i,j-1,output)*coeffyd + data_(i,j+1,output)*coeffyu)/coeff - data_(i,j,output);
                data_(i,j,output) = omega*change + data_(i,j,output);
                maxchange = std::max(std::abs(change),std::abs(maxchange));
            }
        }
        //Update Boundaries
        for(int i=1; i<maxi-1; i++){
            //double Reff = .5; //REPLACE WITH ACTUAL VALUE
            //double Az = 2*data_(i,1,kD)*(1+Reff)/(1-Reff);
            //data_(i,0,output) = data_(i,1,output)*(1-dy/Az);
            
            //Az = 2*data_(i,maxj-2,kD)*(1+Reff)/(1-Reff);
            //data_(i,maxj-1,output) = data_(i,maxj-2,output)*(1-dy/Az);
            
            
            data_(i,0,output) = 0.;
            
            data_(i,maxj-1,output) = 0.;
            
        }
        for(int j=1;j<maxj-1;j++){
            //double Reff = .5; //REPLACE WITH ACTUAL VALUE
            //double Az = 2*data_(1,j,kD)*(1+Reff)/(1-Reff);
            //data_(0,j,output) = data_(1,j,output)*(1-dx/Az);
            
            //Az = 2*data_(maxi-2,j,kD)*(1+Reff)/(1-Reff);
            //data_(maxi-1,j,output) = data_(maxi-2,j,output)*(1-dy/Az);
            data_(0,j,output) = 0.;
            data_(maxi-1,j,output) = 0.;
        }
        
        //Update Corner Points
        data_(0,0,output) = (data_(1,0,output) + data_(0,1,output))/2;
        data_(0,maxj-1,output) = (data_(1,maxj-1,output) + data_(0,maxj-2,output))/2;
        data_(maxi-1,0,output) = (data_(maxi-2,0,output) + data_(maxi-1,1,output))/2;
        data_(maxi-1,maxj-1,output) = (data_(maxi-2,maxj-1,output) + data_(maxi-1,maxj-2,output))/2;
        
        if(maxchange < 10e-8){
            break;
        }
    }
    std::cout << "Forward problem breaking on iteration " << iters << std::endl;
}

void Grid::BackwardDiffusion(int const source, int const output){
    double omega = 1.5;
    int maxiters = 10000;
    int iters;
    for(iters = 0; iters < maxiters; iters++){
        double maxchange = 0.;
        for(int i=1; i<maxi-1; i++){
            for(int j=1;j<maxj-1;j++){
                double coeff = (data_(i-1,j,kD)/2. + data_(i,j,kD) + data_(i+1,j,kD)/2.)/dx/dx
                + (data_(i,j-1,kD)/2. + data_(i,j,kD) + data_(i,j+1,kD)/2.)/dy/dy + data_(i,j,kmu_a);
                double coeffxl = (data_(i-1,j,kD) + data_(i,j,kD))/2./dx/dx;
                double coeffxr = (data_(i+1,j,kD) + data_(i,j,kD))/2./dx/dx;
                double coeffyd = (data_(i,j-1,kD) + data_(i,j,kD))/2./dy/dy;
                double coeffyu = (data_(i,j+1,kD) + data_(i,j,kD))/2./dy/dy;
                
                double change = (data_(i,j,source) + data_(i-1,j,output)*coeffxl + data_(i+1,j,output)*coeffxr + data_(i,j-1,output)*coeffyd + data_(i,j+1,output)*coeffyu)/coeff - data_(i,j,output);
                data_(i,j,output) = omega*change + data_(i,j,output);
                maxchange = std::max(std::abs(change),std::abs(maxchange));
            }
        }
        //Update Boundaries
        for(int i=1; i<maxi-1; i++){
            //double Reff = .5; //REPLACE WITH ACTUAL VALUE
            //double Az = 2*data_(i,1,kD)*(1+Reff)/(1-Reff);
            //data_(i,0,output) = data_(i,1,output)*(1-dy/Az);
            
            //Az = 2*data_(i,maxj-2,kD)*(1+Reff)/(1-Reff);
            //data_(i,maxj-1,output) = data_(i,maxj-2,output)*(1-dy/Az);
            
            data_(i,0,output) = 0.;
            data_(i,maxj-1,output) = 0.;
        }
        for(int j=1;j<maxj-1;j++){
            //double Reff = .5; //REPLACE WITH ACTUAL VALUE
            //double Az = 2*data_(1,j,kD)*(1+Reff)/(1-Reff);
            //data_(0,j,output) = data_(1,j,output)*(1-dx/Az);
            
            //Az = 2*data_(maxi-2,j,kD)*(1+Reff)/(1-Reff);
            //data_(maxi-1,j,output) = data_(maxi-2,j,output)*(1-dy/Az);
            data_(0,j,output) = 0.;
            data_(maxi-1,j,output) =0.;
        }
        
        //Update Corner Points
        data_(0,0,output) = (data_(1,0,output) + data_(0,1,output))/2;
        data_(0,maxj-1,output) = (data_(1,maxj-1,output) + data_(0,maxj-2,output))/2;
        data_(maxi-1,0,output) = (data_(maxi-2,0,output) + data_(maxi-1,1,output))/2;
        data_(maxi-1,maxj-1,output) = (data_(maxi-2,maxj-1,output) + data_(maxi-1,maxj-2,output))/2;
        
        //Break Condition
        if(maxchange < 10e-10){
            break;
        }
    }
    
    std::cout << "Adjoint problem breaking on iteration " << iters << std::endl;
}


void Grid::GetSpeed(double inversionspeedmultiplier){
    
    double epsilon = .00001/inversionspeedmultiplier; //curvature speed
    
    
    ForwardDiffusion(kSurface,kForward,kTemp1);

    for(int i=0; i<maxi;i++)
        for(int j=0; j<maxj;j++)
            data_(i,j,kTemp1) = 0.;
            
    //Store error
    for(int i=1; i<maxi-1;i++)
        for(int j=1; j<maxj-1;j++)
            data_(i,j,kTemp1) = (data_(i,j,kForward) - data_(i,j,kData))*idata_(i, j, kIsData)/dx/dy;
    
    BackwardDiffusion(kTemp1, kBackward);
        
    GetCurvature();
    
    for(int i=1; i<maxi-1;i++)
        for(int j=1; j<maxj-1;j++){
            data_(i,j,kSpeed) = SourceStrength*data_(i,j,kBackward) + epsilon*data_(i,j,kcurv);
        }
    
    //ExtendVelocity(kTemp1, kSpeed, kitemp2, kitemp3, kitemp1);
    
    
    
}

void Grid::GetCurvature(void){
    
    Dx_zero(kSurface,kdx);
    Dxx_zero(kSurface,kdxx);
    Dy_zero(kSurface,kdy);
    Dyy_zero(kSurface,kdyy);
    Dxy_zero(kSurface,kdxy);
    Curvature(kdx, kdy, kdxx, kdxy, kdyy, kcurv);
}


void Grid::InverseStep(double inversionspeedmultiplier){
    
    GetSpeed(inversionspeedmultiplier);
    
    //NormUpwindGrad(kSurface,kSpeed,kNorm);
    //Advance(kSurface,kSpeed,kNorm,omega);
    
    Advance(kSurface,kSpeed,inversionspeedmultiplier);
    
}

void Grid::DoInversion(int maxiter, double inversionspeedmultiplier){
    
    //int maxiters = 100;
    double adjustspeed = inversionspeedmultiplier/maxiter/2.;
    
    for(int iters = 0; iters < maxiter; iters++){
        InverseStep(inversionspeedmultiplier);
        
        for(int i=0; i<maxi;i++)
            for(int j=0; j<maxj;j++)
                data_(i,j,kTemp1) = data_(i,j,kSurface);
        
        Reinitialize(kTemp1, kSurface, kitemp1, kitemp2);
        inversionspeedmultiplier = inversionspeedmultiplier-adjustspeed;
        double maxsurface = *std::max_element(&data_(0,0,kSurface),&data_(maxi,maxj,kSurface));
        if(maxsurface > (pow(dx*nX,2)+pow(dy*nY,2))){
            //std::cout << "Max surface is " << maxsurface << std::endl;
            std::cout << "Source has dissapeared. Try decreasing inversionspeedmultipler." << std::endl;
            exit(9);
        }
    
    }
    

 
    /*
    for(int iters = 0; iters < maxiters; iters++){
        InverseStep(omega);
        
    }
    for(int i=0; i<maxi;i++)
        for(int j=0; j<maxj;j++)
            data_(i,j,kTemp1) = data_(i,j,kSurface);
    
    Reinitialize(kTemp1, kSurface, kitemp1, kitemp2);
    
    for(int iters = 0; iters < maxiters; iters++){
        InverseStep(omega);
        
    }
    
    for(int i=0; i<maxi;i++)
        for(int j=0; j<maxj;j++)
            data_(i,j,kTemp1) = data_(i,j,kSurface);
    
    Reinitialize(kTemp1, kSurface, kitemp1, kitemp2);
*/

    
}


void Grid::OutputToFile(void){
    double* xout = new double[maxi];
    double* yout = new double[maxj];
    double* dataout = new double[maxi*maxj];
    double* sourceout = new double[maxi*maxj];
    double* forwardout = new double[maxi*maxj];
    double* surfout = new double[maxi*maxj];
    double* speed = new double[maxi*maxj];

    
    for(int i=0;i<maxi;i++)
        xout[i] = X(i);
    for(int j=0;j<maxj;j++)
        yout[j] = Y(j);
    
    for(int i=0;i<maxi;i++)
        for(int j=0;j<maxj;j++){
            surfout[cind(i,j)] = data_(i,j,kSurface);
            dataout[cind(i,j)] = data_(i,j,kData);
            sourceout[cind(i,j)] = data_(i,j,kSource);
            speed[cind(i,j)] = data_(i,j,kSpeed);
            forwardout[cind(i,j)] = data_(i,j,kForward);


        }

    std::cout << "Outputting Fields to File " << std::endl;
    FILE *file = fopen("data.bin", "wb");
    fwrite(&maxi,sizeof(int),1,file);
    fwrite(&maxj,sizeof(int),1,file);
    fwrite(xout,sizeof(double),maxi,file);
    fwrite(yout,sizeof(double),maxj,file);
    fwrite(dataout,sizeof(double),(maxi)*(maxj),file);
    fwrite(sourceout,sizeof(double),(maxi)*(maxj),file);
    fwrite(forwardout,sizeof(double),(maxi)*(maxj),file);
    fwrite(surfout,sizeof(double),(maxi)*(maxj),file);
    fwrite(speed,sizeof(double),(maxi)*(maxj),file);
    


    //fwrite(xout,sizeof(double),maxi*maxj,file);

    fclose(file);
    
    delete[] xout;
    delete[] yout;
    delete[] dataout;
    delete[] sourceout;
    delete[] forwardout;
    delete[] surfout;
    delete[] speed;

}
