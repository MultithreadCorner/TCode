/*----------------------------------------------------------------------------
 * 
 *   Copyright (C) 2018-2019 Andrea Contu e Angelo Loi
 *
 *   This file is part of TCode software.
 *
 *   TCode is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   TCode is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with TCode.  If not, see <http://www.gnu.org/licenses/>.
 *
 *---------------------------------------------------------------------------*/
/*
 *  E_deposit.cpp
 *
 *  Created on: 18/02/2019
 *  Author: Angelo Loi
 */


//include gli header di ROOT necessari per la simulazione
#include "E_deposit.h" //classe field
#include <iostream>
#include <vector>
#include <fstream>
#include <TRandom.h>
#include <cmath>

#define PI 3.1415926535897932384626433

using namespace std;

E_deposit::E_deposit(){
}

void E_deposit::G4_E_Dep_gen(TString filename){
    //initialises particle number
    N_tot = 0;              
    shift = 0;
    //open G4 E_dep information file
    ifstream G4_file;
    G4_file.open(filename);
    
    if (G4_file.fail()){
        cout<<"FILE DOES NOT EXIST!\n";
        exit(EXIT_SUCCESS); ;
    }
    else
    {
        while(!G4_file.eof()){ 
            //generate deposit
            G4_file >> X0 >> Y0 >> Z0 >> X1 >> Y1 >> Z1 >> N_eh >> pId;
            
            for(int i=0+shift; i<N_eh+shift ; i++){
                t = (i-shift+1)/(N_eh);
                Pex.push_back(X0+(X1-X0)*t);
                Pey.push_back(Y0+(Y1-Y0)*t);
                Pez.push_back(Z0+(Z1-Z0)*t);
                Phx.push_back(X0+(X1-X0)*t);
                Phy.push_back(Y0+(Y1-Y0)*t);
                Phz.push_back(Z0+(Z1-Z0)*t);
                PIDe.push_back(pId);
                PIDh.push_back(pId);
                N_tot++;
            }
            //shift for next cycle 
            shift = N_tot;
        }
    }
    G4_file.close();
}



void E_deposit::SetCustomEDep(TString filename){

    N_tot = 0;                
    shift = 0;
    
    //open energy deposit file
    ifstream G4_file;
    G4_file.open(filename);
    
    if (G4_file.fail()){
        cout<<"FILE DOES NOT EXIST!\n";
        exit(EXIT_SUCCESS); ;
    }
    else{
        while(!G4_file.eof()){ 
            //file imput: first row is the same as for the G4 deposit    
            G4_file >> X0 >> Y0 >> Z0 >> X1 >> Y1 >> Z1 >> N_eh >> pId;
            //this row simply adds some other parameters
            G4_file >> T0 >> T1 >> sig0 >> sig1 >> Option;
            
            for(int i=0+shift; i<N_eh+shift ; i++){
                //linear evolution parameter
                t = (i-shift+1)/(N_eh);
                
                //evolution of the deposit is handled also lineary
                x_P = X0 + (X1 - X0)*t;
                y_P = Y0 + (Y1 - Y0)*t;
                z_P = Z0 + (Z1 - Z0)*t;
                
                //dispersion is handled lineary...
                sigm = sig0 + (sig1-sig0)*t;
                
                //the shape of the dispersion can be defined in 3 diff. shapes
                if(Option == 0){        //uniform
                    disp = gRandom->Uniform(-sigm/2, sigm/2);
                }
                else if(Option == 1){   //gaussian
                    disp = gRandom->Exp(sigm);
                }
                else if(Option == 2){   //exponential
                    disp = gRandom->Gaus(0,-sigm);
                }
                //... and the final position is still perpendicular to the 
                // direction of the deposit, calc. with the following function
                P_dir_disp();
                
                //now all electron-hole couples are generated with their casual
                //distribution arround the direction of the primary
                Pex.push_back(x_P + dispX);
                Pey.push_back(y_P + dispY);
                Pez.push_back(z_P + dispZ);
                Phx.push_back(x_P + dispX);
                Phy.push_back(y_P + dispY);
                Phz.push_back(z_P + dispZ);
                //their generation time is also included
                Te.push_back(T0+(T1-T0)*t);
                Th.push_back(T0+(T1-T0)*t);
                //the particle ID also
                PIDe.push_back(pId);
                PIDh.push_back(pId);
                N_tot++;
            }
            shift = N_tot;
        }
    }
    G4_file.close();
}

void E_deposit::P_dir_disp(){
    // calcs vector norm related to the deposit track
    mod = sqrt(((X1 - X0)*(X1 - X0) + (Y1 - Y0)*(Y1 - Y0) + (Z1 - Z0)*(Z1 - Z0)));

    // angle of the energy deposit track with respect to all axis
    alpha = acos(((Z1-Z0)/mod));  
    if(Y1 == Y0 && X1 == X0){
        beta = 0;
    }
    else{
        beta = atan(abs(Y1-Y0)/(X1-X0));                                                       
    }
        
    //definition of a perpendicular vector with respect to the energy deposit track
    //module is defined previously with the dispertion parameter sigm.    
    x_Y = disp*sin(alpha + PI/2)*cos(beta);
    y_Y = disp*sin(alpha + PI/2)*sin(beta);    
    z_Y = disp*cos(alpha + PI/2);
       
    
    //random angle generator for final rotation    
    gamma = (rand() % 36000 )*PI/18000;
    
    //casual rotation arround energy deposit track direction of the perpendicular vector 
    dispX = (((X1 - X0)*(X1 - X0)*(1/(mod*mod)) + (1-(X1 - X0)*(X1 - X0)*(1/(mod*mod)))*cos(gamma))*x_Y + ((X1 - X0)*(Y1 - Y0)*(1/(mod*mod))*(1-cos(gamma)) - (Z1 - Z0)*(1/(mod))*sin(gamma))*y_Y + ((X1 - X0)*(Z1 - Z0)*(1/(mod*mod))*(1 - cos(gamma)) + (Y1 - Y0)*(1/(mod))*sin(gamma))*z_Y);
    
    dispY = ((X1 - X0)*(1/(mod*mod))*(Y1 - Y0)*(1 - cos(gamma))*x_Y + (Z1 - Z0)*(1/(mod))*sin(gamma)*x_Y + (Y1 - Y0)*(1/(mod*mod))*(Y1 - Y0)*y_Y + (1-(Y1 - Y0)*(1/(mod*mod))*(Y1 - Y0))*cos(gamma)*y_Y + (Y1 - Y0)*(1/(mod*mod))*(Z1 - Z0)*(1 - cos(gamma))*z_Y - (X1 - X0)*(1/(mod))*sin(gamma)*z_Y);
    
    dispZ = ((X1 - X0)*(Z1 - Z0)*(1/(mod*mod))*(1-cos(gamma))*x_Y - (Y1 - Y0)*(1/(mod))*sin(gamma)*x_Y + (Y1 - Y0)*(1/(mod*mod))*(Z1 - Z0)*(1-cos(gamma))*y_Y + (X1 - X0)*(1/(mod*mod))*sin(gamma)*y_Y + (Z1 - Z0)*(Z1 - Z0)*(1/(mod*mod))*z_Y + (1-(Z1 - Z0)*(1/(mod*mod))*(Z1 - Z0))*cos(gamma)*z_Y);
   
}


//returns the total number of generated couples 
float E_deposit::Return_N_tot(){
    return N_tot;
}
//return all particles coordinates
float E_deposit::GethXPos(int i){
   return Phx[i];
}
float E_deposit::GethYPos(int i){
   return Phy[i];
}
float E_deposit::GethZPos(int i){
   return Phz[i];
}
float E_deposit::GeteXPos(int i){
   return Pex[i];
}
float E_deposit::GeteYPos(int i){
   return Pey[i];
}
float E_deposit::GeteZPos(int i){
   return Pez[i];
}


//simply generates an output to controll all generated particles
void* E_deposit::Print_E_deposit(){
    ofstream myfile;
    myfile.open ("deposit_at_t_zero.txt");
    for(int i=1 ; i<Pex.size() ; i++){
        cout<<Pex[i]<<" "<<Pey[i]<<" "<<Pez[i]<<" "<<PIDe[i]<<" "<<Phx[i]<<" "<<Phy[i]<<" "<<Phz[i]<<" "<<PIDh[i]<<"\n";
        myfile<<Pex[i]<<" "<<Pey[i]<<" "<<Pez[i]<<" "<<PIDe[i]<<" "<<Phx[i]<<" "<<Phy[i]<<" "<<Phz[i]<<" "<<PIDh[i]<<"\n";
    }
    myfile.close();
    cout<<N_tot<<" particles generated\n";
}
 