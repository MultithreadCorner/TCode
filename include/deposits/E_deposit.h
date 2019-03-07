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
 *  E_deposit.h
 *
 *  Created on: 18/02/2019
 *  Author: Angelo Loi
 */

//include gli header di root per gestione file e istogrammi 4D
#include "TH2F.h"
#include "TFile.h"
#include <vector>
#include <iostream>
class E_deposit{
//dichiarazione delle variabili/oggetti pubbliche    
    
protected:

    float X0, X1;
    float Y0, Y1;
    float Z0, Z1;
    
    float T0, T1;
    
    float alpha, beta, gamma;
    
    float sig0, sig1;
    
    float N_eh;
        
    float pId, prim, sec;
    
    float shift,t;
    
    float disp, dispX, dispY, dispZ;
    float versX, versY, versZ;
    float sigm, mod;
    
    int Option;
    
public:
    
    E_deposit(/*TString filename*/);
    
    void G4_E_Dep_gen(TString filename);
    
    float Return_N_tot();
    
    float GeteXPos(int i);
    float GeteYPos(int i);
    float GeteZPos(int i);
    float GethXPos(int i);
    float GethYPos(int i);
    float GethZPos(int i);
    
    float GeteTime(int i);
    float GethTime(int i);

    float GetePID(int i);
    float GethPID(int i);
    
    void* Print_E_deposit();
    
    void P_dir_disp();
    
    void SetCustomEDep(TString filename);
    
    std::vector<float> Pex;
    std::vector<float> Pey;
    std::vector<float> Pez;
    std::vector<float> Phx;
    std::vector<float> Phy;
    std::vector<float> Phz;
    
    std::vector<float> Te ;
    std::vector<float> Th ;
    
    std::vector<float> PIDe;
    std::vector<float> PIDh;
    
    float N_tot;
    
    float x_Y, z_Y , y_Y, x_F, y_F, z_F;
    
    float x_P, y_P, z_P;
    
};
