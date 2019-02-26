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
 * loadfiles.h
 *
 *  Created on: 12/11/2018
 *      Author: Andrea Contu
 */


//Functions to load files

using namespace hydra::placeholders;

namespace loaddata{

    RunningStateHost_t getDummy(size_t nparticles, double length, size_t group=DEFAULT_GROUP){
        INFO_LINE("Generating deposit")
        hydra::Random<> Generator( std::chrono::system_clock::now().time_since_epoch().count());
        size_t np=(nparticles<=0) ? MAXPARTICLES : nparticles;

        if(np==0) ERROR_LINE("Number of particles is zero! Check your configuration.")
        if(np==1){
            if(group >1 && np/group==0) ERROR_LINE("Number of particles (or bunches of particles) is zero! Check your configuration.")
            auto data_d1=RunningStateHost_t(np,RunningTuple_t(std::copysign(1.,length)*group,0.,0.,0., 1.,0.,0.,0.,0.,0.,0.,0.));
            return data_d1;
        }
        size_t halfsize=np/2;
        if(group>1) halfsize/=group;
        RunningStateHost_t data_d(halfsize*2);
        auto data_de=RunningStateHost_t(halfsize,RunningTuple_t(-1.*group,0.,0.,0., 1.,0.,0.,0.,0.,0.,0.,0.));
        auto data_dh=RunningStateHost_t(halfsize,RunningTuple_t(1.*group,0.,0.,0., 1.,0.,0.,0.,0.,0.,0.,0.));
        
        
        //distribute them randomly in a line along Z
        if(length>0.){
            int seed=0;
            Generator.SetSeed(seed); // IMPORTANT, otherwise the seed stays the same
            Generator.Uniform(0., length, data_de.begin(_tc_z), data_de.end(_tc_z));
            Generator.SetSeed(seed);
            Generator.Uniform(0., length, data_dh.begin(_tc_z), data_dh.end(_tc_z));
        }
        hydra::copy(data_de,hydra::make_range(data_d.begin(),data_d.begin()+halfsize));
        hydra::copy(data_dh,hydra::make_range(data_d.begin()+halfsize,data_d.end()));
        
        return data_d;
    }


    RunningStateHost_t loadfile(std::string path, size_t group=DEFAULT_GROUP){
        std::ifstream infile;
        
        RunningStateHost_t data_d;
//         if(path=="dummy"){
//             return getDummy(MAXPARTICLES);
//         }
        
        infile.open(path.c_str());// file containing numbers in 3 columns 
        
        
        if(infile.fail()){ 
            ERROR_LINE("Error loading file " << path << ". Please check!")
            exit(1);
        }
    
        double pix=0;
        double piy=0;
        double piz=0;
        double pfx=0;
        double pfy=0;
        double pfz=0;
        double typ=0;
        
        std::string line;
        
        
        double deposit=0;
        int ct=0;
        while(std::getline(infile,line)){ 
            typ=0;
            std::stringstream lineStream;
            lineStream << line;
            lineStream >> pix >> piy >> piz >> pfx >> pfy >> pfz >> deposit >> typ;
    //         std::cout << line << std::endl;
            if(deposit<1) continue;
            
            unsigned int idep=deposit/group;
            if(idep<1) continue;
    //         std::cout << pix << "\t" << piy << "\t" << piz << "\t" << pfx << "\t" << pfy << "\t" << pfz << "\t" << deposit << std::endl;
            double d=sqrt((pfx-pix)*(pfx-pix)+(pfy-piy)*(pfy-piy)+(pfz-piz)*(pfz-piz));
            
            if(d==0){
                for(unsigned int i=0;i<idep;i++){
                    data_d.push_back(RunningTuple_t(-1.*group,pix, piy, piz, 1.,0.,0.,0.,0.,0.,0.,typ));
                    data_d.push_back(RunningTuple_t(1.*group,pix, piy, piz, 1.,0.,0.,0.,0.,0.,0.,typ));
                }
            }
            else{
                double a=(pfx-pix)/d;
                double b=(pfy-piy)/d;
                double c=(pfz-piz)/d;
                
                std::random_device rd;
                std::mt19937 mt(rd());
                std::uniform_real_distribution<double> dist(0.,d);
                for(unsigned int i=0;i<idep;i++){
                    double vec=dist(mt);
                    data_d.push_back(RunningTuple_t(-1.*group,a*vec+pix, b*vec+piy, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));
                    data_d.push_back(RunningTuple_t(1.*group,a*vec+pix, b*vec+piy, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));
                }
            }

            ct++;
        } 
        infile.close(); 
    //     std::cout << "Size: "<< data_d.size()<< std::endl;
        return data_d;
    }
    
    
    
    // functor to transform input deposit
    struct Translate
    {

        Translate()=delete;                                         // avoid creating objects with undefined data

        // No __hydra__dual__ here. Only allow objects to be constructed in host side. 
        Translate(double x, double y, double z, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax):
        fX(x),
        fY(y),
        fZ(z),
        fxmin(xmin),
        fxmax(xmax),
        fymin(ymin),
        fymax(ymax),
        fzmin(zmin),
        fzmax(zmax)
        {}

        // Copy constructor needs be __hydra_dual__ to be copyable
        // on device and host sides, in the time of thread distribution.
        // hydra, as well as thrust and std pass arguments by value.
        // Obs.: RooFit passes by reference... (crazy design choice!)
        __hydra_dual__
        Translate( Translate const& other):
        fX(other.fX),
        fY(other.fY),
        fZ(other.fZ),
        fxmin(other.fxmin),
        fxmax(other.fxmax),
        fymin(other.fymin),
        fymax(other.fymax),
        fzmin(other.fzmin),
        fzmax(other.fzmax)
        { }


        //function call operator needs to be callable in host and device sides 
        // and be constant
        template<typename Particle>
        __hydra_dual__
        void operator()(Particle p){
            
            hydra::get<_tc_x>(p)=hydra::get<_tc_x>(p)+fX;
            hydra::get<_tc_y>(p)=hydra::get<_tc_y>(p)+fY;
            hydra::get<_tc_z>(p)=hydra::get<_tc_z>(p)+fZ;
            
            if(hydra::get<_tc_x>(p)<fxmin || hydra::get<_tc_x>(p)>fxmax || hydra::get<_tc_y>(p)<fymin || hydra::get<_tc_y>(p)>fymax || hydra::get<_tc_z>(p)<fzmin || hydra::get<_tc_z>(p)>fzmax) hydra::get<_tc_isin>(p)=-1;
            
        }
        
        
    // private:
        double fX;
        double fY;
        double fZ;
        double fxmin;
        double fxmax;
        double fymin;
        double fymax;
        double fzmin;
        double fzmax;

    };
}
