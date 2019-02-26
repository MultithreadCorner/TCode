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
 *  deposit.h
 *
 *  Created on: 13/01/2019
 *  Author: Andrea Contu
 */

#ifndef __DEPOSIT_H__
#define __DEPOSIT_H__

namespace deposit{
    
    //class implementing deposit on host side
    
    class owndeposit{
    private:
        //static are useful not to reload deposits if not needed
        static std::string _static_path;
        static RunningStateHost_t _static_data;
        static size_t _static_nparticles;
        static double _static_length;
        
        std::string _path;
        RunningStateHost_t _data;
        size_t _nparticles;
        double _length=1.;
        
    public:
        owndeposit(std::string path, size_t nparticles, double length){
            this->init(path, nparticles, length);
        }
        
        owndeposit(singleconf cf){
            this->init(std::get<std::string>(cf.get("path")), std::get<size_t>(cf.get("nparticles")),std::get<double>(cf.get("length")));
        }
        //initialise
        void init(std::string path, size_t nparticles, double length){
            if(_path==DEFAULT_PATH){
                if(nparticles!=_static_nparticles || length!=_static_length){
                    _data = this->gendummy(nparticles,length);
                }
                else{
                    _data.resize(_nparticles);
                    hydra::copy(_static_data,_data);
                }
                _nparticles = nparticles;
                _length = length;
            }
            else{
                if(path!=_static_path){
                    _data = loaddata::loadfile(path);
                    _nparticles = _data.size();
                    _static_data.resize(_nparticles);
                    hydra::copy(_data,_static_data);
                }
                else{
                    _data.resize(_nparticles);
                    hydra::copy(_static_data,_data);
                    _nparticles = _data.size();
                }
            }
        }
        //get length
        double getlength(){return _length;}
        
        //get nparticles
        double getnparticles(){return _nparticles;}
        
        //generate dummy deposit
        RunningStateHost_t gendummy(size_t nparticles, double length){
            INFO_LINE("Generating deposit")
            hydra::Random<> Generator( std::chrono::system_clock::now().time_since_epoch().count());
            size_t np=(nparticles<=0) ? MAXPARTICLES : nparticles;
            if(np==1){
                auto data_d1=RunningStateHost_t(np,RunningTuple_t(std::copysign(1.,length),0.,0.,0., 0.,1.,0.,0.,0.,0.,0.,0.));
                return data_d1;
            }
            size_t halfsize=np/2;
            RunningStateHost_t data_d(halfsize*2);
            auto data_de=RunningStateHost_t(halfsize,RunningTuple_t(-1.,0.,0.,0., 1.,0.,0.,0.,0.,0.,0.,0.));
            auto data_dh=RunningStateHost_t(halfsize,RunningTuple_t(1.,0.,0.,0., 1.,0.,0.,0.,0.,0.,0.,0.));
            
            
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
        
        //get deposit
        RunningStateHost_t* getdata(){return &_data;}
        
        //get deposit and change if different
        RunningStateHost_t* getdata(std::string newfile, size_t nparticles, double length=1.){
            if(_path!=newfile){
                _path=newfile;
                if(newfile==DEFAULT_PATH) _data = loaddata::getDummy(nparticles,length);
                else _data=loaddata::loadfile(newfile);
            }
            return &_data;
        }
    };
}

#endif
