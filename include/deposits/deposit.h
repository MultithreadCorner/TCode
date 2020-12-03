/*----------------------------------------------------------------------------
 * 
 *   Copyright (C) 2018-2020 Andrea Contu e Angelo Loi
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
        std::string _cached_path;
        RunningStateHost_t _cached_data;
        size_t _cached_nparticles;
        size_t _cached_group;
        double _cached_length;
        
        std::string _path;
        RunningStateHost_t _data;
        size_t _nparticles;
        size_t _group;
        double _length=1.;
        bool _isloaded;
        
    public:
        
        owndeposit(){}
        
        owndeposit(std::string path, size_t nparticles, double length, size_t group=DEFAULT_GROUP){
            this->init(path, nparticles, length);
        }
        
        owndeposit(singleconf& cf){
            this->init(cf);
        }
        //initialise
        
        void init(singleconf& cf){
            this->init(std::get<std::string>(cf.get("path")), std::get<size_t>(cf.get("nparticles")),std::get<double>(cf.get("length")),std::get<size_t>(cf.get("group")));
        }
        
        void init(std::string path, size_t nparticles, double length, size_t group=DEFAULT_GROUP){
            if(_path==DEFAULT_PATH){
                if(nparticles!=_cached_nparticles || length!=_cached_length || group!=_cached_group){
                    this->gendummy(nparticles,length,group);
                }
                else{
                    _data.resize(_nparticles);
                    hydra::copy(_cached_data,_data);
                }
                _nparticles = nparticles;
                _length = length;
            }
            else{
                if(path!=_cached_path){
                    _data = loaddata::loadfile(path);
                    _nparticles = _data.size();
                    _cached_data.resize(_nparticles);
                    hydra::copy(_data,_cached_data);
                    _cached_nparticles = _nparticles;
                    _cached_length = length;
                }
                else{
                    _data.resize(_nparticles);
                    hydra::copy(_cached_data,_data);
                    _nparticles = _data.size();
                }
            }
            _isloaded=true;
        }
        
        //get _isloaded
        bool isloaded(){return _isloaded;}
        
        //get length
        double getlength(){return _length;}
        
        //get nparticles
        double getnparticles(){return _nparticles;}
        
        //generate dummy deposit
        void gendummy(size_t nparticles, double length, size_t group=DEFAULT_GROUP){
            INFO_LINE("Generating deposit")
            auto A = hydra::Parameter::Create().Name("A").Value(0.0);
            auto B = hydra::Parameter::Create().Name("B").Value(length);
            hydra_thrust::default_random_engine engine;
            auto uniform   = hydra::UniformShape<double>(A,B);
            size_t np=(nparticles<=0) ? MAXPARTICLES : nparticles;
            if(np==1){
                _data=RunningStateHost_t(np,RunningTuple_t(std::copysign(1.,length),0.,0.,0., 0.,1.,0.,0.,0.,0.,0.,0.));
            }
            size_t halfsize=np/2;
            _data.resize(halfsize*2);
            auto data_de=RunningStateHost_t(halfsize,RunningTuple_t(-1.,0.,0.,0., 1.,0.,0.,0.,0.,0.,0.,0.));
            auto data_dh=RunningStateHost_t(halfsize,RunningTuple_t(1.,0.,0.,0., 1.,0.,0.,0.,0.,0.,0.,0.));
            
            
            //distribute them randomly in a line along Z
            if(length>0.){
                int seed=0;
                hydra::fill_random(data_de.begin(_tc_z), data_de.end(_tc_z) , uniform,seed);
                hydra::fill_random(data_dh.begin(_tc_z), data_dh.end(_tc_z) , uniform,seed);
            }
            hydra::copy(data_de,hydra::make_range(_data.begin(),_data.begin()+halfsize));
            hydra::copy(data_dh,hydra::make_range(_data.begin()+halfsize,_data.end()));
            
        }
        
        //get deposit
        RunningStateHost_t getdata(){return _data;}
        
        //get deposit and change if different
//         void SendToDevice(RunningStateDev_t& todev, singleconf)
// 
//         }
    };
}

#endif
