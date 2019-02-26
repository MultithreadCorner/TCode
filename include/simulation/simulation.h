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
 * simulation.h
 *
 *  Created on: 13/01/2019
 *  Author: Andrea Contu
 */

#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include <simulation/Evolve.h>
#include <deposits/deposit.h>

using namespace hydra::placeholders;

namespace sim{
    //class implementing single simulation
    class simulation{
    private:
        size_t _bunchsize;
        size_t _nparticles;
        size_t _nsteps;
        std::string _path;
        std::string _nickname;
        bool _plot;
        bool _extrainfo;
        double _temperature;
        double _timestep;
        double _xshift;
        double _yshift;
        double _zshift;
        double _sig;
        double _amp;
        double _length;
//         RunningStateHost_t _data_dinit;
        bool _prepared;
        RunningState_t _data;
        ReducedDataDev_t _currents;
        CurrentVec_t _tp_currs;
        
    public:
        //constructors
        simulation(singleconf cf){
            _bunchsize      =   std::get<size_t>(cf.get("bunchsize"));
            _nparticles     =   std::get<size_t>(cf.get("nparticles"));
            _nsteps         =   std::get<size_t>(cf.get("nsteps"));
            _path           =   std::get<std::string>(cf.get("path"));
            _nickname       =   std::get<std::string>(cf.get("nickname"));
            _plot           =   std::get<bool>(cf.get("plot"));
            _extrainfo      =   std::get<bool>(cf.get("extrainfo"));
            _temperature    =   std::get<double>(cf.get("temperature"));
            _timestep       =   std::get<double>(cf.get("timestep"));
            _xshift         =   std::get<double>(cf.get("xshift"));
            _yshift         =   std::get<double>(cf.get("yshift"));
            _zshift         =   std::get<double>(cf.get("zshift"));
            _sig            =   std::get<double>(cf.get("sig"));
            _amp            =   std::get<double>(cf.get("amp"));
            _length         =   std::get<double>(cf.get("length"));
            _prepared       =   false;
        }
        
        
        //prepare simulation
        void prepare(){
            deposit::owndeposit mydep(_path,_nparticles,_length);
//             if(_nparticles==0 || _nparticles> mydep.getnparticles()) _nparticles = mydep.getnparticles();
// //             _data.resize(_nparticles);
//             hydra::copy(hydra::make_range(mydep.getdata().begin(),mydep.getdata().begin() + _nparticles),_data); // copy deposit to device
//             _prepared       = true;
//             hydra::for_each(_data, loaddata::Translate(xshift,yshift,zshift,_xmin,_xmax,_ymin,_ymax,_zmin,_zmax));
        }
        
    };
    
}

#endif
