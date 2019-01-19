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
 * datatypes.h
 *
 *  Created on: 12/11/2018
 *  Author: Andrea Contu
 */

#ifndef __DATATYPES_H__
#define __DATATYPES_H__

//defining particle states columns
#define RSDIM 12 // columns in running state
#define RSHOSTDIM 7 // columns in Universe state
#define RSINT 8 // entries in currents tuple

//particle states
using RunningStateHost_t    =   hydra::multiarray<double,RSDIM, hydra::host::sys_t >; //for testing
using RunningState_t        =   hydra::multiarray<double,RSDIM, hydra::device::sys_t >;
using RunningTuple_t        =   hydra::tuple<double,double,double,double,double,double,double,double,double,double,double,double>;
//default initialization
RunningTuple_t RunningTuple_init(-1, 0., 0., 0., 1.,0.,0.,0.,0.,0.,0.,0.);

using StateDev_t            =   hydra::multiarray<double,RSHOSTDIM,  hydra::device::sys_t>; //Particles state in device (e.g. GPU)
using StateHost_t           =   hydra::multiarray<double,RSHOSTDIM,  hydra::host::sys_t>; //Particles state in host
using StateTuple_t          =   hydra::tuple<double,double,double,double,double,double,double>;
StateTuple_t StateTuple_init(0.,0.,0.,0.,0.,0.,0.);


using UniverseDev_t         =   std::vector<StateDev_t>; //All particles in all states, it has to be in the host
using UniverseHost_t        =   std::vector<StateHost_t>; //All particles in all states, it has to be in the host
using ReducedDataDev_t      =   hydra::multiarray<double,RSINT,  hydra::device::sys_t>; //Reduced data container
using ReducedTuple_t        =   hydra::tuple<double,double,double,double,double,double,double,double>;
ReducedTuple_t ReducedTuple_init(0.,0.,0.,0.,0.,0.,0.,0.);

template<typename T>
using VecHost_t             =   hydra::host::vector<T>; //vector container double
template<typename T>
using VecDev_t              =   hydra::device::vector<T>; //vector container double

//currents vector
using CurrentVec_t          =   std::vector<ReducedTuple_t>;

//defining particle state columns meaning
auto _tc_charge     = hydra::placeholders::_0;
auto _tc_x          = hydra::placeholders::_1;
auto _tc_y          = hydra::placeholders::_2;
auto _tc_z          = hydra::placeholders::_3;
auto _tc_isin       = hydra::placeholders::_4;
auto _tc_curr       = hydra::placeholders::_5;
auto _tc_gauss_x    = hydra::placeholders::_6;
auto _tc_gauss_y    = hydra::placeholders::_7;
auto _tc_gauss_z    = hydra::placeholders::_8;
auto _tc_angle_1    = hydra::placeholders::_9;
auto _tc_angle_2    = hydra::placeholders::_10;
auto _tc_issec      = hydra::placeholders::_11;

#endif
