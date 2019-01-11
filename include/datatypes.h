/*----------------------------------------------------------------------------
 * 
 *   Copyright (C) 2018 Andrea Contu e Angelo Loi
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
 *      Author: Andrea Contu
 */

//defining particle states at a given time
#define _tc_charge _0
#define _tc_x _1
#define _tc_y _2
#define _tc_z _3
#define _tc_isin _4
#define _tc_curr _5
#define _tc_gauss_x _6
#define _tc_gauss_y _7
#define _tc_gauss_z _8
#define _tc_angle_1 _9
#define _tc_angle_2 _10
#define _tc_issec _11
#define RSDIM 12
#define RSHOSTDIM 7
#define RSINT 8

typedef hydra::multiarray<double,RSDIM, hydra::host::sys_t > RunningStateHost_t; //for testing
typedef hydra::multiarray<double,RSDIM, hydra::device::sys_t > RunningState_t;
typedef hydra::tuple<double,double,double,double,double,double,double,double,double,double,double,double> RunningTuple_t;
RunningTuple_t RunningTuple_init(-1, 0., 0., 0., 1.,0.,0.,0.,0.,0.,0.,0.);

typedef hydra::multiarray<double,RSHOSTDIM,  hydra::device::sys_t> StateDev_t; //Particles state in device (e.g. GPU)
typedef hydra::multiarray<double,RSHOSTDIM,  hydra::host::sys_t> StateHost_t; //Particles state in host
typedef hydra::tuple<double,double,double,double,double,double,double> StateTuple_t;
StateTuple_t StateTuple_init(0.,0.,0.,0.,0.,0.,0.);


typedef std::vector<StateDev_t> UniverseDev_t; //All particles in all states, it has to be in the host
typedef std::vector<StateHost_t> UniverseHost_t; //All particles in all states, it has to be in the host
typedef hydra::multiarray<double,RSINT,  hydra::device::sys_t> ReducedDataDev_t; //Reduced data container
typedef hydra::tuple<double,double,double,double,double,double,double,double> ReducedTuple_t;
ReducedTuple_t ReducedTuple_init(0.,0.,0.,0.,0.,0.,0.,0.);

typedef hydra::host::vector<double> VecHost_t; //vector container double
typedef hydra::device::vector<double> VecDev_t; //vector container double



typedef std::vector<ReducedTuple_t> CurrentVec_t;
typedef std::array<double, 11> FileLine_t;
typedef hydra::tuple<double,double,double> hydra_tuple3_t;
typedef hydra::tuple<double,double,double,double> hydra_tuple4_t;
typedef hydra::tuple<double,double,double,double,double,double> hydra_tuple6_t;
typedef hydra::tuple<double,double,double,double,double,double,double> hydra_tuple7_t;
typedef hydra::tuple<double,double,double,double,double,double,double,double> hydra_tuple8_t;
typedef hydra::tuple<double,double,double,double,double,double,double,double,double> hydra_tuple9_t;
typedef hydra::tuple<double,double,double,double,double,double,double,double,double,double> hydra_tuple10_t;
typedef hydra::tuple<double,double,double,double,double,double,double,double,double,double,double> hydra_tuple11_t;
typedef hydra::tuple<double,double,double,double,double,double,double,double,double,double,double,double> hydra_tuple12_t;

