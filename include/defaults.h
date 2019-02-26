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
 * defaults.h
 *
 *  Created on: 12/11/2018
 *  Author: Andrea Contu
 */

#ifndef __DEFAULTS_H__
#define __DEFAULTS_H__

//project name
#define PROJECT_NAME "TCoDe"
#define VERSION "0.1-alpha"
#define DEFAULT_FILE_NAME "none"

//map code settings
#define DEFAULT_SC 0                    //default columns shift
#define DEFAULT_SR 0                    //default nuber of rows to ignore
#define DEFAULT_SCALE 1.0               //default scaling of map fields
#define DEFAULT_SCALESPACE 1.0          //default scaling of map space coordinates

//root plot settings
#define COLORELEC 4                     //electron color
#define COLORHOLE 2                     //hole color
#define COLORELECSEC 38                 //secondary electron color
#define COLORHOLESEC 42                 //secondary hole color

//simulation settings
#define MAXPARTICLES 10000              //max number of particles (if unspecified) in dummy deposit
#define MAXVDRIFTE 3.0e12               //max velocity cannot exceed speed of light (not used yet)
#define MAXVDRIFTH 3.0e12               //max velocity cannot exceed speed of light (not used yets)
#define UTIME 1e-12                     //time unit for plotting
#define DEFAULT_PARTICLES 0             //default number of particles
#define DEFAULT_STEPS 600.              //default number of time steps
#define DEFAULT_T 300.                  //default sensor temperature K
#define DEFAULT_TIME 1e-12              //default time step in seconds
#define DEFAULT_BUNCHSIZE 600           //bunch size (if gif animation or extra info requested). Timesteps are processed in bunches to avoid excessive memory usage
#define DEFAULT_AMP 0.                  //not used
#define DEFAULT_SIGMA 1.                //not used
#define DEFAULT_LENGTH 1.               //default deposit lenght for dummy deposit
#define DEFAULT_PATH "dummy"            //default deposit path (DO NOT CHANGE)
#define DEFAULT_DRAW false              //plot gif animation?
#define DEFAULT_EXTRAINFO false         //request full info dumping?
#define DEFAULT_X 0.                    //default deposit shift in x direction
#define DEFAULT_Y 0.                    //default deposit shift in y direction
#define DEFAULT_Z 0.                    //default deposit shift in z direction

#endif
