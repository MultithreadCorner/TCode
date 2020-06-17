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
 * main.inl
 *
 *  Created on: 12/11/2018
 *      Author: Andrea Contu
 */

#ifndef MAIN_INL_
#define MAIN_INL_

#include <details/constants.h>
#include <details/defaults.h>
#include <details/utils.h>
#include <details/datatypes.h>
#include <plotstyle.C>
#include <loadfiles.h>
#include <maps/preparemaps.h>
#include <analysis.h>
//facilities
#include <simulation/HydraSim.h>

int main(int argv, char** argc){
    
    printHeader();
    INFO_LINE("Running sensor response simulation")
    
    cfg::Config Cfg;
    
    std::string cfg_file_name = "settings.cfg";
    std::string EFfile=DEFAULT_FILE_NAME,WFfile=DEFAULT_FILE_NAME,EMOBfile=DEFAULT_FILE_NAME,HMOBfile=DEFAULT_FILE_NAME,PHYSMAPfile=DEFAULT_FILE_NAME;
    int icfg=-1;
    std::string outputdir="output";
    bool isconfig=false;
    bool draw=DEFAULT_DRAW; 
    bool extrainfo=DEFAULT_EXTRAINFO; 
    size_t nparticles=DEFAULT_PARTICLES;
    size_t group=DEFAULT_GROUP;
    size_t nsteps=DEFAULT_STEPS;
    double temperature=DEFAULT_T;
    double timestep=DEFAULT_TIME;
    int bunchsize=DEFAULT_BUNCHSIZE; //depends on machine, it should be as large as possible
    double ampexp = DEFAULT_AMP;
    double sigmaexp = DEFAULT_SIGMA;
    double xshift = DEFAULT_X;
    double yshift = DEFAULT_Y;
    double zshift = DEFAULT_Z;
    double length = DEFAULT_LENGTH;
    std::string file_single=DEFAULT_PATH;
    std::string file_single_name="InlineSim";
    //-------------------------
    // command line arguments
    //-------------------------
    try {
        
        TCLAP::CmdLine cmd("Command line arguments for ");
        
        //get configuration file name
        TCLAP::ValueArg<std::string> Cfg_Arg("c", "config","Configuration file", false, "settings.cfg", "string");
        cmd.add(Cfg_Arg);
        TCLAP::ValueArg<int> ICfg_Arg("i", "iconfig","Choose only one simulation from config file", false, -1 , "int");
        cmd.add(ICfg_Arg);
        //physics maps
        TCLAP::ValueArg<std::string> PMap_Arg("r", "physmap","Physics map", false, PHYSMAPfile.c_str(), "string");
        cmd.add(PMap_Arg);
        TCLAP::ValueArg<std::string> EFfile_Arg("", "efield","Electric field map", false, EFfile.c_str(), "string");
        cmd.add(EFfile_Arg);
        TCLAP::ValueArg<std::string> EMOBfile_Arg("", "emob","Electron mobility map", false, EMOBfile.c_str(), "string");
        cmd.add(EMOBfile_Arg);
        TCLAP::ValueArg<std::string> HMOBfile_Arg("", "hmob","Hole mobility map", false, HMOBfile.c_str(), "string");
        cmd.add(HMOBfile_Arg);
        TCLAP::ValueArg<std::string> WFfile_Arg("", "wfield","Weighting field map", false, WFfile.c_str(), "string");
        cmd.add(WFfile_Arg);
        
        //simulation parameters
        TCLAP::ValueArg<std::string> Fsingle_Arg("f", "path","Input file", false, DEFAULT_PATH, "string");
        cmd.add(Fsingle_Arg);
        TCLAP::ValueArg<std::string> Nsingle_Arg("m", "nickname","Name file", false, "InlineSim", "string");
        cmd.add(Nsingle_Arg);
        TCLAP::ValueArg<std::string> Odir_Arg("o", "OutputDirectory","Output directory", false, "output", "string");
        cmd.add(Odir_Arg);
        
        TCLAP::SwitchArg Plot_Arg("d", "plot","Produce animated gif", DEFAULT_DRAW);
        cmd.add(Plot_Arg);
        TCLAP::SwitchArg Extrainfo_Arg("e", "extrainfo","Store full information", DEFAULT_EXTRAINFO);
        cmd.add(Extrainfo_Arg);
        TCLAP::ValueArg<size_t> Nparts_Arg("p", "nparticles","Number of particles", false, DEFAULT_PARTICLES, "size_t");
        cmd.add(Nparts_Arg);
        TCLAP::ValueArg<size_t> Nsteps_Arg("n", "nsteps","Number of steps", false, DEFAULT_STEPS, "size_t");
        cmd.add(Nsteps_Arg);
        TCLAP::ValueArg<size_t> Sbunch_Arg("b", "bunchsize","Step bunch size", false, DEFAULT_BUNCHSIZE, "size_t");
        cmd.add(Sbunch_Arg);
        TCLAP::ValueArg<double> Temp_Arg("T", "temperature","Sensor temperature in Kelvin", false, DEFAULT_T, "double");
        cmd.add(Temp_Arg);
        TCLAP::ValueArg<double> Timestep_Arg("t", "timestep","Time step in ps", false, DEFAULT_TIME, "double");
        cmd.add(Timestep_Arg);
        TCLAP::ValueArg<double> Aexpl_Arg("a", "Ampexplosion","Explosion amplitude", false, DEFAULT_AMP, "double");
        cmd.add(Aexpl_Arg);
        TCLAP::ValueArg<double> Sexpl_Arg("s", "Sigmaexplosion","Explosion sigma", false, DEFAULT_SIGMA, "double");
        cmd.add(Sexpl_Arg);
        
        TCLAP::ValueArg<double> Xshift_Arg("x", "Xshift","X shift or entry point", false, DEFAULT_X, "double");
        cmd.add(Xshift_Arg);
        TCLAP::ValueArg<double> Yshift_Arg("y", "Yshift","Y shift or entry point", false, DEFAULT_Y, "double");
        cmd.add(Yshift_Arg);
        TCLAP::ValueArg<double> Zshift_Arg("z", "Zshift","Z shift or entry point", false, DEFAULT_Z, "double");
        cmd.add(Zshift_Arg);
        TCLAP::ValueArg<double> Length_Arg("l", "Length","Deposit length", false, DEFAULT_LENGTH, "double");
        cmd.add(Length_Arg);
        TCLAP::ValueArg<size_t> Ngroup_Arg("g", "group","Group particle", false, DEFAULT_GROUP, "size_t");
        cmd.add(Ngroup_Arg);
        
        
        // Parse the argv array.
        cmd.parse(argv, argc);
        
        // Get the value parsed by tclap.
        if(Cfg_Arg.isSet()) isconfig=true;
        cfg_file_name = Cfg_Arg.getValue();
        icfg = ICfg_Arg.getValue();
        file_single = Fsingle_Arg.getValue();
        file_single_name = Nsingle_Arg.getValue();
        draw = Plot_Arg.getValue();
        extrainfo = Extrainfo_Arg.getValue();
        nparticles = Nparts_Arg.getValue();
        nsteps = Nsteps_Arg.getValue();
        temperature = Temp_Arg.getValue();
        timestep = Timestep_Arg.getValue();
        bunchsize = Sbunch_Arg.getValue();
        ampexp = Aexpl_Arg.getValue();
        sigmaexp = Sexpl_Arg.getValue();
        xshift = Xshift_Arg.getValue();
        yshift = Yshift_Arg.getValue();
        zshift = Zshift_Arg.getValue();
        group  = Ngroup_Arg.getValue();
        length = Length_Arg.getValue();
        outputdir = Odir_Arg.getValue();
        PHYSMAPfile = PMap_Arg.getValue();
        
        EFfile = EFfile_Arg.getValue();
        WFfile = WFfile_Arg.getValue(); 
        EMOBfile = EMOBfile_Arg.getValue();
        HMOBfile = HMOBfile_Arg.getValue();
        
    }
    catch (TCLAP::ArgException &e)  {
        
        std::cerr << "Error : "  << e.error()
        << " for arg " << e.argId()
        << std::endl;
        
        return(EXIT_FAILURE);
    }
    
    //-------------------------
    // configuration file
    //-------------------------
    
    //if config is provided, load from config, else load from command line
    if(isconfig){
        //configuration representation
        std::cout<< TOPROW<< std::endl;
        INFO_LINE("  Reading settings from configuration file: "<<cfg_file_name)
        std::cout<< MIDDLEROW<< std::endl;
        
        // Read the file. If there is an error, report it and exit.
        try
        {
            Cfg.readFile(cfg_file_name.c_str());
        }
        catch(const cfg::FileIOException &fioex)
        {
            std::cerr << "I/O error while reading configuration file: "
            << cfg_file_name
            << std::endl;
            
            return(EXIT_FAILURE);
        }
        catch(const cfg::ParseException &pex)
        {
            std::cerr << "Parse error at " << pex.getFile() << " : Line  "
            << pex.getLine() << " - "
            << pex.getError() << std::endl;
            return(EXIT_FAILURE);
        }
        
        //-------------------------
        // actual work
        //-------------------------
        
        cfg::Setting &root = Cfg.getRoot();
        
        
        INFO_LINE("Configuration is valid, ignoring other options (if provided)")
        
        signal_simulation(root, icfg);
    }
    else{
        
        INFO_LINE("Loading inline commands")
        cfg::Setting & root = Cfg.getRoot();
        root.add("OutputDirectory", cfg::Setting::TypeString) = outputdir;
        
        cfg::Setting & root_line_maps = root["PhysicsMaps"];
        if(PHYSMAPfile!="none") root_line_maps.add("map",cfg::Setting::TypeString) = PHYSMAPfile;
        if(EFfile!="none") root_line_maps.add("efield", cfg::Setting::TypeString) = EFfile;
        if(WFfile!="none") root_line_maps.add("wfield", cfg::Setting::TypeString) = WFfile;
        if(EMOBfile!="none") root_line_maps.add("emob", cfg::Setting::TypeString) = EMOBfile;
        if(HMOBfile!="none") root_line_maps.add("hmob", cfg::Setting::TypeString) = HMOBfile;
        
        root.add("InputData",cfg::Setting::TypeGroup);
        cfg::Setting & root_line_simul = root["InputData"];
        root_line_simul.add("InlineSimulation",cfg::Setting::TypeGroup);
        cfg::Setting & root_line_inl = root_line_simul["InlineSimulation"];
        root_line_inl.add("path",cfg::Setting::TypeString) = file_single;
        root_line_inl.add("name",cfg::Setting::TypeString) = file_single_name;
        root_line_inl.add("plot",cfg::Setting::TypeBoolean) = draw;
        root_line_inl.add("process",cfg::Setting::TypeBoolean) = true;
        root_line_inl.add("extrainfo",cfg::Setting::TypeBoolean) = extrainfo;
        root_line_inl.add("particles",cfg::Setting::TypeInt) = (int)nparticles;
        root_line_inl.add("steps",cfg::Setting::TypeInt) = (int)nsteps;
        root_line_inl.add("timestep",cfg::Setting::TypeFloat) = (float)timestep;
        root_line_inl.add("T",cfg::Setting::TypeFloat) = (float)temperature;
        root_line_inl.add("Aexpl",cfg::Setting::TypeFloat) = (float)ampexp;
        root_line_inl.add("Sexpl",cfg::Setting::TypeFloat) = (float)sigmaexp;
        root_line_inl.add("Xshift",cfg::Setting::TypeFloat) = (float)xshift;
        root_line_inl.add("Yshift",cfg::Setting::TypeFloat) = (float)yshift;
        root_line_inl.add("Zshift",cfg::Setting::TypeFloat) = (float)zshift;
        root_line_inl.add("Length",cfg::Setting::TypeFloat) = (float)length;
        root_line_inl.add("BunchSize",cfg::Setting::TypeInt) = bunchsize;
        root_line_inl.add("group",cfg::Setting::TypeInt) = (int)group;
        
        
        
        signal_simulation(root, icfg);

    }
    
    return 0;
}





#endif /* MAIN_INL_ */
