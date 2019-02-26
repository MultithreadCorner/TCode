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
 * main_maps.inl
 *
 *  Created on: 12/11/2018
 *      Author: Andrea Contu
 */

#ifndef MAIN_MAPS_INL_
#define MAIN_MAPS_INL_

#include <details/constants.h>
#include <details/defaults.h>
#include <details/utils.h>
#include <details/datatypes.h>
//facilities
#include <maps/preparemaps.h>

int main(int argv, char** argc){
    
    printHeader();
    INFO_LINE("Running physics map preparation")
    
    cfg::Config Cfg;
    
    
    std::string EFfile=DEFAULT_FILE_NAME,WFfile=DEFAULT_FILE_NAME,EMOBfile=DEFAULT_FILE_NAME,HMOBfile=DEFAULT_FILE_NAME,PHYSMAPfile=DEFAULT_FILE_NAME;
    std::string EFfile_out = (std::string)PROJECT_NAME+(std::string)"_efield.dat";
    std::string WFfile_out = (std::string)PROJECT_NAME+(std::string)"_wfield.dat";
    std::string EMOBfile_out = (std::string)PROJECT_NAME+(std::string)"_emob.dat";
    std::string HMOBfile_out = (std::string)PROJECT_NAME+(std::string)"_hmob.dat";
    std::string PHYSMAPfile_out = (std::string)PROJECT_NAME+(std::string)"_physmap.dat";
    std::string outputdir="Maps";
    int efield_sc=DEFAULT_SC, wfield_sc=DEFAULT_SC, emob_sc=DEFAULT_SC, hmob_sc=DEFAULT_SC, physmap_sc=DEFAULT_SC;
    int efield_sr=DEFAULT_SR, wfield_sr=DEFAULT_SR, emob_sr=DEFAULT_SR, hmob_sr=DEFAULT_SR, physmap_sr=DEFAULT_SR;
    double efield_scalespace=DEFAULT_SCALESPACE, wfield_scalespace=DEFAULT_SCALESPACE, emob_scalespace=DEFAULT_SCALESPACE, hmob_scalespace=DEFAULT_SCALESPACE, physmap_scalespace=DEFAULT_SCALESPACE;
    double efield_scale=DEFAULT_SCALE, wfield_scale=DEFAULT_SCALE, emob_scale=DEFAULT_SCALE, hmob_scale=DEFAULT_SCALE, physmap_scale=DEFAULT_SCALE;

    //-------------------------
    // command line arguments
    //-------------------------
    try {
        
        TCLAP::CmdLine cmd("Command line arguments for ");
        
        //settings to load maps
        TCLAP::ValueArg<std::string> EFfile_Arg("e", "efield","Electric field map", false, EFfile.c_str(), "string");
        cmd.add(EFfile_Arg);
        TCLAP::ValueArg<int> EFfile_sc_Arg("", "efield_sc","Electric field map skip columns", false, efield_sc, "int");
        cmd.add(EFfile_sc_Arg);
        TCLAP::ValueArg<int> EFfile_sr_Arg("", "efield_sr","Electric field map skip rows", false, efield_sr, "int");
        cmd.add(EFfile_sr_Arg);
        TCLAP::ValueArg<double> EFfile_scale_Arg("", "efield_scale","Electric field map scaling factor", false, efield_scale, "double");
        cmd.add(EFfile_scale_Arg);
        TCLAP::ValueArg<double> EFfile_scalespace_Arg("", "efield_scalespace","Electric field map space scaling factor", false, efield_scalespace, "double");
        cmd.add(EFfile_scalespace_Arg);
        TCLAP::ValueArg<std::string> EFfile_output_Arg("", "efield_out","Electric field map output file", false, EFfile_out.c_str(), "string");
        cmd.add(EFfile_output_Arg);
        
        TCLAP::ValueArg<std::string> EMOBfile_Arg("E", "emob","Electron mobility map", false, EMOBfile.c_str(), "string");
        cmd.add(EMOBfile_Arg);
        TCLAP::ValueArg<int> EMOBfile_sc_Arg("", "emob_sc","Electron mobility map skip columns", false, emob_sc, "int");
        cmd.add(EMOBfile_sc_Arg);
        TCLAP::ValueArg<int> EMOBfile_sr_Arg("", "emob_sr","Electron mobility map skip rows", false, emob_sr, "int");
        cmd.add(EMOBfile_sr_Arg);
        TCLAP::ValueArg<double> EMOBfile_scale_Arg("", "emob_scale","Electron mobility map scaling factor", false, emob_scale, "double");
        cmd.add(EMOBfile_scale_Arg);
        TCLAP::ValueArg<double> EMOBfile_scalespace_Arg("", "emob_scalespace","Electron mobility map space scaling factor", false, emob_scalespace, "double");
        cmd.add(EMOBfile_scalespace_Arg);
        TCLAP::ValueArg<std::string> EMOBfile_output_Arg("", "emob_out","Electron mobility map output file", false, EMOBfile_out.c_str(), "string");
        cmd.add(EMOBfile_output_Arg);
        
        TCLAP::ValueArg<std::string> HMOBfile_Arg("H", "hmob","Hole mobility map", false, HMOBfile.c_str(), "string");
        cmd.add(HMOBfile_Arg);
        TCLAP::ValueArg<int> HMOBfile_sc_Arg("", "hmob_sc","Hole mobility map skip columns", false, hmob_sc, "int");
        cmd.add(HMOBfile_sc_Arg);
        TCLAP::ValueArg<int> HMOBfile_sr_Arg("", "hmob_sr","Hole mobility map skip rows", false, hmob_sr, "int");
        cmd.add(HMOBfile_sr_Arg);
        TCLAP::ValueArg<double> HMOBfile_scale_Arg("", "hmob_scale","Hole mobility map scaling factor", false, hmob_scale, "double");
        cmd.add(HMOBfile_scale_Arg);
        TCLAP::ValueArg<double> HMOBfile_scalespace_Arg("", "hmob_scalespace","Hole mobility map space scaling factor", false, hmob_scalespace, "double");
        cmd.add(HMOBfile_scalespace_Arg);
        TCLAP::ValueArg<std::string> HMOBfile_output_Arg("", "hmob_out","Hole mobility map output file", false, HMOBfile_out.c_str(), "string");
        cmd.add(HMOBfile_output_Arg);
        
        TCLAP::ValueArg<std::string> WFfile_Arg("w", "wfield","Weighting field map", false, WFfile.c_str(), "string");
        cmd.add(WFfile_Arg);
        TCLAP::ValueArg<int> WFfile_sc_Arg("", "wfield_sc","Weighting field map skip columns", false, wfield_sc, "int");
        cmd.add(WFfile_sc_Arg);
        TCLAP::ValueArg<int> WFfile_sr_Arg("", "wfield_sr","Weighting field map skip rows", false, wfield_sr, "int");
        cmd.add(WFfile_sr_Arg);
        TCLAP::ValueArg<double> WFfile_scale_Arg("", "wfield_scale","Weighting field map scaling factor", false, wfield_scale, "double");
        cmd.add(WFfile_scale_Arg);
        TCLAP::ValueArg<double> WFfile_scalespace_Arg("", "wfield_scalespace","Weighting field map space scaling factor", false, wfield_scalespace, "double");
        cmd.add(WFfile_scalespace_Arg);
        TCLAP::ValueArg<std::string> WFfile_output_Arg("", "wfield_out","Weighting field map output file", false, WFfile_out.c_str(), "string");
        cmd.add(WFfile_output_Arg);
        
        TCLAP::ValueArg<std::string> PHYSMAPfile_Arg("p", "physmap","Full map of fields", false, PHYSMAPfile.c_str(), "string");
        cmd.add(PHYSMAPfile_Arg);
        TCLAP::ValueArg<int> PHYSMAPfile_sc_Arg("", "physmap_sc","Full map skip columns", false, physmap_sc, "int");
        cmd.add(PHYSMAPfile_sc_Arg);
        TCLAP::ValueArg<int> PHYSMAPfile_sr_Arg("", "physmap_sr","Full map skip rows", false, physmap_sr, "int");
        cmd.add(PHYSMAPfile_sr_Arg);
        TCLAP::ValueArg<double> PHYSMAPfile_scale_Arg("", "physmap_scale","Full map scaling factor", false, physmap_scale, "double");
        cmd.add(PHYSMAPfile_scale_Arg);
        TCLAP::ValueArg<double> PHYSMAPfile_scalespace_Arg("", "physmap_scalespace","Full map space scaling factor", false, physmap_scalespace, "double");
        cmd.add(PHYSMAPfile_scalespace_Arg);
        TCLAP::ValueArg<std::string> PHYSMAPfile_output_Arg("", "physmap_out","Full map of fields output file", false, PHYSMAPfile_out.c_str(), "string");
        cmd.add(PHYSMAPfile_output_Arg);
        
        TCLAP::ValueArg<std::string> Odir_Arg("o", "OutputDirectory","Output directory", false, outputdir, "string");
        cmd.add(Odir_Arg);
        
        // Parse the argv array.
        cmd.parse(argv, argc);
        
        // Get the value parsed by tclap.
        outputdir = Odir_Arg.getValue();
        
        EFfile = EFfile_Arg.getValue();
        EFfile_out = EFfile_output_Arg.getValue();
        efield_sc = EFfile_sc_Arg.getValue();
        efield_sr = EFfile_sr_Arg.getValue();
        efield_scale = EFfile_scale_Arg.getValue();
        efield_scalespace = EFfile_scalespace_Arg.getValue();
        
        WFfile = WFfile_Arg.getValue();
        WFfile_out = WFfile_output_Arg.getValue();
        wfield_sc = WFfile_sc_Arg.getValue();
        wfield_sr = WFfile_sr_Arg.getValue();
        wfield_scale = WFfile_scale_Arg.getValue();
        wfield_scalespace = WFfile_scalespace_Arg.getValue();
        
        EMOBfile = EMOBfile_Arg.getValue();
        EMOBfile_out = EMOBfile_output_Arg.getValue();
        emob_sc = EMOBfile_sc_Arg.getValue();
        emob_sr = EMOBfile_sr_Arg.getValue();
        emob_scale = EMOBfile_scale_Arg.getValue();
        emob_scalespace = EMOBfile_scalespace_Arg.getValue();
        
        HMOBfile = HMOBfile_Arg.getValue();
        HMOBfile_out = HMOBfile_output_Arg.getValue();
        hmob_sc = HMOBfile_sc_Arg.getValue();
        hmob_sr = HMOBfile_sr_Arg.getValue();
        hmob_scale = HMOBfile_scale_Arg.getValue();
        hmob_scalespace = HMOBfile_scalespace_Arg.getValue();
        
        PHYSMAPfile = PHYSMAPfile_Arg.getValue();
        PHYSMAPfile_out = PHYSMAPfile_output_Arg.getValue();
        physmap_sc = PHYSMAPfile_sc_Arg.getValue();
        physmap_sr = PHYSMAPfile_sr_Arg.getValue();
        physmap_scale = PHYSMAPfile_scale_Arg.getValue();
        physmap_scalespace = PHYSMAPfile_scalespace_Arg.getValue();
    }
    catch (TCLAP::ArgException &e)  {
        
        std::cerr << "Error : "  << e.error()
        << " for arg " << e.argId()
        << std::endl;
        
        return(EXIT_FAILURE);
    }
    

    INFO_LINE("Loading inline commands")
    cfg::Setting & root = Cfg.getRoot();
    
    root.add("OutputDirectory", cfg::Setting::TypeString) = outputdir;
    if(PHYSMAPfile!="none"){
        root.add("map",cfg::Setting::TypeGroup);
        cfg::Setting & root_pm = root["map"];
        root_pm.add("path", cfg::Setting::TypeString) = PHYSMAPfile;
        root_pm.add("out", cfg::Setting::TypeString) = PHYSMAPfile_out;
        root_pm.add("sc", cfg::Setting::TypeInt) = physmap_sc;
        root_pm.add("sr", cfg::Setting::TypeInt) = physmap_sr;
        root_pm.add("scale", cfg::Setting::TypeFloat) = physmap_scale;
        root_pm.add("scalespace", cfg::Setting::TypeFloat) = physmap_scalespace;
    }
    else{
        root.add("efield",cfg::Setting::TypeGroup);
        cfg::Setting & root_efield = root["efield"];
        root_efield.add("path", cfg::Setting::TypeString) = EFfile;
        root_efield.add("out", cfg::Setting::TypeString) = EFfile_out;
        root_efield.add("sc", cfg::Setting::TypeInt) = efield_sc;
        root_efield.add("sr", cfg::Setting::TypeInt) = efield_sr;
        root_efield.add("scale", cfg::Setting::TypeFloat) = efield_scale;
        root_efield.add("scalespace", cfg::Setting::TypeFloat) = efield_scalespace;
        
        root.add("wfield",cfg::Setting::TypeGroup);
        cfg::Setting & root_wfield = root["wfield"];
        root_wfield.add("path", cfg::Setting::TypeString) = WFfile;
        root_wfield.add("out", cfg::Setting::TypeString) = WFfile_out;
        root_wfield.add("sc", cfg::Setting::TypeInt) = wfield_sc;
        root_wfield.add("sr", cfg::Setting::TypeInt) = wfield_sr;
        root_wfield.add("scale", cfg::Setting::TypeFloat) = wfield_scale;
        root_wfield.add("scalespace", cfg::Setting::TypeFloat) = wfield_scalespace;
        
        root.add("emob",cfg::Setting::TypeGroup);
        cfg::Setting & root_emob = root["emob"];
        root_emob.add("path", cfg::Setting::TypeString) = EMOBfile;
        root_emob.add("out", cfg::Setting::TypeString) = EMOBfile_out;
        root_emob.add("sc", cfg::Setting::TypeInt) = emob_sc;
        root_emob.add("sr", cfg::Setting::TypeInt) = emob_sr;
        root_emob.add("scale", cfg::Setting::TypeFloat) = emob_scale;
        root_emob.add("scalespace", cfg::Setting::TypeFloat) = emob_scalespace;
        
        root.add("hmob",cfg::Setting::TypeGroup);
        cfg::Setting & root_hmob = root["hmob"];
        root_hmob.add("path", cfg::Setting::TypeString) = HMOBfile;
        root_hmob.add("out", cfg::Setting::TypeString) = HMOBfile_out;
        root_hmob.add("sc", cfg::Setting::TypeInt) = hmob_sc;
        root_hmob.add("sr", cfg::Setting::TypeInt) = hmob_sr;
        root_hmob.add("scale", cfg::Setting::TypeFloat) = hmob_scale;
        root_hmob.add("scalespace", cfg::Setting::TypeFloat) = hmob_scalespace;
    }
    
    maps::preparemaps(root);

    return 0;
}





#endif /* MAIN_MAPS_INL_ */
