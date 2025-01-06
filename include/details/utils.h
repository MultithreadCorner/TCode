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
 * utils.h
 *
 *  Created on: 12/11/2018
 *  Author: Andrea Contu
 *  Updated on: 15/01/2024
 *  Author: Angelo Loi
 */

#ifndef UTILS_H_
#define UTILS_H_

#define VERBOSE_LINE(flag, message)\
if(flag){\
	std::cout<< "\033[1;35mVerbose: \033[0m"\
			 << message << std::endl;\
}

#define INFO_LINE(message)\
	std::cout<< "\033[1;34mInfo: \033[0m"\
			 << message << std::endl;\

#define WARNING_LINE(message)\
	std::cout<< "\033[1;35mWarning: \033[0m"\
			 << message << std::endl;\

#define ERROR_LINE(message)\
	std::cout<< "\033[1;31mError: \033[0m"\
			 << message << std::endl;\
			 
//stdlib
#include <iostream>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <chrono>
#include <random>
#include <algorithm>
#include <tuple>
#include <map>
#include <set>
#include <sys/types.h>
#include <sys/stat.h>
#include <type_traits>

//Hydra
#include <hydra/device/System.h>
#include <hydra/host/System.h>
#include <hydra/Function.h>
// #include <hydra/FunctionWrapper.h>
#include <hydra/Random.h>
#include <hydra/Algorithm.h>
#include <hydra/Tuple.h>
#include <hydra/Distance.h>
#include <hydra/multiarray.h>
#include <hydra/multivector.h>
#include <hydra/Parameter.h>
#include <hydra/UserParameters.h>
#include <hydra/Filter.h>
#include <hydra/DenseHistogram.h>
#include <hydra/Placeholders.h>
#include <hydra/Range.h>
#include <hydra/Zip.h>
// #include <hydra/detail/external/Eigen/Dense>


//ROOT
#ifdef _ROOT_AVAILABLE_
#include "TCanvas.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TApplication.h"
#include "TGraphTime.h"
#include "TMarker.h"
#include "TFile.h"
#include "TPad.h"
#include "TLegend.h"
#include "TTree.h"
#endif //_ROOT_AVAILABLE_

//printouts
#define TOPROW    "===================================================================================================="
#define MIDDLEROW "----------------------------------------------------------------------------------------------------"

//command line
#include <tclap/CmdLine.h>

//configuration file
#include <libconfig.h++>

//typedef for setting
typedef std::tuple<int,size_t,double,bool,std::string> _setting;

//config namespace
namespace cfg = libconfig;


/**
 * Checks if a given setting (path) is present in the configuration file.
 * Exit with "FAILURE" status if the setting is not found.
 *
 * @param Cfg libconfig::Config object representing the
 * @param path Setting path
 */
void checker( const cfg::Config &Cfg, const char* path ){
    
    if( !( Cfg.exists(path)   ))
    {
        std::cerr << path << " is not set. Exiting..." << std::endl ;
        exit(EXIT_FAILURE);
    }
    else
        std::cout<< path <<" is set." << std::endl ;
    
}

inline void filldouble(std::vector<double> *v, cfg::Setting &set){
    int ln=set.getLength();
    double temp=0;
    if(ln==0){
        temp=set;
        v->push_back(temp);
    }
    else{
        for(int l=0;l<ln;l++){
            temp=set[l];
            v->push_back(temp);
        }
    }
}

inline void fillint(std::vector<size_t> *v, cfg::Setting &set){
    int ln=set.getLength();
    int temp=0;
    if(ln==0){
        temp=set;
        if(temp == -1 )v->push_back(temp);
        else v->push_back(abs(temp));
    }
    else{
        for(int l=0;l<ln;l++){
            temp=set[l];
            //             std::cout << temp << std::endl;
            if(temp == -1 )v->push_back(temp);
            else v->push_back(abs(temp));
        }
    }
}

std::string return_current_time_and_date(){
    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
    return ss.str();
}


class singleconf{
    
private:
    std::map<std::string, _setting> m;
    std::map<std::string, size_t> pos;
public:
    singleconf(){
        
    }
    
    void set(std::string nm, size_t a){
        m.insert(std::make_pair(nm,_setting(0,a,0,0,"")));
        pos.insert(std::make_pair(nm,1));
        
    }
    void set(std::string nm, int a){
        m.insert(std::make_pair(nm,_setting(a,0,0,0,"")));
        pos.insert(std::make_pair(nm,0));
    }
    void set(std::string nm, double a){
        m.insert(std::make_pair(nm,_setting(0,0,a,0,"")));
        pos.insert(std::make_pair(nm,2));
    }
    void set(std::string nm, std::string a){
        m.insert(std::make_pair(nm,_setting(0,0,0,0,a)));
        pos.insert(std::make_pair(nm,4));
    }
    void set(std::string nm, bool a){
        m.insert(std::make_pair(nm,_setting(0,0,0,a,"")));
        pos.insert(std::make_pair(nm,3));
    }
    _setting get(std::string nm){
        return m.find(nm)->second;
    }
    void Print(){
        std::cout << "Configuration settings"<< std::endl;
        for(auto _m:m){
            if(pos.find(_m.first)->second==0) std::cout << _m.first << "\t" << std::get<0>(_m.second) << std::endl;
            if(pos.find(_m.first)->second==1) std::cout << _m.first << "\t" << std::get<1>(_m.second) << std::endl;
            if(pos.find(_m.first)->second==2) std::cout << _m.first << "\t" << std::get<2>(_m.second) << std::endl;
            if(pos.find(_m.first)->second==3) std::cout << _m.first << "\t" << std::get<3>(_m.second) << std::endl;
            if(pos.find(_m.first)->second==4) std::cout << _m.first << "\t" << std::get<4>(_m.second) << std::endl;
        }
    }
};

//function to create configuration for one simulation
std::vector<singleconf> FillOneConfig(cfg::Setting const& root, int icfg){
        std::vector<singleconf> sets;
        
        std::string file_single_config, file_single_name_config, custom;
        bool draw_config=false, extrainfo_config=false;
        int bunchsize_config;
        std::vector<size_t> nsteps_config, nparticles_config, group_config;
        std::vector<double> temperature_config, timestep_config, ampexp_config, sigmaexp_config, xshift_config, yshift_config, zshift_config, length_config;
        
	if(root["InputData"][icfg].exists("customdeposit")) root["InputData"][icfg].lookupValue("customdeposit",custom);	//<<----nuovo
        else custom=DEFAULT_DEPOSIT_BUILD;
	
	if(root["InputData"][icfg].exists("path")) root["InputData"][icfg].lookupValue("path",file_single_config);
        else file_single_config=DEFAULT_PATH;
        
        if(root["InputData"][icfg].exists("name")) root["InputData"][icfg].lookupValue("name",file_single_name_config);
        else{
            ERROR_LINE("Setting \"name\" must be specified!")
            exit(1);
        }
        
        if(root["InputData"][icfg].exists("plot")) draw_config=root["InputData"][icfg]["plot"];
        else draw_config=((bool)DEFAULT_DRAW);
        
        if(root["InputData"][icfg].exists("extrainfo")) extrainfo_config=root["InputData"][icfg]["extrainfo"];
        else extrainfo_config=((bool)DEFAULT_EXTRAINFO);
        
        if(root["InputData"][icfg].exists("BunchSize")) bunchsize_config=root["InputData"][icfg]["BunchSize"];
        else bunchsize_config=(DEFAULT_BUNCHSIZE);
        
        if(root["InputData"][icfg].exists("steps")) {
            fillint(&nsteps_config,root["InputData"][icfg]["steps"]);
        }
        else{
            ERROR_LINE("Setting \"steps\" must be specified!")
            exit(1);
        }
        
        if(root["InputData"][icfg].exists("particles")) fillint(&nparticles_config,root["InputData"][icfg]["particles"]);
        else nparticles_config.push_back(DEFAULT_PARTICLES);
        
        if(root["InputData"][icfg].exists("group")) fillint(&group_config,root["InputData"][icfg]["group"]);
        else group_config.push_back(DEFAULT_GROUP);
        
        if(root["InputData"][icfg].exists("timestep")) filldouble(&timestep_config,root["InputData"][icfg]["timestep"]);
        else{
            ERROR_LINE("Setting \"timestep\" must be specified!")
            exit(1);
        }
        
        if(root["InputData"][icfg].exists("T")) filldouble(&temperature_config,root["InputData"][icfg]["T"]);
        else{
            ERROR_LINE("Setting \"T\" must be specified!")
            exit(1);
        }
        
        if(root["InputData"][icfg].exists("Aexpl"))filldouble(&ampexp_config,root["InputData"][icfg]["Aexpl"]);
        else ampexp_config.push_back((double)DEFAULT_AMP);
        if(root["InputData"][icfg].exists("Sexpl"))filldouble(&sigmaexp_config,root["InputData"][icfg]["Sexpl"]);
        else sigmaexp_config.push_back((double)DEFAULT_SIGMA);
        
        if(root["InputData"][icfg].exists("Xshift"))filldouble(&xshift_config,root["InputData"][icfg]["Xshift"]);
        else xshift_config.push_back((double)DEFAULT_X);
        
        if(root["InputData"][icfg].exists("Yshift"))filldouble(&yshift_config,root["InputData"][icfg]["Yshift"]);
        else yshift_config.push_back((double)DEFAULT_Y);
        
        if(root["InputData"][icfg].exists("Zshift"))filldouble(&zshift_config,root["InputData"][icfg]["Zshift"]);
        else zshift_config.push_back((double)DEFAULT_Z);
        
        if(root["InputData"][icfg].exists("Length")) filldouble(&length_config,root["InputData"][icfg]["Length"]);
        else length_config.push_back((double)DEFAULT_LENGTH);
        
        size_t nsims= nsteps_config.size()*nparticles_config.size()*timestep_config.size()*temperature_config.size()*ampexp_config.size()*sigmaexp_config.size()*xshift_config.size()*yshift_config.size()*zshift_config.size();
        INFO_LINE("Input data "<<file_single_name_config<<" configuration loaded, asking for "<<nsims<< " simulation(s)")
        if(nsims==0){
            ERROR_LINE("Please check settings of Input data "<<file_single_name_config<<"!")
            exit(1);
        }
        int count=0;
        for(auto s:nsteps_config){
            for(auto p:nparticles_config){
                for(auto t:temperature_config){
                    for(auto tm:timestep_config){
                        for(auto a:ampexp_config){
                            for(auto si:sigmaexp_config){
                                for(auto xs:xshift_config){
                                    for(auto ys:yshift_config){
                                        for(auto zs:zshift_config){
                                            for(auto ln:length_config){
                                                for(auto gr:group_config){
                                                    singleconf tmpconf;
                                                    tmpconf.set("path",file_single_config);
						    tmpconf.set("customdeposit",custom);
                                                    tmpconf.set("nickname",file_single_name_config+"_"+std::to_string(count));
                                                    tmpconf.set("plot",draw_config);
                                                    tmpconf.set("extrainfo",extrainfo_config);
                                                    tmpconf.set("bunchsize",(size_t)bunchsize_config);
                                                    if(p<=0) tmpconf.set("nparticles",(size_t)0);
                                                    else tmpconf.set("nparticles",p);
                                                    tmpconf.set("nsteps",s);
                                                    tmpconf.set("temperature",t);
                                                    tmpconf.set("timestep",tm);
                                                    tmpconf.set("amp",a);
                                                    tmpconf.set("sig",si);
                                                    tmpconf.set("xshift",xs);
                                                    tmpconf.set("yshift",ys);
                                                    tmpconf.set("zshift",zs);
                                                    tmpconf.set("length",ln);
                                                    tmpconf.set("group",gr);
                                                    sets.push_back(tmpconf);
                                                    count++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        return sets;
}


std::vector<singleconf> FillConfig(cfg::Setting const& root, int icfg=-1){
    
    std::vector<singleconf> sets;
    size_t length = root["InputData"].getLength();
    INFO_LINE("Found "<<length<<" input data")
        
    if(icfg>=0 && abs(icfg)< length){
        INFO_LINE("Processing only InputData "<<icfg)
        auto tmpvec=FillOneConfig(root,icfg);
        sets.insert(sets.end(),tmpvec.begin(),tmpvec.end());
    }
    else if(icfg>0 && abs(icfg)>=length){
        ERROR_LINE("You selected an InputData slot that does not exists!")
        exit(1);
    }
    else{
        for(unsigned int i=0; i<length; i++){
            if(root["InputData"][i].exists("process")){
//                 //skip if process flag is set to false
                if((bool)(root["InputData"][i]["process"])==false) continue;
            }
            auto tmpvec=FillOneConfig(root,i);
            sets.insert(sets.end(),tmpvec.begin(),tmpvec.end());
        }
    }
    return sets;
}

//progress bar
void printProgress (double percentage, std::string text){
    int PBWIDTH=50;
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r %s %3d%% [%.*s%*s]", text.c_str(), val, lpad, "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||", rpad, "");
    fflush (stdout);
}

//header
void printHeader(){
    std::cout << TOPROW << std::endl;
    std::cout << "                 Welcome to "<<PROJECT_NAME<< " "<< VERSION <<", a fast silicon sensor simulation               " << std::endl;
    std::cout << MIDDLEROW << std::endl;
    std::cout << "     Developed by Andrea Contu (andrea.contu@ca.infn.it) and Angelo Loi (angelo.loi@ca.infn.it)     " << std::endl;
//    std::cout << MIDDLEROW << std::endl;
//    std::cout << "      Powered by HYDRA (Header only library for data analysis in massively parallel platforms.)     " << std::endl;
//     std::cout << "          Developed by Antonio Augusto Alves Junior, https://github.com/MultithreadCorner           " << std::endl;
    std::cout << TOPROW << std::endl;
}


size_t countlines(std::string _file){
    int numLines = 0;
    std::ifstream in(_file);
    std::string unused;
    while ( std::getline(in, unused) )
    ++numLines;
    in.close();
    return numLines;
}

// #include <tuple>
// #include <iostream>
// #include <type_traits>

template <size_t n, typename... T>
typename std::enable_if<(n >= sizeof...(T))>::type
    print_tuple(std::ostream&, const hydra::tuple<T...>&)
{}

template <size_t n, typename... T>
typename std::enable_if<(n < sizeof...(T))>::type
    print_tuple(std::ostream& os, const hydra::tuple<T...>& tup)
{
    if (n != 0)
        os << "\t";
    os << hydra::get<n>(tup);
    print_tuple<n+1>(os, tup);
}

template <typename... T>
std::ostream& operator<<(std::ostream& os, const hydra::tuple<T...>& tup)
{
    os << "";
    print_tuple<0>(os, tup);
    return os << "";
}

//get sign
template <typename T> 
__hydra_dual__
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

struct vectorLength
{
    
    template<typename Particle>
    __hydra_dual__
    double operator()(Particle p){
        return sqrt(hydra::get<0>(p)*hydra::get<0>(p) + hydra::get<1>(p)*hydra::get<1>(p) +hydra::get<2>(p)*hydra::get<2>(p));
                                                                                
    }
    
};

#endif /* UTILS_H_ */
