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
 * HydraSim.h
 *
 *  Created on: 12/11/2018
 *      Author: Andrea Contu
 */
#include <simulation/Evolve.h>
#include <simulation/simulation.h>

#include <hydra/functions/Gaussian.h>

void signal_simulation(cfg::Setting const& root, int icfg=-1){
    
    //configuration list
    std::vector<singleconf> configlist=FillConfig(root, icfg);
    
    INFO_LINE("Requested "<<configlist.size() << " simulation(s)")
    
    //take time
    auto start_all = std::chrono::high_resolution_clock::now();

    //check if output directory exists, if it does not it is created
    std::string outputdir=root["OutputDirectory"];
    mkdir(outputdir.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IROTH);
    
    //getting map configuration
    maps::mapconfig _maps(root);
    
    #ifdef _ROOT_AVAILABLE_
    //3d hist of efield if necessary
    TH3D *h3 = new TH3D();
    for(auto cf:configlist){
        if(std::get<bool>(cf.get("plot"))){
            _maps.getHist3D(h3,0);
            break;
        }
    }
    #endif //_ROOT_AVAILABLE_
    
    
    INFO_LINE("Starting signal simulation...")
    std::cout << MIDDLEROW << std::endl;
    double simulationtime=0;
    
    std::string tmpfile=std::get<std::string>(configlist[0].get("path"));
    
    RunningStateHost_t data_dinit;
    if(tmpfile==DEFAULT_PATH) data_dinit = loaddata::getDummy(  std::get<size_t>(configlist[0].get("nparticles")), 
                                                                std::get<double>(configlist[0].get("length")),std::get<size_t>(configlist[0].get("group")));
    else data_dinit=loaddata::loadfile(tmpfile);
    
    size_t sim=0;
    for(auto cf:configlist){
    {
        std::cout << MIDDLEROW << std::endl;
        //get configuration
        auto bunchsize      =   std::get<size_t>(cf.get("bunchsize"));
        auto nparticles     =   std::get<size_t>(cf.get("nparticles"));
        auto nsteps         =   std::get<size_t>(cf.get("nsteps"));
        auto group          =   std::get<size_t>(cf.get("group"));
        auto path           =   std::get<std::string>(cf.get("path"));
        auto nickname       =   std::get<std::string>(cf.get("nickname"));
        auto plot           =   std::get<bool>(cf.get("plot"));
        auto extrainfo      =   std::get<bool>(cf.get("extrainfo"));
        auto temperature    =   std::get<double>(cf.get("temperature"));
        auto timestep       =   std::get<double>(cf.get("timestep"));
        auto xshift         =   std::get<double>(cf.get("xshift"));
        auto yshift         =   std::get<double>(cf.get("yshift"));
        auto zshift         =   std::get<double>(cf.get("zshift"));
        auto sig            =   std::get<double>(cf.get("sig"));
        auto amp            =   std::get<double>(cf.get("amp"));
        auto length         =   std::get<double>(cf.get("length"));
        
        auto start = std::chrono::high_resolution_clock::now();
        if(path!=tmpfile){
            tmpfile=path;
            if(tmpfile==DEFAULT_PATH) data_dinit=loaddata::getDummy(nparticles,length,group);
            else data_dinit=loaddata::loadfile(tmpfile,group);
        }
        
        //create running state

        if(nparticles==0 || nparticles>data_dinit.size()) nparticles=data_dinit.size();

        //create running state
        RunningState_t data_d(nparticles,RunningTuple_init); // running state
        data_d.reserve(nparticles);
        //create reduced vector for plotting
        ReducedDataDev_t currents(nsteps);
        currents.reserve(nsteps);
        CurrentVec_t tp_currs;
        
        
        //fill running state with initial particles
        hydra::copy(hydra::make_range(data_dinit.begin(),data_dinit.begin()+nparticles),data_d); // with ranges

        hydra::for_each(data_d, loaddata::Translate(xshift,yshift,zshift,std::get<0>(_maps.Boundaries()),std::get<1>(_maps.Boundaries()),
                                                                    std::get<2>(_maps.Boundaries()),std::get<3>(_maps.Boundaries()),
                                                                    std::get<4>(_maps.Boundaries()),std::get<5>(_maps.Boundaries())));
        
        //creating large host vector to copy bunches of dev vectors
        size_t stp=1;
        if(extrainfo || plot) stp=nsteps;
        UniverseHost_t states_host(stp, StateHost_t(nparticles,StateTuple_init));
        states_host.reserve(stp);
        for(auto& st:states_host){
            st.reserve(nparticles);
        }
        
        size_t nbunches = nsteps/bunchsize;
        size_t rest = nsteps%bunchsize;
        
        INFO_LINE("Simulation "<<sim<<" ("<<nickname<<")")
        if(group>1) INFO_LINE("Particles are grouped in "<< nparticles <<" bunches of "<<group);
        INFO_LINE("Nparticles(dep. charge): "<<nparticles*group<<"("<< (double)(nparticles*group)*HCHARGE/2. <<" C)"<< ", Temperature "<< temperature << ", timestep "<<timestep)
        INFO_LINE("Nsteps: "<<nsteps)
//         if(ampexp[sim]!=0.0) INFO_LINE("Amp_exp "<< cf.get("amp").getd() << ", Sigma_exp "<<cf.get("sig").getd())
        if(xshift!=0.0) INFO_LINE("shifted in x by "<<xshift<< " micron")
        if(yshift!=0.0) INFO_LINE("shifted in y by "<<yshift<< " micron")
        if(zshift!=0.0) INFO_LINE("shifted in z by "<<zshift<< " micron")
        if((extrainfo || plot)) INFO_LINE("Requested "<<nsteps<< " steps. Dividing into "<< nbunches << " bunches of "<<bunchsize<<" steps plus one of "<< rest << " steps")
        
        
        auto A = hydra::Parameter::Create().Name("A").Value(0.0);
        auto B = hydra::Parameter::Create().Name("B").Value(PI);
        auto C = hydra::Parameter::Create().Name("C").Value(2*PI);
        
        auto uniform_theta   = hydra::UniformShape<double>(A,B);
        auto uniform_phi   = hydra::UniformShape<double>(A,C);
        
        auto mean  = hydra::Parameter::Create("mean" ).Value(0.0);
        auto sigma = hydra::Parameter::Create("sigma").Value(1.0);
        
        auto gauss     = hydra::Gaussian<double>(mean, sigma);
        
//         hydra::Random<> Generator( std::chrono::system_clock::now().time_since_epoch().count());
        
        
        size_t niter=0;
        //creating running vectors of states and reserving memory for them and their elements
        size_t bnc=1;
        if((extrainfo || plot)) bnc=bunchsize; 
        UniverseDev_t states(bnc, StateDev_t(nparticles,StateTuple_init));
        states.reserve(bnc);
        for(auto& final_state:states){
            final_state.reserve(nparticles);
        }
        
        size_t cb=0; //count bunches
        
        {
        
            //for loop on all states in the universe
            if(_maps.ismulti()){
                for(size_t id=0; id<=nsteps; id++){
                    double temperature_nmob = (((KB*temperature))/HCHARGE)*(1 + amp*exp(-0.5*pow((timestep*niter - sig)/sig,2)));
                    
                    hydra::fill_random(data_d.begin(_tc_gauss_x), data_d.end(_tc_gauss_x) , gauss,std::rand());
                    hydra::fill_random(data_d.begin(_tc_angle_1), data_d.end(_tc_angle_1) , uniform_phi,std::rand());
                    hydra::fill_random(data_d.begin(_tc_angle_2), data_d.end(_tc_angle_2) , uniform_theta,std::rand());
                    
                    
                    hydra::for_each(data_d, evolve::make_RamoCurrent(niter,
                                                                    timestep,
                                                                    temperature_nmob,
                                                                    _maps.VecXEField(), _maps.VecYEField(), _maps.VecZEField(), _maps.EField(),
                                                                    _maps.VecXEMob(), _maps.VecYEMob(), _maps.VecZEMob(), _maps.EMob(),
                                                                    _maps.VecXHMob(), _maps.VecYHMob(), _maps.VecZHMob(), _maps.HMob(),
                                                                    _maps.VecXWField(), _maps.VecYWField(), _maps.VecZWField(), _maps.WField()));
                    
                    //reducing data
                    
//                     auto int_curr = hydra::reduce ( hydra::make_range(data_d.begin(_tc_charge,_tc_isin,_tc_curr,_tc_issec),data_d.end(_tc_charge,_tc_isin,_tc_curr,_tc_issec)) | (analysis::_SelectChargeAndSec), ReducedTuple_init,analysis::_SumTuples);
//                     auto int_curr = hydra::reduce ( hydra::transform (hydra::make_range(data_d.begin(_tc_charge,_tc_isin,_tc_curr,_tc_issec),data_d.end(_tc_charge,_tc_isin,_tc_curr,_tc_issec)), analysis::SelectChargeAndSec()), ReducedTuple_init,analysis::SumTuples());
                    auto int_curr=hydra_thrust::transform_reduce(data_d.begin(_tc_charge,_tc_isin,_tc_curr,_tc_issec),data_d.end(_tc_charge,_tc_isin,_tc_curr,_tc_issec),analysis::SelectChargeAndSec(), ReducedTuple_init,analysis::SumTuples());
                    tp_currs.push_back(int_curr);
                    
//                     auto int_curr=hydra::reduce(hydra::make_range(data_d.begin(_tc_charge,_tc_isin,_tc_curr,_tc_issec),data_d.end(_tc_charge,_tc_isin,_tc_curr,_tc_issec)),ReducedTuple_init, analysis::SumCarriers());
//                     
                    if(extrainfo || plot){
                        hydra::copy(hydra::make_range(data_d.begin(_tc_charge,_tc_x,_tc_y,_tc_z,_tc_isin,_tc_curr,_tc_issec), data_d.end(_tc_charge,_tc_x,_tc_y,_tc_z,_tc_isin,_tc_curr,_tc_issec)), states[id%bunchsize]);
                        if(id!=0 && (id+1)%bunchsize==0){
                            hydra::copy(states,hydra::make_range(states_host.begin()+cb*bunchsize,states_host.end()));
                            cb++;
                        }
                    }
                    
                    printProgress((double)(niter+1)/(double)nsteps,"Simulation progress:");
                    
                    
                    if(niter==nsteps-1) break;
                    
                    hydra::for_each(data_d, evolve::make_Evolution(niter,
                                                                    timestep,
                                                                    _maps.VecXEField(), _maps.VecYEField(), _maps.VecZEField(), _maps.EField(),
                                                                    _maps.VecXEMob(), _maps.VecYEMob(), _maps.VecZEMob(), _maps.EMob(),
                                                                    _maps.VecXHMob(), _maps.VecYHMob(), _maps.VecZHMob(), _maps.HMob(),
                                                                    std::get<0>(_maps.Boundaries()),std::get<1>(_maps.Boundaries()),
                                                                    std::get<2>(_maps.Boundaries()),std::get<3>(_maps.Boundaries()),
                                                                    std::get<4>(_maps.Boundaries()),std::get<5>(_maps.Boundaries())));
                    
                    niter++;
                    
                    
                }
            }
            else{
                for(size_t id=0; id<=nsteps; id++){
                    double temperature_nmob = (((KB*temperature))/HCHARGE)*(1 + amp*exp(-0.5*pow((timestep*niter - sig)/sig,2)));
                    
                    hydra::fill_random(data_d.begin(_tc_gauss_x), data_d.end(_tc_gauss_x) , gauss,std::rand());
                    hydra::fill_random(data_d.begin(_tc_angle_1), data_d.end(_tc_angle_1) , uniform_phi,std::rand());
                    hydra::fill_random(data_d.begin(_tc_angle_2), data_d.end(_tc_angle_2) , uniform_theta,std::rand());
                    
                    
                    hydra::for_each(data_d, evolve::make_RamoCurrent(niter,
                                                                    timestep,
                                                                    temperature_nmob,
                                                                    _maps.VecX(), _maps.VecY(), _maps.VecZ(), _maps.PhysMap()));
                    
                    //reducing data
//                     auto int_curr = hydra::reduce ( hydra::make_range(data_d.begin(_tc_charge,_tc_isin,_tc_curr,_tc_issec),data_d.end(_tc_charge,_tc_isin,_tc_curr,_tc_issec)) | (analysis::_SelectChargeAndSec), ReducedTuple_init,analysis::_SumTuples);
// //                     auto int_curr = hydra::reduce ( hydra::transform(hydra::make_range(data_d.begin(_tc_charge,_tc_isin,_tc_curr,_tc_issec),data_d.end(_tc_charge,_tc_isin,_tc_curr,_tc_issec)),analysis::SelectChargeAndSec()), ReducedTuple_init,analysis::SumTuples());
// //                     auto r = reduce ( transform (data, functorT), init, functorR)
                    auto int_curr=hydra_thrust::transform_reduce(data_d.begin(_tc_charge,_tc_isin,_tc_curr,_tc_issec),data_d.end(_tc_charge,_tc_isin,_tc_curr,_tc_issec),analysis::SelectChargeAndSec(), ReducedTuple_init,analysis::SumTuples());
                    
//                     auto int_curr=hydra::reduce(hydra::make_range(data_d.begin(_tc_charge,_tc_isin,_tc_curr,_tc_issec),data_d.end(_tc_charge,_tc_isin,_tc_curr,_tc_issec)),ReducedTuple_init,analysis::SumCarriers());
                    
                    tp_currs.push_back(int_curr);
                    
                    if(extrainfo || plot){
                        hydra::copy(hydra::make_range(data_d.begin(_tc_charge,_tc_x,_tc_y,_tc_z,_tc_isin,_tc_curr,_tc_issec), data_d.end(_tc_charge,_tc_x,_tc_y,_tc_z,_tc_isin,_tc_curr,_tc_issec)), states[id%bunchsize]);
                        if(id!=0 && (id+1)%bunchsize==0){
                            hydra::copy(states,hydra::make_range(states_host.begin()+cb*bunchsize,states_host.end()));
                            cb++;
                        }
                    }
                    
                    printProgress((double)(niter+1)/(double)nsteps,"Simulation progress:");
                    
                    
                    if(niter==nsteps-1) break;

                    hydra::for_each(data_d, evolve::make_Evolution(niter,
                                                                    timestep,
                                                                    _maps.VecX(), _maps.VecY(), _maps.VecZ(), _maps.PhysMap()));
                    
                    niter++;
                    
                    
                }
            }
            if(extrainfo || plot) {
                hydra::copy(hydra::make_range(states.begin(),states.begin()+rest),hydra::make_range(states_host.begin()+cb*bunchsize,states_host.end()));
            }
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;
        std::cout << std::endl;
        INFO_LINE("Simulation time : "<<elapsed.count()<<" ms\n")
        simulationtime+=elapsed.count();
        
        //Analysis and plot results
        setstyle(); //set ROOT style
        std::map<std::string,std::string> settings;
        settings["name"]=nickname;
        settings["outputdir"]=outputdir;
        settings["pairs"]=Form("%i",(int)nparticles/2);
        settings["timestep"]=Form("%g s",timestep);
        settings["nstep"]=Form("%i s",(int)nsteps);
        settings["T"]=Form("%0.0f K",temperature);
        if(amp) settings["expl. ampl."]=Form("%g",amp);
        if(sig) settings["expl. #sigma"]=Form("%g",sig);
        settings["x shift"]=Form("%g #mu m",xshift);
        settings["y shift"]=Form("%g #mu m",yshift);
        settings["z shift"]=Form("%g #mu m",zshift);
        if(!(extrainfo || plot)) analysis::AnalyseSim(tp_currs,timestep,settings);
        else{
            analysis::AnalyseSim(tp_currs,timestep,settings);
            analysis::ExtraPlots(states_host, tp_currs, h3, timestep, settings, extrainfo, plot);
        }
        
    }
        sim++;
    }
    
    std::cout << TOPROW << std::endl;

    INFO_LINE("Total simulation time (excluding animated gif generation and/or extra plots) : "<<simulationtime<<" ms")
    
    //calculate total time
    auto end_all = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_all = end_all - start_all;
    std::cout << TOPROW << std::endl;
    INFO_LINE("Job time : "<<elapsed_all.count()<<" s")
    std::cout << TOPROW << std::endl;
    
}
