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
 * HydraSim.h
 *
 *  Created on: 12/11/2018
 *      Author: Andrea Contu
 */
#include <Evolve.h>
#include <simulation.h>

void signal_simulation(cfg::Setting const& root, int icfg=-1){
    
    //configuration list
    std::vector<singleconf> configlist=FillConfig(root, icfg);
    
    INFO_LINE("Requested "<<configlist.size() << " simulation(s)")
    
    //take time
    auto start_all = std::chrono::high_resolution_clock::now();

    //check if output directory exists, if it does not it is created
    std::string outputdir=root["OutputDirectory"];
    mkdir(outputdir.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IROTH);
    
    auto startmapload = std::chrono::high_resolution_clock::now();
    
//     size_t nl=countlines(physmap_path);
    
    //3d hist of efield
    TH3D* h3=NULL;
    
    VecDev_t<double> _xvec, _yvec, _zvec;
    //efield vecs
    VecDev_t<double> _xvec_ef, _yvec_ef, _zvec_ef;
    VecDev_t<double> efieldvecx, efieldvecy, efieldvecz;
    //wfield vecs
    VecDev_t<double> _xvec_wf, _yvec_wf, _zvec_wf;
    VecDev_t<double> wfieldvecx, wfieldvecy, wfieldvecz;
    
    //emob vecs
    VecDev_t<double> _xvec_emob, _yvec_emob, _zvec_emob;
    VecDev_t<double> emobvec;
    
    //hmob vecs
    VecDev_t<double> _xvec_hmob, _yvec_hmob, _zvec_hmob;
    VecDev_t<double> hmobvec;
    
    double _xmin=0, _xmax=0, _ymin=0, _ymax=0, _zmin=0, _zmax=0;
    //load physics maps
    
    
    bool multimap=true;
    //get physicsmap path
    if(root["PhysicsMaps"].exists("map")){  
        
        INFO_LINE("Loading single map")
        multimap=false;
        //get physicsmap path
        std::string physmap_path=root["PhysicsMaps"]["map"];
        maps::physmap<double,11> pm(physmap_path);
        pm.prepare();
        
        pm.fillspacepoints<VecDev_t<double> >(_xvec,_yvec,_zvec);
        pm.fillallquantities<VecDev_t<double> >(efieldvecx,efieldvecy,efieldvecz,emobvec,hmobvec,wfieldvecx,wfieldvecy,wfieldvecz);
        
        _xmin = *(pm.getx().begin());
        _xmax = *(pm.getx().rbegin());
        _ymin = *(pm.gety().begin());
        _ymax = *(pm.gety().rbegin());
        _zmin = *(pm.getz().begin());
        _zmax = *(pm.getz().rbegin());
        
        for(auto cf:configlist){
            if(std::get<bool>(cf.get("plot"))){
                VecHost_t<double> efieldvecx_host, efieldvecy_host, efieldvecz_host;
                pm.fillvector<VecHost_t<double> >(efieldvecx_host,efieldvecy_host,efieldvecz_host);
                h3=analysis::getHist3DDraw( efieldvecx_host, efieldvecy_host, efieldvecz_host, pm.getx(),pm.gety(),pm.getz());
                break;
            }
        }
    }
    else{
        if(root["PhysicsMaps"].exists("efield") && root["PhysicsMaps"].exists("wfield") && root["PhysicsMaps"].exists("emob") && root["PhysicsMaps"].exists("hmob")){
            
            INFO_LINE("Loading separate maps")
            
            std::string efieldmap_path=root["PhysicsMaps"]["efield"];
            
            maps::physmap<double,6> efieldm(efieldmap_path); //E field is a vector
            efieldm.prepare();
            
            efieldm.fillspacepoints<VecDev_t<double> >(_xvec_ef,_yvec_ef,_zvec_ef);
            efieldm.fillvector<VecDev_t<double> >(efieldvecx,efieldvecy,efieldvecz);
            
            for(auto cf:configlist){
                if(std::get<bool>(cf.get("plot"))){
                    VecHost_t<double> efieldvecx_host, efieldvecy_host, efieldvecz_host;
                    efieldm.fillvector<VecHost_t<double> >(efieldvecx_host,efieldvecy_host,efieldvecz_host);
                    h3=analysis::getHist3DDraw( efieldvecx_host, efieldvecy_host, efieldvecz_host, efieldm.getx(),efieldm.gety(),efieldm.getz());
                    break;
                }
            }

            
            std::string wfieldmap_path=root["PhysicsMaps"]["wfield"];
            maps::physmap<double,6> wfieldm(wfieldmap_path); // Weighting field is a vector
            wfieldm.prepare();
            wfieldm.fillspacepoints<VecDev_t<double> >(_xvec_wf,_yvec_wf,_zvec_wf);
            wfieldm.fillvector<VecDev_t<double> >(wfieldvecx,wfieldvecy,wfieldvecz);
            
            std::string emobmap_path=root["PhysicsMaps"]["emob"];
            maps::physmap<double,4> emobm(emobmap_path); // Electron mobility is a scalar
            emobm.prepare();
            emobm.fillspacepoints<VecDev_t<double> >(_xvec_emob,_yvec_emob,_zvec_emob);
            emobm.fillscalar<VecDev_t<double> >(emobvec);
            
            std::string hmobmap_path=root["PhysicsMaps"]["hmob"];
            maps::physmap<double,4> hmobm(hmobmap_path); // Hole mobility is a scalar
            hmobm.prepare();
            hmobm.fillspacepoints<VecDev_t<double> >(_xvec_hmob,_yvec_hmob,_zvec_hmob);
            hmobm.fillscalar<VecDev_t<double> >(hmobvec);
            
            
            _xmin = std::max(*(efieldm.getx().begin()),std::max(*(wfieldm.getx().begin()),std::max(*(emobm.getx().begin()),*(hmobm.getx().begin()))));
            _xmax = std::min(*(efieldm.getx().rbegin()),std::min(*(wfieldm.getx().rbegin()),std::min(*(emobm.getx().rbegin()),*(hmobm.getx().rbegin()))));
            _ymin = std::max(*(efieldm.gety().begin()),std::max(*(wfieldm.gety().begin()),std::max(*(emobm.gety().begin()),*(hmobm.gety().begin()))));
            _ymax = std::min(*(efieldm.gety().rbegin()),std::min(*(wfieldm.gety().rbegin()),std::min(*(emobm.gety().rbegin()),*(hmobm.gety().rbegin()))));
            _zmin = std::max(*(efieldm.getz().begin()),std::max(*(wfieldm.getz().begin()),std::max(*(emobm.getz().begin()),*(hmobm.getz().begin()))));
            _zmax = std::min(*(efieldm.getz().rbegin()),std::min(*(wfieldm.getz().rbegin()),std::min(*(emobm.getz().rbegin()),*(hmobm.getz().rbegin()))));
            
            std::cout << MIDDLEROW << std::endl;
            INFO_LINE("Volume boundaries: [ "<< _xmin << " < x < "<< _xmax << " ] micron, "<<"[ "<< _ymin << " < y < "<< _ymax <<" ] micron, "<<"[ "<< _zmin << " < z < "<< _zmax <<" ] micron")
        }
    }
    
    auto endmapload = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedmapload = endmapload - startmapload;
    INFO_LINE("Map(s) loaded in "<<elapsedmapload.count()<<" ms\n")
    

    INFO_LINE("Starting signal simulation...")
    std::cout << MIDDLEROW << std::endl;
    double simulationtime=0;
    
    std::string tmpfile=std::get<std::string>(configlist[0].get("path"));
    
    RunningStateHost_t data_dinit;
    if(tmpfile==DEFAULT_PATH) data_dinit = loaddata::getDummy(  std::get<size_t>(configlist[0].get("nparticles")), 
                                                                std::get<double>(configlist[0].get("length")));
    else data_dinit=loaddata::loadfile(tmpfile);
    
    
    
    size_t sim=0;
    for(auto cf:configlist){
    {
        std::cout << MIDDLEROW << std::endl;
        //get configuration
        auto bunchsize      =   std::get<size_t>(cf.get("bunchsize"));
        auto nparticles     =   std::get<size_t>(cf.get("nparticles"));
        auto nsteps         =   std::get<size_t>(cf.get("nsteps"));
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
        
//         cf.Print();
        
        auto start = std::chrono::high_resolution_clock::now();
        if(path!=tmpfile){
            tmpfile=path;
            if(tmpfile==DEFAULT_PATH) data_dinit=loaddata::getDummy(nparticles,length);
            else data_dinit=loaddata::loadfile(tmpfile);
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

        hydra::for_each(data_d, loaddata::Translate(xshift,yshift,zshift,_xmin,_xmax,_ymin,_ymax,_zmin,_zmax));
        
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
        INFO_LINE("Nparticles(dep. charge): "<<nparticles<<"("<< (double)(nparticles)*HCHARGE/2. <<" C)"<< ", Temperature "<< temperature << ", timestep "<<timestep)
        INFO_LINE("Nsteps: "<<nsteps)
//         if(ampexp[sim]!=0.0) INFO_LINE("Amp_exp "<< cf.get("amp").getd() << ", Sigma_exp "<<cf.get("sig").getd())
        if(xshift!=0.0) INFO_LINE("shifted in x by "<<xshift<< " micron")
        if(yshift!=0.0) INFO_LINE("shifted in y by "<<yshift<< " micron")
        if(zshift!=0.0) INFO_LINE("shifted in z by "<<zshift<< " micron")
        if((extrainfo || plot)) INFO_LINE("Requested "<<nsteps<< " steps. Dividing into "<< nbunches << " bunches of "<<bunchsize<<" steps plus one of "<< rest << " steps")
        
    
        hydra::Random<> Generator( std::chrono::system_clock::now().time_since_epoch().count());
        
        
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
            if(multimap){
                for(size_t id=0; id<=nsteps; id++){
                    double temperature_nmob = (((KB*temperature))/HCHARGE)*(1 + amp*exp(-0.5*pow((timestep*niter - sig)/sig,2)));
                    
                    Generator.SetSeed(std::rand());
                    Generator.Gauss(0.0, 1.0, data_d.begin(_tc_gauss_x), data_d.end(_tc_gauss_x));
                    
                    Generator.SetSeed(std::rand());
                    Generator.Uniform(0.0, 2*PI, data_d.begin(_tc_angle_1), data_d.end(_tc_angle_1));
                    Generator.SetSeed(std::rand());
                    Generator.Uniform(0.0, PI, data_d.begin(_tc_angle_2), data_d.end(_tc_angle_2));
                    
                    
                    hydra::for_each(data_d, evolve::ApplyRamo_multi(niter,
                                                                    timestep,
                                                                    temperature_nmob,
                                                                    _xvec_ef, _yvec_ef, _zvec_ef, efieldvecx, efieldvecy, efieldvecz,
                                                                    _xvec_emob, _yvec_emob, _zvec_emob, emobvec,
                                                                    _xvec_hmob, _yvec_hmob, _zvec_hmob, hmobvec,
                                                                    _xvec_wf, _yvec_wf, _zvec_wf, wfieldvecx, wfieldvecy, wfieldvecz));
                    
                    //reducing data
                    auto int_curr=HYDRA_EXTERNAL_NS::thrust::transform_reduce(data_d.begin(_tc_charge,_tc_isin,_tc_curr,_tc_issec),data_d.end(_tc_charge,_tc_isin,_tc_curr,_tc_issec),analysis::SelectChargeAndSec(), ReducedTuple_init,analysis::SumTuples());
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

                    hydra::for_each(data_d, evolve::Evolve_multi(niter,
                                                            timestep,
                                                            _xvec_ef, _yvec_ef, _zvec_ef, efieldvecx, efieldvecy, efieldvecz,
                                                            _xvec_emob, _yvec_emob, _zvec_emob, emobvec,
                                                            _xvec_hmob, _yvec_hmob, _zvec_hmob, hmobvec,
                                                            _xmin,_xmax,_ymin,_ymax,_zmin,_zmax));
                    niter++;
                    
                    
                }
            }
            else{
                for(size_t id=0; id<=nsteps; id++){
                    double temperature_nmob = (((KB*temperature))/HCHARGE)*(1 + amp*exp(-0.5*pow((timestep*niter - sig)/sig,2)));
                    
                    Generator.SetSeed(std::rand());
                    Generator.Gauss(0.0, 1.0, data_d.begin(_tc_gauss_x), data_d.end(_tc_gauss_x));
                    
                    Generator.SetSeed(std::rand());
                    Generator.Uniform(0.0, 2*PI, data_d.begin(_tc_angle_1), data_d.end(_tc_angle_1));
                    Generator.SetSeed(std::rand());
                    Generator.Uniform(0.0, PI, data_d.begin(_tc_angle_2), data_d.end(_tc_angle_2));
                    
                    
                    hydra::for_each(data_d, evolve::ApplyRamo(niter,
                                                                    timestep,
                                                                    temperature_nmob,
                                                                    _xvec, _yvec, _zvec,
                                                                    efieldvecx, efieldvecy, efieldvecz,
                                                                    emobvec,hmobvec,
                                                                    wfieldvecx, wfieldvecy, wfieldvecz));
                    
                    //reducing data
                    auto int_curr=HYDRA_EXTERNAL_NS::thrust::transform_reduce(data_d.begin(_tc_charge,_tc_isin,_tc_curr,_tc_issec),data_d.end(_tc_charge,_tc_isin,_tc_curr,_tc_issec),analysis::SelectChargeAndSec(), ReducedTuple_init,analysis::SumTuples());
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

                    hydra::for_each(data_d, evolve::Evolve(niter,
                                                            timestep,
                                                            _xvec, _yvec, _zvec,
                                                            efieldvecx, efieldvecy, efieldvecz,
                                                            emobvec,hmobvec));
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
