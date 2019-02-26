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
 * preparemaps.h
 *
 *  Created on: 12/11/2018
 *  Author: Andrea Contu
 */

#ifndef __PREPAREMAPS_H__
#define __PREPAREMAPS_H__

#include <hydra/detail/Hash.h>
using namespace hydra::placeholders;
namespace maps{
//     typedef hydra::tuple<double,double,double,double,double,double> tuplemap_t;
    typedef hydra::tuple<double,double,double> spacepoint_t;
    
    struct hashfunc{
        template<typename T>
        auto operator()(T t){
            return hydra::detail::hash_tuple(t);
        }
    };
    
    //class for single physics map
    template<typename MP, size_t N>
    class physmap{
    private:
        //data
        bool _isown=false;
        std::string _filename="none";
        bool _sparse=false;
        size_t _cskip=0;
        size_t _rskip=0;
        bool _gridyfied=false;
        double _scalespace=1.;
        double _scalefield=1.;
        //data points container
        hydra::multiarray<MP,3,hydra::host::sys_t> _extremes;
        hydra::multiarray<MP,N,hydra::host::sys_t> _data;
        
        //own grid
        std::set<MP> _gridx;
        std::set<MP> _gridy;
        std::set<MP> _gridz;
        
        //methods
        void checkfile(); // load data 
        void loaddata(); // load data 
        void removeduplicates(); //data
        void orderdata();// order data
        void calcgrid(); //calculate own grid
        void selfgridify(){this->gridify(_gridx,_gridy,_gridz);}
        std::vector<std::string> genheader();
        

    public:
        
        physmap(){
        }
        
        physmap(std::string s, size_t skipcolumns=0, size_t skiprows=0, double scale=1., double scalespace=1., bool sparse=false){
            this->init(s,skipcolumns,skiprows,scale,scalespace,sparse);
        }
        
        void init(std::string s, size_t skipcolumns=0, size_t skiprows=0, double scale=1., double scalespace=1., bool sparse=false){
            _filename=s;
            _cskip=skipcolumns;
            _rskip=skiprows;
            _sparse=false;
            _gridyfied=false;
            _scalefield=scale;
            _scalespace=scalespace;
        }
        
        void prepare(){
            this->checkfile();
            this->loaddata();
            if(_isown) this->calcgrid();
            else{
                this->orderdata();
                this->removeduplicates();
                this->calcgrid();
            }
        }
        
        void set_scalespace(double A){_scalespace=A;}
        
        void set_scale(double A){_scalefield=A;}
        
        double get_scalespace(double A){return _scalespace;}
        
        double get_scale(double A){return _scalefield;}
        
        void gridify(std::set<MP> xs, std::set<MP> ys, std::set<MP> zs); //get in a grid
        
        void write(std::string outputpath);
        
        hydra::multiarray<MP,N,hydra::host::sys_t> getdata(){return _data;}
        std::set<MP> getx(){return _gridx;}
        std::set<MP> gety(){return _gridy;}
        std::set<MP> getz(){return _gridz;}
        template <typename Container>
        void fillspacepoints(Container &vx, Container &vy, Container &vz){
            vx.resize(_gridx.size());
            vy.resize(_gridy.size());
            vz.resize(_gridz.size());
            hydra::copy(_gridx,vx);
            hydra::copy(_gridy,vy);
            hydra::copy(_gridz,vz);
        }
        
        template <typename Container>
        void fillscalar(Container &v){
            v.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_3),_data.end(_3)),v);
        }
        
        template <typename Container>
        void fillvector(Container &vx, Container &vy, Container &vz){
            vx.resize(_data.size());
            vy.resize(_data.size());
            vz.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_3),_data.end(_3)),vx);
            hydra::copy(hydra::make_range(_data.begin(_4),_data.end(_4)),vy);
            hydra::copy(hydra::make_range(_data.begin(_5),_data.end(_5)),vz);
        }
        
        template <typename Container>
        void fillallquantities(Container &ex, Container &ey, Container &ez, Container &emob, Container &hmob, Container &wx, Container &wy, Container &wz){
            ex.resize(_data.size());
            ey.resize(_data.size());
            ez.resize(_data.size());
            wx.resize(_data.size());
            wy.resize(_data.size());
            wz.resize(_data.size());
            emob.resize(_data.size());
            hmob.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_3),_data.end(_3)),ex);
            hydra::copy(hydra::make_range(_data.begin(_4),_data.end(_4)),ey);
            hydra::copy(hydra::make_range(_data.begin(_5),_data.end(_5)),ez);
            hydra::copy(hydra::make_range(_data.begin(_6),_data.end(_6)),emob);
            hydra::copy(hydra::make_range(_data.begin(_7),_data.end(_7)),hmob);
            hydra::copy(hydra::make_range(_data.begin(_8),_data.end(_8)),wx);
            hydra::copy(hydra::make_range(_data.begin(_9),_data.end(_9)),wy);
            hydra::copy(hydra::make_range(_data.begin(_10),_data.end(_10)),wz);
        }
        
    };
    
    //load data from file
    template<typename MP, size_t N>
    void physmap<MP,N>::checkfile(){
        std::cout << TOPROW << std::endl;
        INFO_LINE("Loading map "<< _filename)
        size_t checks=6;
        size_t ck=0;
        std::ifstream infile;
        infile.open(_filename);// file containing numbers in 3 columns 
        if(infile.fail()){ 
            ERROR_LINE("error loading map"+_filename+"!")
            exit(1); // no point continuing if the file didn't open...
        }
        
        std::vector<std::string> comp;
        comp.push_back((std::string)"<"+(std::string)PROJECT_NAME+(std::string)"_MAP>");
        comp.push_back((std::string)PROJECT_NAME + (std::string)" v");
        comp.push_back((std::string)"Date:");
        comp.push_back((std::string)"Npoints:");
        comp.push_back((std::string)"NPX:");
        comp.push_back((std::string)"</"+(std::string)PROJECT_NAME+(std::string)"_MAP>");
        
        std::string line;
        std::vector<std::string> lines;
        
        for(size_t i=0;i<checks;i++){
            std::getline(infile,line);
            if(line.find(comp[i])!=std::string::npos) ck++;
//             else std::cout << "CN: " << line << "   "<<comp[i]<<std::endl;
        }
        if(checks==ck){
            INFO_LINE("Map was generated with "<<PROJECT_NAME<<". It will be imported directly. ")
            _isown=true;
            _cskip=0;
            _rskip=6;
            _gridyfied=true;
        }

        infile.close();
    }
    
    //load data from file
    template<typename MP, size_t N>
    void physmap<MP,N>::loaddata(){
        //load
        std::ifstream infile;
        infile.open(_filename);// file containing numbers in 3 columns 
        if(infile.fail()){ 
            ERROR_LINE("error loading map"+_filename+"!")
            exit(1); // no point continuing if the file didn't open...
        }
        std::string line;
        for(size_t i=0; i<_rskip; i++) std::getline(infile,line);
        
        
        std::string tmpstring;
        std::array<MP,N> temparray;
        while(std::getline(infile,line)){
            tmpstring.clear();
            std::stringstream lineStream;
            std::replace(line.begin(), line.end(), ',', ' ');
            lineStream << line;
            for(size_t k=0; k<_cskip; k++) lineStream >> tmpstring;
            
            for(size_t i=0; i<N; i++) lineStream >> temparray[i];
            temparray[0]=temparray[0]*_scalespace;
            temparray[1]=temparray[1]*_scalespace;
            temparray[2]=temparray[2]*_scalespace;
            for(size_t i=3; i<N; i++) temparray[i]=temparray[i]*_scalefield;
            auto temp=hydra::detail::arrayToTuple(temparray);
            _data.push_back(temp);
        }
        
        infile.close(); 
        INFO_LINE("Map loaded with "<<_data.size() << " entries")
    }
    
    //remove duplicate entries
    template<typename MP, size_t N>
    void physmap<MP,N>::removeduplicates(){
        hydra::multiarray<MP,N,hydra::host::sys_t> _datandp;
        size_t loadsize=_data.size();
        bool ct=true;
        typename hydra::multiarray<MP,N,hydra::host::sys_t>::value_type tmpt;
        for(auto e:_data){
            if(ct){
                _datandp.push_back(e); 
                ct=false;
            }
            else{
                if(hydra::get<0>(e)!=hydra::get<0>(tmpt) || hydra::get<1>(e)!=hydra::get<1>(tmpt) || hydra::get<2>(e)!=hydra::get<2>(tmpt)) _datandp.push_back(e);
            }
            tmpt=e;
        }
        _data.resize(_datandp.size());
        hydra::copy(_datandp,_data);
        if(loadsize!=_data.size()) INFO_LINE("Removed "<< loadsize - _data.size() << " duplicate points, new number of entries  is " << _data.size())
        
    }
    
    //calculate grid
    template<typename MP, size_t N>
    void physmap<MP,N>::calcgrid(){
        for(auto e:_data){
            _gridx.insert(hydra::get<0>(e));
            _gridy.insert(hydra::get<1>(e));
            _gridz.insert(hydra::get<2>(e));
        }
        
        //set extremes
        _extremes.push_back(spacepoint_t(*_gridx.begin(),*_gridy.begin(),*_gridz.begin()));
        _extremes.push_back(spacepoint_t(*_gridx.rbegin(),*_gridy.begin(),*_gridz.begin()));
        _extremes.push_back(spacepoint_t(*_gridx.begin(),*_gridy.rbegin(),*_gridz.begin()));
        _extremes.push_back(spacepoint_t(*_gridx.begin(),*_gridy.begin(),*_gridz.rbegin()));
        _extremes.push_back(spacepoint_t(*_gridx.rbegin(),*_gridy.rbegin(),*_gridz.begin()));
        _extremes.push_back(spacepoint_t(*_gridx.begin(),*_gridy.rbegin(),*_gridz.rbegin()));
        _extremes.push_back(spacepoint_t(*_gridx.rbegin(),*_gridy.begin(),*_gridz.rbegin()));
        _extremes.push_back(spacepoint_t(*_gridx.rbegin(),*_gridy.rbegin(),*_gridz.rbegin()));
        
        INFO_LINE("Range in x/y/z [microns]: [ "<< *_gridx.begin() << " , "<< *_gridx.rbegin() << " ]/" 
                                         << "[ "<< *_gridy.begin() << " , "<< *_gridy.rbegin() << " ]/"
                                         << "[ "<< *_gridz.begin() << " , "<< *_gridz.rbegin() << " ]")
        
        
        INFO_LINE("Single entries in x/y/z: " << _gridx.size() << "/" << _gridy.size() << "/" << _gridz.size())
        
//         for(auto e:_gridy) std::cout << e << std::endl;
        
        
        if(!_sparse && (_gridx.size()* _gridy.size()* _gridz.size())!=_data.size()){
            if((_gridx.size()* _gridy.size()* _gridz.size())<_data.size()){
                ERROR_LINE("Calculated grid has less points than map entries! It could be corrupted...")
                exit(1);
            }
            WARNING_LINE("Points are not in a grid! Missing "<< size_t(_gridx.size()* _gridy.size()* _gridz.size()) - _data.size() << " points, filling them with zeroes.")
            
            hydra::host::vector<size_t> _datahash_host(_data.size());
            hydra::device::vector<size_t> _datahash(_data.size());
            _datahash.reserve(_data.size());
            hydra::transform(hydra::make_range(_data.begin(_0,_1,_2),_data.end(_0,_1,_2)),_datahash_host, hashfunc());
            hydra::copy(_datahash_host,_datahash);
            
            std::array<MP,N> fd;
            size_t ct=0;
            for(auto px:_gridx){
                for(auto py:_gridy){
                    for(auto pz:_gridz){
                        auto myhash=hydra::detail::hash_tuple(spacepoint_t(px,py,pz));
                        
                        if((HYDRA_EXTERNAL_NS::thrust::find(_datahash.begin(),_datahash.end(), myhash))==_datahash.end()){
                        
                            fd[0]=px;
                            fd[1]=py;
                            fd[2]=pz;
                            for(size_t n=3;n<N;n++) fd[n]=0.;
                            _data.push_back(hydra::detail::arrayToTuple(fd));
//                             std::cout << hydra::detail::arrayToTuple(fd) << std::endl;
//                             std::cout << myhash << std::endl;
                            ct++;
                        }
                    }
                }
            }
            INFO_LINE("added "<<ct<<" points");
            this->removeduplicates();
            this->orderdata();
            _gridyfied = true;
            INFO_LINE("New map size is " << _data.size())
        }
    }
    
    //order the table
    template<typename MP, size_t N>
    void physmap<MP,N>::orderdata(){
        
        struct sortfunc{
            __hydra_dual__
            bool operator()( typename hydra::multiarray<MP,N,hydra::host::sys_t>::value_type left,  typename hydra::multiarray<MP,N,hydra::host::sys_t>::value_type right){
                return (hydra::get<2>(left) < hydra::get<2>(right) || 
                        (hydra::get<2>(left) == hydra::get<2>(right) && hydra::get<0>(left) < hydra::get<0>(right))  ||
                        ((hydra::get<2>(left) == hydra::get<2>(right) && hydra::get<0>(left) <= hydra::get<0>(right)) && (hydra::get<1>(left) < hydra::get<1>(right))));
            }
        };
        
        hydra::sort(_data, sortfunc() );

    }
    
    
    template<typename MP, size_t N>
    void physmap<MP,N>::gridify(std::set<MP> xs, std::set<MP> ys, std::set<MP> z){
        //TO BE IMPLEMENTED
    }
    
    template<typename MP, size_t N>
    std::vector<std::string> physmap<MP,N>::genheader(){
        std::vector<std::string> toout;
        toout.push_back((std::string)"<"+(std::string)PROJECT_NAME+(std::string)"_MAP>");
        toout.push_back((std::string)"Map generated with " + (std::string)PROJECT_NAME + (std::string)" v"+(std::string)VERSION + (std::string)" by importing from "+_filename);
        toout.push_back((std::string)"Date: "+return_current_time_and_date());
        toout.push_back((std::string)"Npoints: "+ std::to_string(_data.size()));
        toout.push_back((std::string)"NPX: "+ std::to_string(_gridx.size()) + (std::string)"  NPY: "+ std::to_string(_gridy.size()) + (std::string)"  NPZ: " +std::to_string(_gridz.size()));
        toout.push_back((std::string)"</"+(std::string)PROJECT_NAME+(std::string)"_MAP>");
        return toout;
    }
    
    template<typename MP, size_t N>
    void physmap<MP,N>::write(std::string outputpath){
        
        if(!_gridyfied){
            ERROR_LINE("Map is not in a grid yet, something is wrong!")
            exit(1);
        }
        INFO_LINE("Saving map as "<<outputpath)
        std::ofstream outfile(outputpath);
        //write header
        auto header=this->genheader();
        for(auto s:header) outfile << s << std::endl;
        //write data
        for(auto d:_data){
            outfile << d <<std::endl;
        }
        outfile.close();
    }
    
    
    
    //class to manage map configuration
    class mapconfig{
    private:
        bool _multimap = false;
        std::tuple<double,double,double,double,double,double> _borders;
        
        maps::physmap<double,11> _pm;
        maps::physmap<double,6> _efieldm;
        maps::physmap<double,6> _wfieldm;
        maps::physmap<double,4> _emobm;
        maps::physmap<double,4> _hmobm;
        
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
        
    public:
        mapconfig(){}
        
        void load(std::string map_path){
            INFO_LINE("Loading single map")
            _pm.init(map_path);
            _pm.prepare();
            _pm.fillspacepoints<VecDev_t<double> >(_xvec,_yvec,_zvec);
            _pm.fillallquantities<VecDev_t<double> >(efieldvecx,efieldvecy,efieldvecz,emobvec,hmobvec,wfieldvecx,wfieldvecy,wfieldvecz);
            std::get<0>(_borders) = *(_pm.getx().begin());
            std::get<1>(_borders) = *(_pm.getx().rbegin());
            std::get<2>(_borders) = *(_pm.gety().begin());
            std::get<3>(_borders) = *(_pm.gety().rbegin());
            std::get<4>(_borders) = *(_pm.getz().begin());
            std::get<5>(_borders) = *(_pm.getz().rbegin());
            _multimap = false;
            std::cout << MIDDLEROW << std::endl;
            INFO_LINE("Volume boundaries: [ "<< std::get<0>(_borders) << " < x < "<< std::get<1>(_borders) << " ] micron, "<<"[ "<< std::get<2>(_borders) << " < y < "<< std::get<3>(_borders) <<" ] micron, "<<"[ "<< std::get<4>(_borders) << " < z < "<< std::get<5>(_borders) <<" ] micron")
        }
        
        void load(std::string efield_path, std::string emob_path, std::string hmob_path, std::string wfield_path){
            INFO_LINE("Loading separate maps")
            
            
            _efieldm.init(efield_path); //E field is a vector
            _efieldm.prepare();
            _efieldm.fillspacepoints<VecDev_t<double> >(_xvec_ef,_yvec_ef,_zvec_ef);
            _efieldm.fillvector<VecDev_t<double> >(efieldvecx,efieldvecy,efieldvecz);
            
//             for(auto cf:configlist){
//                 if(std::get<bool>(cf.get("plot"))){
//                     VecHost_t<double> efieldvecx_host, efieldvecy_host, efieldvecz_host;
//                     efieldm.fillvector<VecHost_t<double> >(efieldvecx_host,efieldvecy_host,efieldvecz_host);
//                     h3=analysis::getHist3DDraw( efieldvecx_host, efieldvecy_host, efieldvecz_host, efieldm.getx(),efieldm.gety(),efieldm.getz());
//                     break;
//                 }
//             }

        
            _wfieldm.init(wfield_path); // Weighting field is a vector
            _wfieldm.prepare();
            _wfieldm.fillspacepoints<VecDev_t<double> >(_xvec_wf,_yvec_wf,_zvec_wf);
            _wfieldm.fillvector<VecDev_t<double> >(wfieldvecx,wfieldvecy,wfieldvecz);
            
            
            _emobm.init(emob_path); // Electron mobility is a scalar
            _emobm.prepare();
            _emobm.fillspacepoints<VecDev_t<double> >(_xvec_emob,_yvec_emob,_zvec_emob);
            _emobm.fillscalar<VecDev_t<double> >(emobvec);
            
            
            _hmobm.init(hmob_path); // Hole mobility is a scalar
            _hmobm.prepare();
            _hmobm.fillspacepoints<VecDev_t<double> >(_xvec_hmob,_yvec_hmob,_zvec_hmob);
            _hmobm.fillscalar<VecDev_t<double> >(hmobvec);
            
            
            std::get<0>(_borders) = std::min(*(_efieldm.getx().begin()),std::min(*(_wfieldm.getx().begin()),std::min(*(_emobm.getx().begin()),*(_hmobm.getx().begin()))));
            std::get<1>(_borders) = std::max(*(_efieldm.getx().rbegin()),std::max(*(_wfieldm.getx().rbegin()),std::max(*(_emobm.getx().rbegin()),*(_hmobm.getx().rbegin()))));
            std::get<2>(_borders) = std::min(*(_efieldm.gety().begin()),std::min(*(_wfieldm.gety().begin()),std::min(*(_emobm.gety().begin()),*(_hmobm.gety().begin()))));
            std::get<3>(_borders) = std::max(*(_efieldm.gety().rbegin()),std::max(*(_wfieldm.gety().rbegin()),std::max(*(_emobm.gety().rbegin()),*(_hmobm.gety().rbegin()))));
            std::get<4>(_borders) = std::min(*(_efieldm.getz().begin()),std::min(*(_wfieldm.getz().begin()),std::min(*(_emobm.getz().begin()),*(_hmobm.getz().begin()))));
            std::get<5>(_borders) = std::max(*(_efieldm.getz().rbegin()),std::max(*(_wfieldm.getz().rbegin()),std::max(*(_emobm.getz().rbegin()),*(_hmobm.getz().rbegin()))));
            
            _multimap = true;
            
            std::cout << MIDDLEROW << std::endl;
            INFO_LINE("Volume boundaries: [ "<< std::get<0>(_borders) << " < x < "<< std::get<1>(_borders) << " ] micron, "<<"[ "<< std::get<2>(_borders) << " < y < "<< std::get<3>(_borders) <<" ] micron, "<<"[ "<< std::get<4>(_borders) << " < z < "<< std::get<5>(_borders) <<" ] micron")
        }
        
        //is multimap or not?
        bool ismulti(){return _multimap;}
    };
    
    //prepare maps
    void preparemaps(cfg::Setting const& root){
        
        //check if output directory exists, if it does not it is created
        std::string outputdir=root["OutputDirectory"];
        mkdir(outputdir.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRGRP | S_IROTH);
        
        auto start = std::chrono::high_resolution_clock::now();
        if(root.exists("pm")){
            INFO_LINE("Loading single map")
            
            physmap<double,11> pm(root["map"]["path"], (int)root["map"]["sc"], (int)root["map"]["sr"], (double)root["map"]["scale"], (double)root["map"]["scalespace"]);
            pm.prepare();
            std::string out=root["pm"]["out"];
            pm.write(outputdir+(std::string)"/"+out);
            
        }
        else{
            INFO_LINE("Loading separate maps")
            
            {
                physmap<double,6> ef(root["efield"]["path"], (int)root["efield"]["sc"], (int)root["efield"]["sr"], (double)root["efield"]["scale"], (double)root["efield"]["scalespace"]);
                ef.prepare();
                std::string out_ef=root["efield"]["out"];
                ef.write(outputdir+(std::string)"/"+out_ef);
            }
            
            {
                physmap<double,6> wf(root["wfield"]["path"], (int)root["wfield"]["sc"], (int)root["wfield"]["sr"], (double)root["wfield"]["scale"], (double)root["wfield"]["scalespace"]);
                wf.prepare();
                std::string out_wf=root["wfield"]["out"];
                wf.write(outputdir+(std::string)"/"+out_wf);
            }
            
            {
                physmap<double,4> emob(root["emob"]["path"], (int)root["emob"]["sc"], (int)root["emob"]["sr"], (double)root["emob"]["scale"], (double)root["emob"]["scalespace"]);
                emob.prepare();
                std::string out_emob=root["emob"]["out"];
                emob.write(outputdir+(std::string)"/"+out_emob);
            }
            
            {
                physmap<double,4> hmob(root["hmob"]["path"], (int)root["hmob"]["sc"], (int)root["hmob"]["sr"], (double)root["hmob"]["scale"], (double)root["hmob"]["scalespace"]);
                hmob.prepare();
                std::string out_hmob=root["hmob"]["out"];
                hmob.write(outputdir+(std::string)"/"+out_hmob);
            }
        }
        
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        INFO_LINE("Time to process the maps : "<<elapsed.count()<<" s")
        std::cout << TOPROW << std::endl;
    }
}

#endif
