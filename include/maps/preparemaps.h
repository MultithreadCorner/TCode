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
        bool _isloaded=false;
        std::string _filename="none";
        bool _sparse=false;
        size_t _cskip=0;
        size_t _rskip=0;
        bool _gridyfied=false;
        double _scalespace=1.;
        double _scalefield=1.;
        bool _has_active_flag=false;
        std::string _active_flag="";
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
        
        physmap(std::string s, size_t skipcolumns=0, size_t skiprows=0, double scale=1., double scalespace=1., bool sparse=false, std::string search_active=""){
            this->init(s,skipcolumns,skiprows,scale,scalespace,sparse,search_active);
        }
        
        void init(std::string s, size_t skipcolumns=0, size_t skiprows=0, double scale=1., double scalespace=1., bool sparse=false, std::string search_active=""){
            _filename=s;
            _cskip=skipcolumns;
            _rskip=skiprows;
            _sparse=false;
            _gridyfied=false;
            _scalefield=scale;
            _scalespace=scalespace;
            if(search_active!=""){
                _has_active_flag = true;
                _active_flag = search_active;
            }
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
            _isloaded = true;
        }
        
        void set_scalespace(double A){_scalespace=A;}
        
        void set_scale(double A){_scalefield=A;}
        
        double get_scalespace(double A){return _scalespace;}
        
        double get_scale(double A){return _scalefield;}
        
        void gridify(std::set<MP> xs, std::set<MP> ys, std::set<MP> zs); //get in a grid
        
        void write(std::string outputpath);
        
        void getHist3D(TH3D *h3, size_t tp=0);
        
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
        
        
        template <typename Container, typename ElemType = decltype(*(std::declval<Container>().begin()))>
        std::enable_if_t< hydra::tuple_size< ElemType >::value == 2 && N==5, void >
        fillscalar_and_active(Container &v){
            v.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_3,_4),_data.end(_3,_4)),v);
        }
        
        template <typename Container>
        std::enable_if<N==11 || N==12, void >
        fillemobility(Container &v){
            v.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_6),_data.end(_6)),v);
        }
        
        template <typename Container, typename ElemType = decltype(*(std::declval<Container>().begin()))>
        std::enable_if_t< hydra::tuple_size< ElemType >::value == 2 && N==12, void >
        fillemobility_and_active(Container &v){
            v.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_6,_11),_data.end(_6,_11)),v);
        }
        
        template <typename Container>
        std::enable_if<N==11 || N==12, void >
        fillhmobility(Container &v){
            v.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_7),_data.end(_7)),v);
        }
        
        template <typename Container, typename ElemType = decltype(*(std::declval<Container>().begin()))>
        std::enable_if_t< hydra::tuple_size< ElemType >::value == 2 && N==11, void >
        fillhmobility_and_active(Container &v){
            v.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_7,_11),_data.end(_7,_11)),v);
        }
        
        template <typename Container, typename ElemType = decltype(*(std::declval<Container>().begin()))>
        std::enable_if_t< hydra::tuple_size< ElemType >::value == 3 && (N==6 || N==7), void >
        fillvector(Container &vect){
            vect.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_3,_4,_5),_data.end(_3,_4,_5)),vect);
        }
        
        template <typename Container, typename ElemType = decltype(*(std::declval<Container>().begin()))>
        std::enable_if_t< hydra::tuple_size< ElemType >::value == 4 && N==7, void >
        fillvector_and_active(Container &vect){
            vect.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_3,_4,_5,_6),_data.end(_3,_4,_5,_6)),vect);
        }
        
        template <typename Container, typename ElemType = decltype(*(std::declval<Container>().begin()))>
        std::enable_if_t< hydra::tuple_size< ElemType >::value == 3 && (N==11 || N==12), void >
        fillvectorEfield(Container &vect){
            vect.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_3,_4,_5),_data.end(_3,_4,_5)),vect);
        }
        
        template <typename Container, typename ElemType = decltype(*(std::declval<Container>().begin()))>
        std::enable_if_t< hydra::tuple_size< ElemType >::value == 4 && N==12, void >
        fillvectorEfield_and_active(Container &vect){
            vect.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_3,_4,_5,_11),_data.end(_3,_4,_5,_11)),vect);
        }
        
        template <typename Container, typename ElemType = decltype(*(std::declval<Container>().begin()))>
        std::enable_if_t< hydra::tuple_size< ElemType >::value == 3 && (N==11 || N==12), void >
        fillvectorWfield(Container &vect){
            vect.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_8,_9,_10),_data.end(_8,_9,_10)),vect);
        }
        
        template <typename Container, typename ElemType = decltype(*(std::declval<Container>().begin()))>
        std::enable_if_t< hydra::tuple_size< ElemType >::value == 4 && N==12, void >
        fillvectorWfield_and_active(Container &vect){
            vect.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_8,_9,_10,_11),_data.end(_8,_9,_10,_11)),vect);
        }
        
        template <typename Container, typename ElemType = decltype(*(std::declval<Container>().begin()))>
        std::enable_if_t< hydra::tuple_size< ElemType >::value == 8 && (N==11 || N==12), void >
        fillallquantities(Container &fields){
            fields.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_3,_4,_5,_6,_7,_8,_9,_10),_data.end(_3,_4,_5,_6,_7,_8,_9,_10)),fields);
        }
            
        template <typename Container, typename ElemType = decltype(*(std::declval<Container>().begin()))>
        std::enable_if_t< hydra::tuple_size< ElemType >::value == 9 && N==12, void >
        fillallquantities_and_active(Container &fields){
            fields.resize(_data.size());
            hydra::copy(hydra::make_range(_data.begin(_3,_4,_5,_6,_7,_8,_9,_10,_11),_data.end(_3,_4,_5,_6,_7,_8,_9,_10,_11)),fields);
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
    
    
#ifdef _ROOT_AVAILABLE_
    template<typename MP, size_t N>
    void physmap<MP,N>::getHist3D(TH3D *h3, size_t tp){
        
        if(!_isloaded){
            ERROR_LINE("Load map before asking for the electric field histogram!");
            exit(1);
        }
        
        std::vector<double> vx(_gridx.begin(),_gridx.end());
        std::vector<double> vy(_gridy.begin(),_gridy.end());
        std::vector<double> vz(_gridz.begin(),_gridz.end());
        std::vector<double> vex(vx.size()+1);
        std::vector<double> vey(vy.size()+1);
        std::vector<double> vez(vz.size()+1);
        
        VectorMapHost_t _map(_data.size());
        
        std::string nameh3 = "efield";
        std::string titleh3 = "Electric Field [V #mum^{-1}];x [#mum];y [#mum];z [#mum]";
        
        if(N<6){
            ERROR_LINE("Cannot get 3D histogram from this map!")
            exit(1);
        }
        
        std::vector<double> _mdata(_data.size());
        
        if(tp==0 || N==6) hydra::transform(hydra::make_range(_data.begin(_3,_4,_5),_data.end(_3,_4,_5)), _mdata, vectorLength());
        
        if(tp==1){
            if(N<6){
                ERROR_LINE("Cannot get 3D weigthing field histogram from this map!")
                exit(1);
            }
            hydra::transform(hydra::make_range(_data.begin(_3,_4,_5),_data.end(_3,_4,_5)), _mdata, vectorLength());
            nameh3 = "wfield";
            titleh3 = "Weighting Field [#mum^{-1}];x [#mum];y [#mum];z [#mum]";
        }
        
        vex[0]=vx[0]-fabs(vx[0]-vx[1])/10.;
        vey[0]=vy[0]-fabs(vy[0]-vy[1])/10.;
        vez[0]=vz[0]-fabs(vz[0]-vz[1])/10.;
        
        for(size_t ix=0;ix<vx.size()-1;ix++){
            vex[ix+1]=(vx[ix]+vx[ix+1])/2.;
            
        }
        for(size_t iy=0;iy<vy.size()-1;iy++){
            vey[iy+1]=(vy[iy]+vy[iy+1])/2.;
            
        }
        for(size_t iz=0;iz<vz.size()-1;iz++){
            vez[iz+1]=(vz[iz]+vz[iz+1])/2.;
        }
        
        vex[vx.size()]=vx[vx.size()-1]+fabs(vx[vx.size()-2]-vx[vx.size()-1])/10.;
        vey[vy.size()]=vy[vy.size()-1]+fabs(vy[vy.size()-2]-vy[vy.size()-1])/10.;
        vez[vz.size()]=vz[vz.size()-1]+fabs(vz[vz.size()-2]-vz[vz.size()-1])/10.;
        
        
        h3->SetName(nameh3.c_str());
        h3->SetTitle(titleh3.c_str());
        h3->SetBins(vx.size(),vex.data(),vy.size(),vey.data(),vz.size(),vez.data());

//         for(auto& n: _mdata) std::cout << n << std::endl;
        size_t count=0;
        for(size_t iz=0;iz<vz.size();iz++){
            for(size_t ix=0;ix<vx.size();ix++){
                for(size_t iy=0;iy<vy.size();iy++){
                    h3->SetBinContent(ix+1,iy+1,iz+1, _mdata[count]);
                    count++;
                }
            }
        }
    }
#endif //_ROOT_AVAILABLE_
    
    
    //class to manage map configuration
    class mapconfig{
    private:
        bool _multimap = false;
        bool _isloaded = false;
        bool _isondevice = false;
        std::tuple<double,double,double,double,double,double> _borders;
        
        maps::physmap<double,11> _pm;
        maps::physmap<double,6> _efieldm;
        maps::physmap<double,6> _wfieldm;
        maps::physmap<double,4> _emobm;
        maps::physmap<double,4> _hmobm;
        
        VecDev_t<double> _xvec, _yvec, _zvec;
        PhysMapDev_t _vec_map;
        
        //efield vecs
        VecDev_t<double> _xvec_ef, _yvec_ef, _zvec_ef;
        VectorMapDev_t _vec_ef;
        
        //wfield vecs
        VecDev_t<double> _xvec_wf, _yvec_wf, _zvec_wf;
        VectorMapDev_t _vec_wf;
        
        //emob vecs
        VecDev_t<double> _xvec_emob, _yvec_emob, _zvec_emob;
        VecDev_t<double> _emobvec;
        
        //hmob vecs
        VecDev_t<double> _xvec_hmob, _yvec_hmob, _zvec_hmob;
        VecDev_t<double> _hmobvec;
        
    public:
        mapconfig() {
            _isloaded = false;
        };
        
        mapconfig(std::string map_path) {
            this->load(map_path);
        }
        
        mapconfig(cfg::Setting const& root){
            if(root["PhysicsMaps"].exists("map")) {
                INFO_LINE("Loading single map")
                this->load(root["PhysicsMaps"]["map"]);
            }
            else{
                if(root["PhysicsMaps"].exists("efield") && root["PhysicsMaps"].exists("wfield") && root["PhysicsMaps"].exists("emob") && root["PhysicsMaps"].exists("hmob")){
                    INFO_LINE("Loading separate maps")
                    this->load(root["PhysicsMaps"]["efield"],root["PhysicsMaps"]["emob"],root["PhysicsMaps"]["hmob"],root["PhysicsMaps"]["wfield"]);
                }
                else{
                    ERROR_LINE("Not all the maps exist, check your configuration!")
                    exit(1);
                }
            }
        }
        
        mapconfig(std::string efield_path, std::string emob_path, std::string hmob_path, std::string wfield_path){
            this->load(efield_path, emob_path, hmob_path, wfield_path);
        }
        
        void load(std::string map_path){
            _pm.init(map_path);
            _pm.prepare();
            std::get<0>(_borders) = *(_pm.getx().begin());
            std::get<1>(_borders) = *(_pm.getx().rbegin());
            std::get<2>(_borders) = *(_pm.gety().begin());
            std::get<3>(_borders) = *(_pm.gety().rbegin());
            std::get<4>(_borders) = *(_pm.getz().begin());
            std::get<5>(_borders) = *(_pm.getz().rbegin());
            _multimap = false;
            std::cout << MIDDLEROW << std::endl;
            INFO_LINE("Volume boundaries: [ "<< std::get<0>(_borders) << " < x < "<< std::get<1>(_borders) << " ] micron, "<<"[ "<< std::get<2>(_borders) << " < y < "<< std::get<3>(_borders) <<" ] micron, "<<"[ "<< std::get<4>(_borders) << " < z < "<< std::get<5>(_borders) <<" ] micron")
            _isloaded = true;
        }
        
        void load(std::string efield_path, std::string emob_path, std::string hmob_path, std::string wfield_path){
            //get efield
            _efieldm.init(efield_path); //E field is a vector
            _efieldm.prepare();
        
            _wfieldm.init(wfield_path); // Weighting field is a vector
            _wfieldm.prepare();
            
            _emobm.init(emob_path); // Electron mobility is a scalar
            _emobm.prepare();
            
            _hmobm.init(hmob_path); // Hole mobility is a scalar
            _hmobm.prepare();
            
            std::get<0>(_borders) = std::min(*(_efieldm.getx().begin()),std::min(*(_wfieldm.getx().begin()),std::min(*(_emobm.getx().begin()),*(_hmobm.getx().begin()))));
            std::get<1>(_borders) = std::max(*(_efieldm.getx().rbegin()),std::max(*(_wfieldm.getx().rbegin()),std::max(*(_emobm.getx().rbegin()),*(_hmobm.getx().rbegin()))));
            std::get<2>(_borders) = std::min(*(_efieldm.gety().begin()),std::min(*(_wfieldm.gety().begin()),std::min(*(_emobm.gety().begin()),*(_hmobm.gety().begin()))));
            std::get<3>(_borders) = std::max(*(_efieldm.gety().rbegin()),std::max(*(_wfieldm.gety().rbegin()),std::max(*(_emobm.gety().rbegin()),*(_hmobm.gety().rbegin()))));
            std::get<4>(_borders) = std::min(*(_efieldm.getz().begin()),std::min(*(_wfieldm.getz().begin()),std::min(*(_emobm.getz().begin()),*(_hmobm.getz().begin()))));
            std::get<5>(_borders) = std::max(*(_efieldm.getz().rbegin()),std::max(*(_wfieldm.getz().rbegin()),std::max(*(_emobm.getz().rbegin()),*(_hmobm.getz().rbegin()))));
            
            _multimap = true;
            
            std::cout << MIDDLEROW << std::endl;
            INFO_LINE("Volume boundaries: [ "<< std::get<0>(_borders) << " < x < "<< std::get<1>(_borders) << " ] micron, "<<"[ "<< std::get<2>(_borders) << " < y < "<< std::get<3>(_borders) <<" ] micron, "<<"[ "<< std::get<4>(_borders) << " < z < "<< std::get<5>(_borders) <<" ] micron")
            
            _isloaded = true;
        }
        
        //is multimap or not?
        bool ismulti(){return _multimap;}
        bool isloaded(){return _isloaded;}
        
        void todevice(){
            if(_multimap){
                _efieldm.fillspacepoints(_xvec_ef,_yvec_ef,_zvec_ef);
                _efieldm.fillvector(_vec_ef);
                
                _emobm.fillspacepoints(_xvec_emob,_yvec_emob,_zvec_emob);
                _emobm.fillscalar(_emobvec);
                
                _hmobm.fillspacepoints(_xvec_hmob,_yvec_hmob,_zvec_hmob);
                _hmobm.fillscalar(_hmobvec);
                
                _wfieldm.fillspacepoints(_xvec_wf,_yvec_wf,_zvec_wf);
                _wfieldm.fillvector(_vec_wf);
            }
            else{
                _pm.fillspacepoints(_xvec,_yvec,_zvec);
                _pm.fillallquantities(_vec_map);
            }
            _isondevice = true;
        }
        
        //getters for maps
        PhysMapDev_t& PhysMap()   { 
            if(!_isondevice) this->todevice();
            return _vec_map; 
        }
        
        VecDev_t<double>& VecX()   { 
            if(!_isondevice) this->todevice();
            return _xvec; 
        }
        
        VecDev_t<double>& VecY()   { 
            if(!_isondevice) this->todevice();
            return _yvec; 
        }
        
        VecDev_t<double>& VecZ()   { 
            if(!_isondevice) this->todevice();
            return _zvec; 
        }
        
        VectorMapDev_t& EField()   { 
            if(!_isondevice) this->todevice();
            return _vec_ef; 
        }
        
        VecDev_t<double>& VecXEField()   { 
            if(!_isondevice) this->todevice();
            return _xvec_ef; 
        }
        
        VecDev_t<double>& VecYEField()   { 
            if(!_isondevice) this->todevice();
            return _yvec_ef; 
        }
        
        VecDev_t<double>& VecZEField()   { 
            if(!_isondevice) this->todevice();
            return _zvec_ef; 
        }
        
        VectorMapDev_t& WField()   { 
            if(!_isondevice) this->todevice();
            return _vec_wf; 
        }
        
        VecDev_t<double>& VecXWField()   { 
            if(!_isondevice) this->todevice();
            return _xvec_wf; 
        }
        
        VecDev_t<double>& VecYWField()   { 
            if(!_isondevice) this->todevice();
            return _yvec_wf; 
        }
        
        VecDev_t<double>& VecZWField()   { 
            if(!_isondevice) this->todevice();
            return _zvec_wf; 
        }
        
        VecDev_t<double>& EMob()   { 
            if(!_isondevice) this->todevice();
            return _emobvec; 
        }
        
        VecDev_t<double>& VecXEMob()   { 
            if(!_isondevice) this->todevice();
            return _xvec_emob; 
        }
        
        VecDev_t<double>& VecYEMob()   { 
            if(!_isondevice) this->todevice();
            return _yvec_emob; 
        }
        
        VecDev_t<double>& VecZEMob()   { 
            if(!_isondevice) this->todevice();
            return _zvec_emob; 
        }
        
        VecDev_t<double>& HMob()   { 
            if(!_isondevice) this->todevice();
            return _hmobvec; 
        }
        
        VecDev_t<double>& VecXHMob()   { 
            if(!_isondevice) this->todevice();
            return _xvec_hmob; 
        }
        
        VecDev_t<double>& VecYHMob()   { 
            if(!_isondevice) this->todevice();
            return _yvec_hmob; 
        }
        
        VecDev_t<double>& VecZHMob()   { 
            if(!_isondevice) this->todevice();
            return _zvec_hmob; 
        }
        
        std::tuple<double,double,double,double,double,double>& Boundaries(){
            if(!_isloaded){
                ERROR_LINE("The map(s)  was not loaded!")
                exit(1);
            }
            return _borders;
        }
        
        #ifdef _ROOT_AVAILABLE_
        void getHist3D(TH3D *h3, size_t tp=0){
            
            if(!_isloaded){
                ERROR_LINE("Load map before asking for the electric field histogram!");
                exit(1);
            }
            
            if(_multimap) _efieldm.getHist3D(h3,tp);
            else _pm.getHist3D(h3,tp);
        }
        #endif //_ROOT_AVAILABLE_
        
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
