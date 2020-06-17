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
 * Evolve.h
 *
 *  Created on: 12/11/2018
 *  Author: Andrea Contu
 */

#ifndef __EVOLVE_H__
#define __EVOLVE_H__

namespace evolve{    
    
    
    //functions to calculates trilinear interpolation
    __hydra_dual__
    inline double cfun(double d, double t1, double t2){
        return t1*(1.-d)+t2*d;
    }
    
    __hydra_dual__
    inline double val(double xd, double yd, double zd, double p0, double px, double py, double pz, double pxy, double pyz, double pxz, double pxyz){     
        return cfun(zd,cfun(yd, cfun(xd,p0,px),cfun(xd,py,pxy)), cfun(yd,cfun(xd,pz,pxz), cfun(xd,pyz,pxyz)));
    }
    
    
    
    //function to get index in one dimension (binary search)
    template<typename Iterator>
    __hydra_dual__
    inline size_t getindexIT(Iterator v, size_t n, double f, size_t first=0){
        if(f<= *(v+first)){
            return first+1;

        }
        if(f>= *(v+n-1)){
            return n-1;
        }
        size_t lo = first;
        size_t hi = n - 1;
        size_t mid=0;
        while(lo<hi){

             mid=(hi+lo)/2;
             
             if( *(v+mid)==f) return mid+1;

             if(f< *(v+mid)){
                 if(mid>0 && f> *(v+mid-1)) return mid;
                 hi=mid;
             }
             else{
                 if(mid<n-1 && f < *(v+mid+1)) return mid+1;
                 lo=mid+1;
             }
         }
         return 0;
    }
    
    //function to get points indices in a given cartesian dimension
    template<typename Iterator>
    __hydra_dual__
    inline hydra::tuple<size_t,size_t,size_t> getindicesIT(Iterator vx, size_t sx, Iterator vy, size_t sy, Iterator vz, size_t sz, hydra::Vector3R &pos){
        
        return hydra::make_tuple(
            getindexIT(vx, sx, pos.get(0)),
            getindexIT(vy, sy, pos.get(1)),
            getindexIT(vz, sz, pos.get(2))
        );
    }
    
    //function to get vectors and indices in physics tables
    template<typename Iterator>
    __hydra_dual__
    inline mapindices_t getallindicesIT(hydra::Vector3R &vvec, Iterator vx, size_t sx, Iterator vy, size_t sy, Iterator vz, size_t sz, hydra::Vector3R &pos){
        hydra::tuple<size_t,size_t,size_t> ind = getindicesIT(vx, sx, vy, sy, vz, sz, pos);
        vvec.set(   (pos.get(0)- *(vx+hydra::get<0>(ind)-1))/( *(vx+hydra::get<0>(ind))- *(vx+hydra::get<0>(ind)-1)),
                    (pos.get(1)- *(vy+hydra::get<1>(ind)-1))/( *(vy+hydra::get<1>(ind))- *(vy+hydra::get<1>(ind)-1)),
                    (pos.get(2)- *(vz+hydra::get<2>(ind)-1))/( *(vz+hydra::get<2>(ind))- *(vz+hydra::get<2>(ind)-1))
                );
        
        
        double fNPXY=sx*sy;
        
        return hydra::make_tuple(
            hydra::get<1>(ind)-1 + (hydra::get<0>(ind)-1)*sy+(hydra::get<2>(ind)-1)*fNPXY,
            hydra::get<1>(ind)-1 + (hydra::get<0>(ind))*sy+(hydra::get<2>(ind)-1)*fNPXY,
            hydra::get<1>(ind) + (hydra::get<0>(ind)-1)*sy+(hydra::get<2>(ind)-1)*fNPXY,
            hydra::get<1>(ind)-1 + (hydra::get<0>(ind)-1)*sy+(hydra::get<2>(ind))*fNPXY,
            hydra::get<1>(ind) + (hydra::get<0>(ind))*sy+(hydra::get<2>(ind)-1)*fNPXY,
            hydra::get<1>(ind) + (hydra::get<0>(ind)-1)*sy+(hydra::get<2>(ind))*fNPXY,
            hydra::get<1>(ind)-1 +(hydra::get<0>(ind))*sy+(hydra::get<2>(ind))*fNPXY,
            hydra::get<1>(ind) + (hydra::get<0>(ind))*sy+(hydra::get<2>(ind))*fNPXY
        );
    }
    
    template<typename Iterator>
    __hydra_dual__
    inline double interpolIT(hydra::Vector3R &vvec, Iterator v, mapindices_t &ind8){     
        return val(vvec.get(0),vvec.get(1), vvec.get(2), *(v+hydra::get<0>(ind8)), *(v+hydra::get<1>(ind8)), *(v+hydra::get<2>(ind8)), *(v+hydra::get<3>(ind8)), *(v+hydra::get<4>(ind8)), *(v+hydra::get<5>(ind8)), *(v+hydra::get<6>(ind8)), *(v+hydra::get<7>(ind8)));
    }
    
    template<typename Iterator>
    __hydra_dual__
    inline hydra::Vector3R interpolITvec(hydra::Vector3R &vvec, Iterator v, mapindices_t &ind8){     
        
        return hydra::Vector3R(val(vvec.get(0),vvec.get(1), vvec.get(2), hydra::get<0>(*(v+hydra::get<0>(ind8))), hydra::get<0>(*(v+hydra::get<1>(ind8))), hydra::get<0>(*(v+hydra::get<2>(ind8))), hydra::get<0>(*(v+hydra::get<3>(ind8))), hydra::get<0>(*(v+hydra::get<4>(ind8))),         hydra::get<0>(*(v+hydra::get<5>(ind8))), hydra::get<0>(*(v+hydra::get<6>(ind8))), hydra::get<0>(*(v+hydra::get<7>(ind8)))), 
                               val(vvec.get(0),vvec.get(1), vvec.get(2), hydra::get<1>(*(v+hydra::get<0>(ind8))), hydra::get<1>(*(v+hydra::get<1>(ind8))), hydra::get<1>(*(v+hydra::get<2>(ind8))), hydra::get<1>(*(v+hydra::get<3>(ind8))), hydra::get<1>(*(v+hydra::get<4>(ind8))), hydra::get<1>(*(v+hydra::get<5>(ind8))), hydra::get<1>(*(v+hydra::get<6>(ind8))), hydra::get<1>(*(v+hydra::get<7>(ind8)))),
                                val(vvec.get(0),vvec.get(1), vvec.get(2), hydra::get<2>(*(v+hydra::get<0>(ind8))), hydra::get<2>(*(v+hydra::get<1>(ind8))), hydra::get<2>(*(v+hydra::get<2>(ind8))), hydra::get<2>(*(v+hydra::get<3>(ind8))), hydra::get<2>(*(v+hydra::get<4>(ind8))), hydra::get<2>(*(v+hydra::get<5>(ind8))), hydra::get<2>(*(v+hydra::get<6>(ind8))), hydra::get<2>(*(v+hydra::get<7>(ind8)))));
        
    }
    
    template<typename Iterator>
    __hydra_dual__
    inline mapentry_t interpolITtuple(hydra::Vector3R &vvec, Iterator v, mapindices_t &ind8){     
        
        return hydra::make_tuple(val(vvec.get(0),vvec.get(1), vvec.get(2), hydra::get<0>(*(v+hydra::get<0>(ind8))), hydra::get<0>(*(v+hydra::get<1>(ind8))), hydra::get<0>(*(v+hydra::get<2>(ind8))), hydra::get<0>(*(v+hydra::get<3>(ind8))), hydra::get<0>(*(v+hydra::get<4>(ind8))),         hydra::get<0>(*(v+hydra::get<5>(ind8))), hydra::get<0>(*(v+hydra::get<6>(ind8))), hydra::get<0>(*(v+hydra::get<7>(ind8)))), 
                               val(vvec.get(0),vvec.get(1), vvec.get(2), hydra::get<1>(*(v+hydra::get<0>(ind8))), hydra::get<1>(*(v+hydra::get<1>(ind8))), hydra::get<1>(*(v+hydra::get<2>(ind8))), hydra::get<1>(*(v+hydra::get<3>(ind8))), hydra::get<1>(*(v+hydra::get<4>(ind8))), hydra::get<1>(*(v+hydra::get<5>(ind8))), hydra::get<1>(*(v+hydra::get<6>(ind8))), hydra::get<1>(*(v+hydra::get<7>(ind8)))),
                                val(vvec.get(0),vvec.get(1), vvec.get(2), hydra::get<2>(*(v+hydra::get<0>(ind8))), hydra::get<2>(*(v+hydra::get<1>(ind8))), hydra::get<2>(*(v+hydra::get<2>(ind8))), hydra::get<2>(*(v+hydra::get<3>(ind8))), hydra::get<2>(*(v+hydra::get<4>(ind8))), hydra::get<2>(*(v+hydra::get<5>(ind8))), hydra::get<2>(*(v+hydra::get<6>(ind8))), hydra::get<2>(*(v+hydra::get<7>(ind8)))),
                                  val(vvec.get(0),vvec.get(1), vvec.get(2), hydra::get<3>(*(v+hydra::get<0>(ind8))), hydra::get<3>(*(v+hydra::get<1>(ind8))), hydra::get<3>(*(v+hydra::get<2>(ind8))), hydra::get<3>(*(v+hydra::get<3>(ind8))), hydra::get<3>(*(v+hydra::get<4>(ind8))), hydra::get<3>(*(v+hydra::get<5>(ind8))), hydra::get<3>(*(v+hydra::get<6>(ind8))), hydra::get<3>(*(v+hydra::get<7>(ind8)))),
                                  val(vvec.get(0),vvec.get(1), vvec.get(2), hydra::get<4>(*(v+hydra::get<0>(ind8))), hydra::get<4>(*(v+hydra::get<1>(ind8))), hydra::get<4>(*(v+hydra::get<2>(ind8))), hydra::get<4>(*(v+hydra::get<3>(ind8))), hydra::get<4>(*(v+hydra::get<4>(ind8))), hydra::get<4>(*(v+hydra::get<5>(ind8))), hydra::get<4>(*(v+hydra::get<6>(ind8))), hydra::get<4>(*(v+hydra::get<7>(ind8)))),
                                  val(vvec.get(0),vvec.get(1), vvec.get(2), hydra::get<5>(*(v+hydra::get<0>(ind8))), hydra::get<5>(*(v+hydra::get<1>(ind8))), hydra::get<5>(*(v+hydra::get<2>(ind8))), hydra::get<5>(*(v+hydra::get<3>(ind8))), hydra::get<5>(*(v+hydra::get<4>(ind8))), hydra::get<5>(*(v+hydra::get<5>(ind8))), hydra::get<5>(*(v+hydra::get<6>(ind8))), hydra::get<5>(*(v+hydra::get<7>(ind8)))),
                                  val(vvec.get(0),vvec.get(1), vvec.get(2), hydra::get<6>(*(v+hydra::get<0>(ind8))), hydra::get<6>(*(v+hydra::get<1>(ind8))), hydra::get<6>(*(v+hydra::get<2>(ind8))), hydra::get<6>(*(v+hydra::get<3>(ind8))), hydra::get<6>(*(v+hydra::get<4>(ind8))), hydra::get<6>(*(v+hydra::get<5>(ind8))), hydra::get<6>(*(v+hydra::get<6>(ind8))), hydra::get<6>(*(v+hydra::get<7>(ind8)))),
                                val(vvec.get(0),vvec.get(1), vvec.get(2), hydra::get<7>(*(v+hydra::get<0>(ind8))), hydra::get<7>(*(v+hydra::get<1>(ind8))), hydra::get<7>(*(v+hydra::get<2>(ind8))), hydra::get<7>(*(v+hydra::get<3>(ind8))), hydra::get<7>(*(v+hydra::get<4>(ind8))), hydra::get<7>(*(v+hydra::get<5>(ind8))), hydra::get<7>(*(v+hydra::get<6>(ind8))), hydra::get<7>(*(v+hydra::get<7>(ind8))))
                                
            
        );
        
    }
    
    
    //calculate current with Ramo's theorem
    template<typename Iteratordim, typename Iteratormap>
    struct RamoCurrent
    {

        RamoCurrent()=delete;                                         // avoid creating objects with undefined data

        RamoCurrent(size_t n, double timestep, double diffusion, Iteratordim begindimx, size_t sizedimx, Iteratordim begindimy, size_t sizedimy, Iteratordim begindimz, size_t sizedimz, Iteratormap beginmap, size_t sizemap):
        _fN(n),
        _fStep(timestep),
        _fDiff(diffusion),
        _begindimx(begindimx),
        _sizedimx(sizedimx),
        _begindimy(begindimy),
        _sizedimy(sizedimy),
        _begindimz(begindimz),
        _sizedimz(sizedimz),
        _beginmap(beginmap),
        _sizemap(sizemap)
        {}
        
        // Copy constructor needs be __hydra_dual__ to be copyable
        // on device and host sides, in the time of thread distribution.
        // hydra, as well as thrust and std pass arguments by value.
        // Obs.: RooFit passes by reference... (crazy design choice!)
        __hydra_dual__
        RamoCurrent( RamoCurrent const& other):
        _fN(other._fN),
        _fStep(other._fStep),
        _fDiff(other._fDiff),
        _begindimx(other._begindimx),
        _sizedimx(other._sizedimx),
        _begindimy(other._begindimy),
        _sizedimy(other._sizedimy),
        _begindimz(other._begindimz),
        _sizedimz(other._sizedimz),
        _beginmap(other._beginmap),
        _sizemap(other._sizemap)
        {}
        
        
        //function call operator needs to be callable in host and device sides 
        template<typename Particle>
        __hydra_dual__
        void operator()(Particle p)
        {
            //check if particle went outside the volume
            if(hydra::get<_tc_isin>(p)<0) return;
            
            //get particle charge and initial position
            double charge=sgn(hydra::get<_tc_charge>(p));
            double mobility=0;
            hydra::Vector3R position(hydra::get<_tc_x>(p), hydra::get<_tc_y>(p), hydra::get<_tc_z>(p));
            
            
            hydra::Vector3R vvec(0.,0.,0.);
            mapindices_t ind8;
            mapentry_t entry8;
            
            ind8 = getallindicesIT(vvec, _begindimx, _sizedimx, _begindimy, _sizedimy, _begindimz, _sizedimz, position);
            entry8=interpolITtuple(vvec, _beginmap, ind8);
            
            if(charge<0.){
                mobility = hydra::get<3>(entry8);
            }
            else{
                mobility = hydra::get<4>(entry8);
            }
            //get electric field vector
            hydra::Vector3R efvec(hydra::get<0>(entry8),hydra::get<1>(entry8),hydra::get<2>(entry8));
            
            //calculate velocity vector
            auto k1 =charge*mobility*efvec;
            
            //get weighting field
            hydra::Vector3R wfvec(hydra::get<5>(entry8),hydra::get<6>(entry8),hydra::get<7>(entry8));
            
            //calculate diffusion velocity
            //
            //
            double s = sqrt(2*(_fDiff*mobility)/_fStep/fabs(hydra::get<_tc_charge>(p)))*hydra::get<_tc_gauss_x>(p); //divide by square root of group size
            
            hydra::Vector3R ref(s*sin(hydra::get<_tc_angle_2>(p))*cos(hydra::get<_tc_angle_1>(p)),s*sin(hydra::get<_tc_angle_2>(p))*sin(hydra::get<_tc_angle_1>(p)),s*cos(hydra::get<_tc_angle_2>(p)));
            
            //set diffusion velocity to be used in evolution step
            hydra::get<_tc_gauss_x>(p) = ref.get(0);
            hydra::get<_tc_gauss_y>(p) = ref.get(1);
            hydra::get<_tc_gauss_z>(p) = ref.get(2);
            
            //calculate total velocity as sum of drift and diffusion
            auto Vdrift = k1 +ref;
            
            //calculate current with Ramo's theorem
            hydra::get<5>(p)=hydra::get<_tc_charge>(p)*HCHARGE*(Vdrift.get(0)*wfvec.get(0)+Vdrift.get(1)*wfvec.get(1)+Vdrift.get(2)*wfvec.get(2));
        }
        
        size_t _fN;
        double _fStep;
        double _fDiff;
        Iteratordim _begindimx;
        size_t _sizedimx;
        Iteratordim _begindimy;
        size_t _sizedimy;
        Iteratordim _begindimz;
        size_t _sizedimz;
        Iteratormap _beginmap;
        size_t _sizemap;

    };
    
    
    template<typename Iteratordim, typename Iteratorvec, typename Iteratorscal>
    struct RamoCurrent_multi
    {

        RamoCurrent_multi()=delete;                                         // avoid creating objects with undefined data

        RamoCurrent_multi(size_t n, double timestep, double diffusion, 
                          Iteratordim begindimx_ef, size_t sizedimx_ef, Iteratordim begindimy_ef, size_t sizedimy_ef, Iteratordim begindimz_ef, size_t sizedimz_ef, 
                          Iteratorvec beginmap_ef, size_t sizemap_ef,
                          Iteratordim begindimx_em, size_t sizedimx_em, Iteratordim begindimy_em, size_t sizedimy_em, Iteratordim begindimz_em, size_t sizedimz_em, 
                          Iteratorscal beginmap_em, size_t sizemap_em,
                          Iteratordim begindimx_hm, size_t sizedimx_hm, Iteratordim begindimy_hm, size_t sizedimy_hm, Iteratordim begindimz_hm, size_t sizedimz_hm, 
                          Iteratorscal beginmap_hm, size_t sizemap_hm,
                          Iteratordim begindimx_wf, size_t sizedimx_wf, Iteratordim begindimy_wf, size_t sizedimy_wf, Iteratordim begindimz_wf, size_t sizedimz_wf, 
                          Iteratorvec beginmap_wf, size_t sizemap_wf
                         ):
        _fN(n), _fStep(timestep), _fDiff(diffusion),
        _begindimx_ef(begindimx_ef), _sizedimx_ef(sizedimx_ef),
        _begindimy_ef(begindimy_ef), _sizedimy_ef(sizedimy_ef),
        _begindimz_ef(begindimz_ef), _sizedimz_ef(sizedimz_ef),
        _beginmap_ef(beginmap_ef), _sizemap_ef(sizemap_ef),
        _begindimx_em(begindimx_em), _sizedimx_em(sizedimx_em),
        _begindimy_em(begindimy_em), _sizedimy_em(sizedimy_em),
        _begindimz_em(begindimz_em), _sizedimz_em(sizedimz_em),
        _beginmap_em(beginmap_em), _sizemap_em(sizemap_em),
        _begindimx_hm(begindimx_hm), _sizedimx_hm(sizedimx_hm),
        _begindimy_hm(begindimy_hm), _sizedimy_hm(sizedimy_hm),
        _begindimz_hm(begindimz_hm), _sizedimz_hm(sizedimz_hm),
        _beginmap_hm(beginmap_hm), _sizemap_hm(sizemap_hm),
        _begindimx_wf(begindimx_wf), _sizedimx_wf(sizedimx_wf),
        _begindimy_wf(begindimy_wf), _sizedimy_wf(sizedimy_wf),
        _begindimz_wf(begindimz_wf), _sizedimz_wf(sizedimz_wf),
        _beginmap_wf(beginmap_wf), _sizemap_wf(sizemap_wf)
        {}
        
        // Copy constructor needs be __hydra_dual__ to be copyable
        // on device and host sides, in the time of thread distribution.
        // hydra, as well as thrust and std pass arguments by value.
        // Obs.: RooFit passes by reference... (crazy design choice!)
        __hydra_dual__
        RamoCurrent_multi( RamoCurrent_multi const& other):
        _fN(other._fN),  _fStep(other._fStep), _fDiff(other._fDiff),
        _begindimx_ef(other._begindimx_ef), _sizedimx_ef(other._sizedimx_ef),
        _begindimy_ef(other._begindimy_ef), _sizedimy_ef(other._sizedimy_ef),
        _begindimz_ef(other._begindimz_ef), _sizedimz_ef(other._sizedimz_ef),
        _beginmap_ef(other._beginmap_ef), _sizemap_ef(other._sizemap_ef),
        _begindimx_em(other._begindimx_em), _sizedimx_em(other._sizedimx_em),
        _begindimy_em(other._begindimy_em), _sizedimy_em(other._sizedimy_em),
        _begindimz_em(other._begindimz_em), _sizedimz_em(other._sizedimz_em),
        _beginmap_em(other._beginmap_em), _sizemap_em(other._sizemap_em),
        _begindimx_hm(other._begindimx_hm), _sizedimx_hm(other._sizedimx_hm),
        _begindimy_hm(other._begindimy_hm), _sizedimy_hm(other._sizedimy_hm),
        _begindimz_hm(other._begindimz_hm), _sizedimz_hm(other._sizedimz_hm),
        _beginmap_hm(other._beginmap_hm), _sizemap_hm(other._sizemap_hm),
        _begindimx_wf(other._begindimx_wf), _sizedimx_wf(other._sizedimx_wf),
        _begindimy_wf(other._begindimy_wf), _sizedimy_wf(other._sizedimy_wf),
        _begindimz_wf(other._begindimz_wf), _sizedimz_wf(other._sizedimz_wf),
        _beginmap_wf(other._beginmap_wf), _sizemap_wf(other._sizemap_wf)
        {}
        
        
        //function call operator needs to be callable in host and device sides 
        template<typename Particle>
        __hydra_dual__
        void operator()(Particle p)
        {        
            
            //check if particle went outside the volume
            if(hydra::get<_tc_isin>(p)<0) return;
            
            //get particle charge and initial position
            double charge=sgn(hydra::get<_tc_charge>(p));
            double mobility=0;
            hydra::Vector3R position(hydra::get<_tc_x>(p), hydra::get<_tc_y>(p), hydra::get<_tc_z>(p));
            
            
            hydra::Vector3R vvec(0.,0.,0.);
            mapindices_t ind8;
            
            if(charge<0.){
                ind8 = getallindicesIT(vvec, _begindimx_em, _sizedimx_em, _begindimy_em, _sizedimy_em, _begindimz_em, _sizedimz_em, position);
                mobility=interpolIT(vvec, _beginmap_em, ind8);
            }
            else{
                ind8 = getallindicesIT(vvec, _begindimx_hm, _sizedimx_hm, _begindimy_hm, _sizedimy_hm, _begindimz_hm, _sizedimz_hm, position);
                mobility=interpolIT(vvec, _beginmap_hm, ind8);
            }
            
            
            //get electric field
            ind8 = getallindicesIT(vvec, _begindimx_ef, _sizedimx_ef, _begindimy_ef, _sizedimy_ef, _begindimz_ef, _sizedimz_ef, position);
            auto efvec = interpolITvec(vvec, _beginmap_ef, ind8);
            
            //calculate velocity vector
            auto k1 =charge*mobility*efvec;
            
            //get weighting field
            ind8 = getallindicesIT(vvec, _begindimx_wf, _sizedimx_wf, _begindimy_wf, _sizedimy_wf, _begindimz_wf, _sizedimz_wf, position);
            auto wfvec = interpolITvec(vvec, _beginmap_wf, ind8);
            
            //calculate diffusion velocity
            double s = sqrt(2*(_fDiff*mobility)/_fStep/fabs(hydra::get<_tc_charge>(p)))*hydra::get<_tc_gauss_x>(p); //divide by square root of group size
            
            hydra::Vector3R ref(s*sin(hydra::get<_tc_angle_2>(p))*cos(hydra::get<_tc_angle_1>(p)),s*sin(hydra::get<_tc_angle_2>(p))*sin(hydra::get<_tc_angle_1>(p)),s*cos(hydra::get<_tc_angle_2>(p)));
            
            //set diffusion velocity to be used in evolution step
            hydra::get<_tc_gauss_x>(p) = ref.get(0);
            hydra::get<_tc_gauss_y>(p) = ref.get(1);
            hydra::get<_tc_gauss_z>(p) = ref.get(2);
            
            //calculate total velocity as sum of drift and diffusion
            auto Vdrift = k1 +ref;
            
            //calculate current with Ramo's theorem
            hydra::get<5>(p)=hydra::get<_tc_charge>(p)*HCHARGE*(Vdrift.get(0)*wfvec.get(0)+Vdrift.get(1)*wfvec.get(1)+Vdrift.get(2)*wfvec.get(2));
            
        }
        
        size_t _fN;
        double _fStep;
        double _fDiff;
        Iteratordim _begindimx_ef;
        size_t _sizedimx_ef;
        Iteratordim _begindimy_ef;
        size_t _sizedimy_ef;
        Iteratordim _begindimz_ef;
        size_t _sizedimz_ef;
        Iteratorvec _beginmap_ef;
        size_t _sizemap_ef;
        Iteratordim _begindimx_em;
        size_t _sizedimx_em;
        Iteratordim _begindimy_em;
        size_t _sizedimy_em;
        Iteratordim _begindimz_em;
        size_t _sizedimz_em;
        Iteratorscal _beginmap_em;
        size_t _sizemap_em;
        Iteratordim _begindimx_hm;
        size_t _sizedimx_hm;
        Iteratordim _begindimy_hm;
        size_t _sizedimy_hm;
        Iteratordim _begindimz_hm;
        size_t _sizedimz_hm;
        Iteratorscal _beginmap_hm;
        size_t _sizemap_hm;
        Iteratordim _begindimx_wf;
        size_t _sizedimx_wf;
        Iteratordim _begindimy_wf;
        size_t _sizedimy_wf;
        Iteratordim _begindimz_wf;
        size_t _sizedimz_wf;
        Iteratorvec _beginmap_wf;
        size_t _sizemap_wf;
    };
    
    //functors for time evolution of carrier positions
    template<typename Iteratordim, typename Iteratormap>
    struct Evolution
    {

        Evolution()=delete;                                         // avoid creating objects with undefined data
        
        // No __hydra__dual__ here. Only allow objects to be constructed in host side. 
        Evolution(size_t n, double timestep, Iteratordim begindimx, size_t sizedimx, Iteratordim begindimy, size_t sizedimy, Iteratordim begindimz, size_t sizedimz, Iteratormap beginmap, size_t sizemap):
        _fN(n), _fStep(timestep),
        _begindimx(begindimx), _sizedimx(sizedimx),
        _begindimy(begindimy), _sizedimy(sizedimy),
        _begindimz(begindimz), _sizedimz(sizedimz),
        _beginmap(beginmap), _sizemap(sizemap)
        {}
        
        // Copy constructor needs be __hydra_dual__ to be copyable
        // on device and host sides, in the time of thread distribution.
        // hydra, as well as thrust and std pass arguments by value.
        // Obs.: RooFit passes by reference... (crazy design choice!)
        __hydra_dual__
        Evolution( Evolution const& other):
        _fN(other._fN),  _fStep(other._fStep),
        _begindimx(other._begindimx), _sizedimx(other._sizedimx),
        _begindimy(other._begindimy), _sizedimy(other._sizedimy),
        _begindimz(other._begindimz), _sizedimz(other._sizedimz),
        _beginmap(other._beginmap), _sizemap(other._sizemap)
        {  }
        
        //function call operator needs to be callable in host and device sides 
        // and be constant
        template<typename Particle>
        __hydra_dual__
        void operator()(Particle p){
            //check if particle went outside the volume
            if(hydra::get<_tc_isin>(p)<0) return;
            
            //get particle charge and initial position
            double charge=sgn(hydra::get<_tc_charge>(p));
            double mobility=0;
            hydra::Vector3R position(hydra::get<_tc_x>(p), hydra::get<_tc_y>(p), hydra::get<_tc_z>(p));
            hydra::Vector3R position_init(hydra::get<_tc_x>(p), hydra::get<_tc_y>(p), hydra::get<_tc_z>(p));
            
            hydra::Vector3R vvec(0.,0.,0.);
            mapindices_t ind8;
            mapentry_t entry8;
            
            ind8 = getallindicesIT(vvec, _begindimx, _sizedimx, _begindimy, _sizedimy, _begindimz, _sizedimz, position);
            entry8=interpolITtuple(vvec, _beginmap, ind8);
            
            if(charge<0.){
                mobility = hydra::get<3>(entry8);
            }
            else{
                mobility = hydra::get<4>(entry8);
            }
            //get electric field vector
            hydra::Vector3R efvec(hydra::get<0>(entry8),hydra::get<1>(entry8),hydra::get<2>(entry8));
            
            //calculate velocity vector
            auto k1 =charge*mobility*efvec;
            
            position = position + k1*_fStep*0.5;
            
            //beginning 2nd RK step
            
            ind8 = getallindicesIT(vvec, _begindimx, _sizedimx, _begindimy, _sizedimy, _begindimz, _sizedimz, position);
            entry8=interpolITtuple(vvec, _beginmap, ind8);
            
            if(charge<0.){
                mobility = hydra::get<3>(entry8);
            }
            else{
                mobility = hydra::get<4>(entry8);
            }
            
            //get electric field
            efvec.set(hydra::get<0>(entry8),hydra::get<1>(entry8),hydra::get<2>(entry8));
            
            //calculate velocity vector
            auto k2 =charge*mobility*efvec;
            
            position = position + k2*_fStep*0.5;
            
            //beginning 3rd RK step
            
            ind8 = getallindicesIT(vvec, _begindimx, _sizedimx, _begindimy, _sizedimy, _begindimz, _sizedimz, position);
            entry8=interpolITtuple(vvec, _beginmap, ind8);
            
            if(charge<0.){
                mobility = hydra::get<3>(entry8);
            }
            else{
                mobility = hydra::get<4>(entry8);
            }
            
            //get electric field
            efvec.set(hydra::get<0>(entry8),hydra::get<1>(entry8),hydra::get<2>(entry8));
            
            //calculate velocity vector
            auto k3 =charge*mobility*efvec;
            
            position = position + k3*_fStep;
            
            
            //beginning 4th RK step
            
            ind8 = getallindicesIT(vvec, _begindimx, _sizedimx, _begindimy, _sizedimy, _begindimz, _sizedimz, position);
            entry8=interpolITtuple(vvec, _beginmap, ind8);
            
            if(charge<0.){
                mobility = hydra::get<3>(entry8);
            }
            else{
                mobility = hydra::get<4>(entry8);
            }
            
            //get electric field
            efvec.set(hydra::get<0>(entry8),hydra::get<1>(entry8),hydra::get<2>(entry8));
            
            //calculate velocity vector
            auto k4 =charge*mobility*efvec;
            

            //get diffusion velocity components
            hydra::Vector3R Vdiff(hydra::get<_tc_gauss_x>(p), hydra::get<_tc_gauss_y>(p), hydra::get<_tc_gauss_z>(p));

            //final velocity
            auto Vfinal = 0.1666666667*(k1+2.*k2+2.*k3+k4) + Vdiff;
            
            position = position_init + Vfinal*_fStep;
            
            
            //check if particle is outside defined volume
            if(position.get(0)< *(_begindimx) || position.get(0)> *(_begindimx+_sizedimx-1) || position.get(1)< *(_begindimy) || position.get(1)> *(_begindimy+_sizedimy-1) || position.get(2)< *(_begindimz) || position.get(2)> *(_begindimz+_sizedimz-1)){
                hydra::get<_tc_isin>(p)=-1;
            }
            
            hydra::get<_tc_x>(p) = position.get(0);
            hydra::get<_tc_y>(p) = position.get(1);
            hydra::get<_tc_z>(p) = position.get(2);
        }
        
        size_t _fN;
        double _fStep;
        Iteratordim _begindimx;
        size_t _sizedimx;
        Iteratordim _begindimy;
        size_t _sizedimy;
        Iteratordim _begindimz;
        size_t _sizedimz;
        Iteratormap _beginmap;
        size_t _sizemap;
    };
    
    
    template<typename Iteratordim, typename Iteratorvec, typename Iteratorscal>
    struct Evolution_multi
    {

        Evolution_multi()=delete;                                         // avoid creating objects with undefined data
        
        // No __hydra__dual__ here. Only allow objects to be constructed in host side. 
        Evolution_multi(size_t n, double timestep, Iteratordim begindimx_ef, size_t sizedimx_ef, Iteratordim begindimy_ef, size_t sizedimy_ef, Iteratordim begindimz_ef, size_t sizedimz_ef, 
                          Iteratorvec beginmap_ef, size_t sizemap_ef,
                          Iteratordim begindimx_em, size_t sizedimx_em, Iteratordim begindimy_em, size_t sizedimy_em, Iteratordim begindimz_em, size_t sizedimz_em, 
                          Iteratorscal beginmap_em, size_t sizemap_em,
                          Iteratordim begindimx_hm, size_t sizedimx_hm, Iteratordim begindimy_hm, size_t sizedimy_hm, Iteratordim begindimz_hm, size_t sizedimz_hm, 
                          Iteratorscal beginmap_hm, size_t sizemap_hm,
                          double xmin, double xmax, double ymin, double ymax, double zmin, double zmax):
        _fN(n), _fStep(timestep),
        _xmin(xmin), _xmax(xmax), _ymin(ymin), _ymax(ymax), _zmin(zmin), _zmax(zmax),
        _begindimx_ef(begindimx_ef), _sizedimx_ef(sizedimx_ef),
        _begindimy_ef(begindimy_ef), _sizedimy_ef(sizedimy_ef),
        _begindimz_ef(begindimz_ef), _sizedimz_ef(sizedimz_ef),
        _beginmap_ef(beginmap_ef), _sizemap_ef(sizemap_ef),
        _begindimx_em(begindimx_em), _sizedimx_em(sizedimx_em),
        _begindimy_em(begindimy_em), _sizedimy_em(sizedimy_em),
        _begindimz_em(begindimz_em), _sizedimz_em(sizedimz_em),
        _beginmap_em(beginmap_em), _sizemap_em(sizemap_em),
        _begindimx_hm(begindimx_hm), _sizedimx_hm(sizedimx_hm),
        _begindimy_hm(begindimy_hm), _sizedimy_hm(sizedimy_hm),
        _begindimz_hm(begindimz_hm), _sizedimz_hm(sizedimz_hm),
        _beginmap_hm(beginmap_hm), _sizemap_hm(sizemap_hm)
        {}
        
        // Copy constructor needs be __hydra_dual__ to be copyable
        // on device and host sides, in the time of thread distribution.
        // hydra, as well as thrust and std pass arguments by value.
        // Obs.: RooFit passes by reference... (crazy design choice!)
        __hydra_dual__
        Evolution_multi( Evolution_multi const& other):
        _fN(other._fN),  _fStep(other._fStep),
        _xmin(other._xmin), _xmax(other._xmax), _ymin(other._ymin), _ymax(other._ymax), _zmin(other._zmin), _zmax(other._zmax),
        _begindimx_ef(other._begindimx_ef), _sizedimx_ef(other._sizedimx_ef),
        _begindimy_ef(other._begindimy_ef), _sizedimy_ef(other._sizedimy_ef),
        _begindimz_ef(other._begindimz_ef), _sizedimz_ef(other._sizedimz_ef),
        _beginmap_ef(other._beginmap_ef), _sizemap_ef(other._sizemap_ef),
        _begindimx_em(other._begindimx_em), _sizedimx_em(other._sizedimx_em),
        _begindimy_em(other._begindimy_em), _sizedimy_em(other._sizedimy_em),
        _begindimz_em(other._begindimz_em), _sizedimz_em(other._sizedimz_em),
        _beginmap_em(other._beginmap_em), _sizemap_em(other._sizemap_em),
        _begindimx_hm(other._begindimx_hm), _sizedimx_hm(other._sizedimx_hm),
        _begindimy_hm(other._begindimy_hm), _sizedimy_hm(other._sizedimy_hm),
        _begindimz_hm(other._begindimz_hm), _sizedimz_hm(other._sizedimz_hm),
        _beginmap_hm(other._beginmap_hm), _sizemap_hm(other._sizemap_hm)
        { }
        
        //function call operator needs to be callable in host and device sides 
        // and be constant
        template<typename Particle>
        __hydra_dual__
        void operator()(Particle p){
            
            //check if particle went outside the volume
            if(hydra::get<_tc_isin>(p)<0) return;
            
            //get particle charge and initial position
            double charge=sgn(hydra::get<_tc_charge>(p));
            double mobility=0;
            hydra::Vector3R position(hydra::get<_tc_x>(p), hydra::get<_tc_y>(p), hydra::get<_tc_z>(p));
            hydra::Vector3R position_init(hydra::get<_tc_x>(p), hydra::get<_tc_y>(p), hydra::get<_tc_z>(p));
            
            hydra::Vector3R vvec(0.,0.,0.);
            mapindices_t ind8;
            
            //beginning 1st RK step
            
            if(charge<0.){
                ind8 = getallindicesIT(vvec, _begindimx_em, _sizedimx_em, _begindimy_em, _sizedimy_em, _begindimz_em, _sizedimz_em, position);
                mobility=interpolIT(vvec, _beginmap_em, ind8);
            }
            else{
                ind8 = getallindicesIT(vvec, _begindimx_hm, _sizedimx_hm, _begindimy_hm, _sizedimy_hm, _begindimz_hm, _sizedimz_hm, position);
                mobility=interpolIT(vvec, _beginmap_hm, ind8);
            }
            
            //get electric field
            ind8 = getallindicesIT(vvec, _begindimx_ef, _sizedimx_ef, _begindimy_ef, _sizedimy_ef, _begindimz_ef, _sizedimz_ef, position);
            auto efvec = interpolITvec(vvec, _beginmap_ef, ind8);
            
            //calculate velocity vector
            auto k1 =charge*mobility*efvec;
            
            position = position + k1*_fStep*0.5;
            
            //beginning 2nd RK step
             
            if(charge<0.){
                ind8 = getallindicesIT(vvec, _begindimx_em, _sizedimx_em, _begindimy_em, _sizedimy_em, _begindimz_em, _sizedimz_em, position);
                mobility=interpolIT(vvec, _beginmap_em, ind8);
            }
            else{
                ind8 = getallindicesIT(vvec, _begindimx_hm, _sizedimx_hm, _begindimy_hm, _sizedimy_hm, _begindimz_hm, _sizedimz_hm, position);
                mobility=interpolIT(vvec, _beginmap_hm, ind8);
            }
            
            //get electric field
            ind8 = getallindicesIT(vvec, _begindimx_ef, _sizedimx_ef, _begindimy_ef, _sizedimy_ef, _begindimz_ef, _sizedimz_ef, position);
            efvec = interpolITvec(vvec, _beginmap_ef, ind8);
            
            //calculate velocity vector
            auto k2 =charge*mobility*efvec;
            
            position = position + k2*_fStep*0.5;
            
            //beginning 3rd RK step
            
            if(charge<0.){
                ind8 = getallindicesIT(vvec, _begindimx_em, _sizedimx_em, _begindimy_em, _sizedimy_em, _begindimz_em, _sizedimz_em, position);
                mobility=interpolIT(vvec, _beginmap_em, ind8);
            }
            else{
                ind8 = getallindicesIT(vvec, _begindimx_hm, _sizedimx_hm, _begindimy_hm, _sizedimy_hm, _begindimz_hm, _sizedimz_hm, position);
                mobility=interpolIT(vvec, _beginmap_hm, ind8);
            }
            
            //get electric field
            ind8 = getallindicesIT(vvec, _begindimx_ef, _sizedimx_ef, _begindimy_ef, _sizedimy_ef, _begindimz_ef, _sizedimz_ef, position);
            efvec = interpolITvec(vvec, _beginmap_ef, ind8);
            
            //calculate velocity vector
            auto k3 =charge*mobility*efvec;
            
            position = position + k3*_fStep;
            
            
            //beginning 4th RK step
            
            if(charge<0.){
                ind8 = getallindicesIT(vvec, _begindimx_em, _sizedimx_em, _begindimy_em, _sizedimy_em, _begindimz_em, _sizedimz_em, position);
                mobility=interpolIT(vvec, _beginmap_em, ind8);
            }
            else{
                ind8 = getallindicesIT(vvec, _begindimx_hm, _sizedimx_hm, _begindimy_hm, _sizedimy_hm, _begindimz_hm, _sizedimz_hm, position);
                mobility=interpolIT(vvec, _beginmap_hm, ind8);
            }
            
            //get electric field
            ind8 = getallindicesIT(vvec, _begindimx_ef, _sizedimx_ef, _begindimy_ef, _sizedimy_ef, _begindimz_ef, _sizedimz_ef, position);
            efvec = interpolITvec(vvec, _beginmap_ef, ind8);
            
            //calculate velocity vector
            auto k4 =charge*mobility*efvec;
            

            //get diffusion velocity components
            hydra::Vector3R Vdiff(hydra::get<_tc_gauss_x>(p), hydra::get<_tc_gauss_y>(p), hydra::get<_tc_gauss_z>(p));

            //final velocity
            auto Vfinal = 0.1666666667*(k1+2.*k2+2.*k3+k4) + Vdiff;
            
            position = position_init + Vfinal*_fStep;
            
            
            //check if particle is outside defined volume
            if(position.get(0)<_xmin || position.get(0)>_xmax || position.get(1)<_ymin || position.get(1)>_ymax || position.get(2)<_zmin || position.get(2)>_zmax){
                hydra::get<_tc_isin>(p)=-1;
            }
            
            hydra::get<_tc_x>(p) = position.get(0);
            hydra::get<_tc_y>(p) = position.get(1);
            hydra::get<_tc_z>(p) = position.get(2);
            
        }
        
        
        
        size_t _fN;
        double _fStep;
        Iteratordim _begindimx_ef;
        size_t _sizedimx_ef;
        Iteratordim _begindimy_ef;
        size_t _sizedimy_ef;
        Iteratordim _begindimz_ef;
        size_t _sizedimz_ef;
        Iteratorvec _beginmap_ef;
        size_t _sizemap_ef;
        Iteratordim _begindimx_em;
        size_t _sizedimx_em;
        Iteratordim _begindimy_em;
        size_t _sizedimy_em;
        Iteratordim _begindimz_em;
        size_t _sizedimz_em;
        Iteratorscal _beginmap_em;
        size_t _sizemap_em;
        Iteratordim _begindimx_hm;
        size_t _sizedimx_hm;
        Iteratordim _begindimy_hm;
        size_t _sizedimy_hm;
        Iteratordim _begindimz_hm;
        size_t _sizedimz_hm;
        Iteratorscal _beginmap_hm;
        size_t _sizemap_hm;
        double _xmin;
        double _xmax;
        double _ymin;
        double _ymax;
        double _zmin;
        double _zmax;
    };
    
    
    //helper functions to initialize functors
    template<typename Iterabledim, typename Iterablemap>
    RamoCurrent< typename Iterabledim::iterator, typename Iterablemap::iterator> 
    make_RamoCurrent(size_t n, double timestep, double diffusion,  Iterabledim& datadimx,  Iterabledim& datadimy,  Iterabledim& datadimz, Iterablemap& datamap ) {
    
        return RamoCurrent< typename Iterabledim::iterator, typename Iterablemap::iterator>( n, timestep, diffusion, 
                                                                                std::forward<Iterabledim>(datadimx).begin(), std::forward<Iterabledim>(datadimx).size(), 
                                                                                std::forward<Iterabledim>(datadimy).begin(), std::forward<Iterabledim>(datadimy).size(),    
                                                                                std::forward<Iterabledim>(datadimz).begin(), std::forward<Iterabledim>(datadimz).size(), 
                                                                                std::forward<Iterablemap>(datamap).begin(),  std::forward<Iterablemap>(datamap).size());
    }
    
    template<typename Iterabledim, typename Iterablevec, typename Iterablescal>
    RamoCurrent_multi< typename Iterabledim::iterator, typename Iterablevec::iterator, typename Iterablescal::iterator>
    make_RamoCurrent(size_t n, double timestep, double diffusion,  
                          Iterabledim& datadimx_ef,  Iterabledim& datadimy_ef,  Iterabledim& datadimz_ef, Iterablevec& datamap_ef,
                          Iterabledim& datadimx_em,  Iterabledim& datadimy_em,  Iterabledim& datadimz_em, Iterablescal& datamap_em,
                          Iterabledim& datadimx_hm,  Iterabledim& datadimy_hm,  Iterabledim& datadimz_hm, Iterablescal& datamap_hm,
                          Iterabledim& datadimx_wf,  Iterabledim& datadimy_wf,  Iterabledim& datadimz_wf, Iterablevec& datamap_wf
                         ) {
    
        return RamoCurrent_multi< typename Iterabledim::iterator, typename Iterablevec::iterator, typename Iterablescal::iterator>( n, timestep, diffusion, 
                                                                                std::forward<Iterabledim>(datadimx_ef).begin(), std::forward<Iterabledim>(datadimx_ef).size(), 
                                                                                std::forward<Iterabledim>(datadimy_ef).begin(), std::forward<Iterabledim>(datadimy_ef).size(),    
                                                                                std::forward<Iterabledim>(datadimz_ef).begin(), std::forward<Iterabledim>(datadimz_ef).size(), 
                                                                                std::forward<Iterablevec>(datamap_ef).begin(),  std::forward<Iterablevec>(datamap_ef).size(),
                                                                                std::forward<Iterabledim>(datadimx_em).begin(), std::forward<Iterabledim>(datadimx_em).size(), 
                                                                                std::forward<Iterabledim>(datadimy_em).begin(), std::forward<Iterabledim>(datadimy_em).size(),    
                                                                                std::forward<Iterabledim>(datadimz_em).begin(), std::forward<Iterabledim>(datadimz_em).size(), 
                                                                                std::forward<Iterablescal>(datamap_em).begin(), std::forward<Iterablescal>(datamap_em).size(),
                                                                                std::forward<Iterabledim>(datadimx_hm).begin(), std::forward<Iterabledim>(datadimx_hm).size(), 
                                                                                std::forward<Iterabledim>(datadimy_hm).begin(), std::forward<Iterabledim>(datadimy_hm).size(),    
                                                                                std::forward<Iterabledim>(datadimz_hm).begin(), std::forward<Iterabledim>(datadimz_hm).size(), 
                                                                                std::forward<Iterablescal>(datamap_hm).begin(), std::forward<Iterablescal>(datamap_hm).size(),
                                                                                std::forward<Iterabledim>(datadimx_wf).begin(), std::forward<Iterabledim>(datadimx_wf).size(), 
                                                                                std::forward<Iterabledim>(datadimy_wf).begin(), std::forward<Iterabledim>(datadimy_wf).size(),    
                                                                                std::forward<Iterabledim>(datadimz_wf).begin(), std::forward<Iterabledim>(datadimz_wf).size(), 
                                                                                std::forward<Iterablevec>(datamap_wf).begin(),  std::forward<Iterablevec>(datamap_wf).size());
    }
    
    template<typename Iterabledim, typename Iterablemap>
    Evolution< typename Iterabledim::iterator, typename Iterablemap::iterator> 
    make_Evolution(size_t n, double timestep,  Iterabledim& datadimx,  Iterabledim& datadimy,  Iterabledim& datadimz, Iterablemap& datamap ) {
    
        return Evolution< typename Iterabledim::iterator, typename Iterablemap::iterator>( n, timestep,
                                                                                std::forward<Iterabledim>(datadimx).begin(), std::forward<Iterabledim>(datadimx).size(), 
                                                                                std::forward<Iterabledim>(datadimy).begin(), std::forward<Iterabledim>(datadimy).size(),    
                                                                                std::forward<Iterabledim>(datadimz).begin(), std::forward<Iterabledim>(datadimz).size(), 
                                                                                std::forward<Iterablemap>(datamap).begin(),  std::forward<Iterablemap>(datamap).size());
    }
    
    template<typename Iterabledim, typename Iterablevec, typename Iterablescal>
    Evolution_multi< typename Iterabledim::iterator, typename Iterablevec::iterator, typename Iterablescal::iterator>
    make_Evolution(size_t n, double timestep, 
                          Iterabledim& datadimx_ef,  Iterabledim& datadimy_ef,  Iterabledim& datadimz_ef, Iterablevec& datamap_ef,
                          Iterabledim& datadimx_em,  Iterabledim& datadimy_em,  Iterabledim& datadimz_em, Iterablescal& datamap_em,
                          Iterabledim& datadimx_hm,  Iterabledim& datadimy_hm,  Iterabledim& datadimz_hm, Iterablescal& datamap_hm,
                          double xmin, double xmax, double ymin, double ymax, double zmin, double zmax
                         ) {
    
        return Evolution_multi< typename Iterabledim::iterator, typename Iterablevec::iterator, typename Iterablescal::iterator>( n, timestep, 
                                                                                std::forward<Iterabledim>(datadimx_ef).begin(), std::forward<Iterabledim>(datadimx_ef).size(), 
                                                                                std::forward<Iterabledim>(datadimy_ef).begin(), std::forward<Iterabledim>(datadimy_ef).size(),    
                                                                                std::forward<Iterabledim>(datadimz_ef).begin(), std::forward<Iterabledim>(datadimz_ef).size(), 
                                                                                std::forward<Iterablevec>(datamap_ef).begin(),  std::forward<Iterablevec>(datamap_ef).size(),
                                                                                std::forward<Iterabledim>(datadimx_em).begin(), std::forward<Iterabledim>(datadimx_em).size(), 
                                                                                std::forward<Iterabledim>(datadimy_em).begin(), std::forward<Iterabledim>(datadimy_em).size(),    
                                                                                std::forward<Iterabledim>(datadimz_em).begin(), std::forward<Iterabledim>(datadimz_em).size(), 
                                                                                std::forward<Iterablescal>(datamap_em).begin(), std::forward<Iterablescal>(datamap_em).size(),
                                                                                std::forward<Iterabledim>(datadimx_hm).begin(), std::forward<Iterabledim>(datadimx_hm).size(), 
                                                                                std::forward<Iterabledim>(datadimy_hm).begin(), std::forward<Iterabledim>(datadimy_hm).size(),    
                                                                                std::forward<Iterabledim>(datadimz_hm).begin(), std::forward<Iterabledim>(datadimz_hm).size(), 
                                                                                std::forward<Iterablescal>(datamap_hm).begin(), std::forward<Iterablescal>(datamap_hm).size(),
                                                                                xmin, xmax, ymin, ymax, zmin, zmax);
    }
    
}
    
#endif
