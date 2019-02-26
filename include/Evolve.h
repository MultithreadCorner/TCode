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
    
    __hydra_dual__
    inline double interpol(double vvec[], double v[], size_t ind8[8]){     
        return val(vvec[0], vvec[1], vvec[2], v[ind8[0]], v[ind8[1]], v[ind8[2]], v[ind8[3]], v[ind8[4]], v[ind8[5]], v[ind8[6]], v[ind8[7]]);
        
    }
    
    /*
    review A.A. 
    why not using std::binary_search (https://en.cppreference.com/w/cpp/algorithm/binary_search) ?
    if V is too large (O(6) or higher) and it becomes too slow, we can wrapp a vectorized version 
    in hydra...  
    */
    
    //function to get index in one dimension (binary search)
    __hydra_dual__
    inline size_t getindex(double v[], size_t n, double f, size_t first=0){
        if(f<=v[first]){
            return first+1;

        }
        if(f>=v[n-1]){
            return n-1;
        }
        size_t lo = first;
        size_t hi = n - 1;
        size_t mid=0;
        while(lo<hi){

             mid=(hi+lo)/2;
             
             if(v[mid]==f) return mid+1;

             if(f<v[mid]){
                 if(mid>0 && f>v[mid-1]) return mid;
                 hi=mid;
             }
             else{
                 if(mid<n-1 && f < v[mid+1]) return mid+1;
                 lo=mid+1;
             }
         }
         return 0;
    }
    
    //function to get dim indices
    __hydra_dual__
    inline void getindices(size_t ind[3], double vx[], size_t sx, double vy[], size_t sy, double vz[], size_t sz, double x, double y, double z){
        ind[0]=getindex(vx, sx, x);
        ind[1]=getindex(vy, sy, y);
        ind[2]=getindex(vz, sz, z);
    }
    
    //function to get vectors and indices in in phys tables
    __hydra_dual__
    inline void getallindices(double vvec[3], size_t ind8[8], double vx[], size_t sx, double vy[], size_t sy, double vz[], size_t sz, double x, double y, double z){
        size_t ind[3];
        getindices(ind, vx, sx, vy, sy, vz, sz, x, y, z);
        vvec[0] = (x-vx[ind[0]-1])/(vx[ind[0]]-vx[ind[0]-1]);
        vvec[1] = (y-vy[ind[1]-1])/(vy[ind[1]]-vy[ind[1]-1]);
        vvec[2] = (z-vz[ind[2]-1])/(vz[ind[2]]-vz[ind[2]-1]);
        
        double fNPXY=sx*sy;
        
        ind8[0] = ind[1]-1 + (ind[0]-1)*sy+(ind[2]-1)*fNPXY;
        ind8[1] = ind[1]-1 + (ind[0])*sy+(ind[2]-1)*fNPXY;
        ind8[2] = ind[1] + (ind[0]-1)*sy+(ind[2]-1)*fNPXY;
        ind8[3] = ind[1]-1 + (ind[0]-1)*sy+(ind[2])*fNPXY;
        ind8[4] = ind[1] + (ind[0])*sy+(ind[2]-1)*fNPXY;
        ind8[5] = ind[1] + (ind[0]-1)*sy+(ind[2])*fNPXY;
        ind8[6] = ind[1]-1 +(ind[0])*sy+(ind[2])*fNPXY;
        ind8[7] = ind[1] + (ind[0])*sy+(ind[2])*fNPXY;
    }
    /*
    Review A.A.
    
    ApplyRamo  is  receiving and using pointers without testing the validity. 
    
    */
    //calculate currents
    struct ApplyRamo
    {

        ApplyRamo()=delete;                                         // avoid creating objects with undefined data

         /*
         Review A.A.
    
         ApplyRamo ctor is  receiving and using pointers without testing the validity. 
         */
        // No __hydra__dual__ here. Only allow objects to be constructed in host side. 
        ApplyRamo(size_t n, double timestep, double diffusion, double *in_x, size_t dim_x,double *in_y, size_t dim_y,double *in_z, size_t dim_z, double *in_xef,double *in_yef, double *in_zef, double *in_em, double *in_hm, double *in_xwf,double *in_ywf, double *in_zwf):
        fN(n),
        fStep(timestep),
        fDiff(diffusion),
        fNPX(dim_x),
        fNPY(dim_y),
        fNPZ(dim_z),
        xvec(in_x),
        yvec(in_y),
        zvec(in_z),
        xvecef(in_xef),
        yvecef(in_yef),
        zvecef(in_zef),
        vecem(in_em),
        vechm(in_hm),
        xvecwf(in_xwf),
        yvecwf(in_ywf),
        zvecwf(in_zwf)
        {}
        
        //
        // No __hydra__dual__ here. Only allow objects to be constructed in host side. 
        ApplyRamo(size_t n, double timestep, double diffusion, VecDev_t<double> &in_x, VecDev_t<double> &in_y, VecDev_t<double> &in_z, VecDev_t<double> &in_xef, VecDev_t<double> &in_yef, VecDev_t<double> &in_zef, VecDev_t<double> &in_em, VecDev_t<double> &in_hm, VecDev_t<double> &in_xwf, VecDev_t<double> &in_ywf, VecDev_t<double> &in_zwf):
        fN(n),
        fStep(timestep),
        fDiff(diffusion),
        fNPX(in_x.size()),
        fNPY(in_y.size()),
        fNPZ(in_z.size()),
        xvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_x.data())),
        yvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_y.data())),
        zvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_z.data())),
        xvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_xef.data())),
        yvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_yef.data())),
        zvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_zef.data())),
        vecem(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_em.data())),
        vechm(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_hm.data())),
        xvecwf(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_xwf.data())),
        yvecwf(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_ywf.data())),
        zvecwf(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_zwf.data()))
        {}

        // Copy constructor needs be __hydra_dual__ to be copyable
        // on device and host sides, in the time of thread distribution.
        // hydra, as well as thrust and std pass arguments by value.
        // Obs.: RooFit passes by reference... (crazy design choice!)
        __hydra_dual__
        ApplyRamo( ApplyRamo const& other):
        fN(other.fN),
        fStep(other.fStep),
        fDiff(other.fDiff),
        fNPX(other.fNPX),
        fNPY(other.fNPY),
        fNPZ(other.fNPZ),
        xvec(other.xvec),
        yvec(other.yvec),
        zvec(other.zvec),
        xvecef(other.xvecef),
        yvecef(other.yvecef),
        zvecef(other.zvecef),
        vecem(other.vecem),
        vechm(other.vechm),
        xvecwf(other.xvecwf),
        yvecwf(other.yvecwf),
        zvecwf(other.zvecwf)
        { }

        __hydra_dual__
        ApplyRamo operator=( ApplyRamo const& other){
            if(this == &other) return *this;

            fN=other.fN;
            fStep=other.fStep;
            fDiff=other.fDiff;
            fNPX=other.fNPX;
            fNPY=other.fNPY;
            fNPZ=other.fNPZ;
            xvec=other.xvec;
            yvec=other.yvec;
            zvec=other.zvec;
            xvecef=other.xvecef;
            yvecef=other.yvecef;
            zvecef=other.zvecef;
            vecem=other.vecem;
            vechm=other.vechm;
            xvecwf=other.xvecwf;
            yvecwf=other.yvecwf;
            zvecwf=other.zvecwf;

            return *this;
        }

        // FIXME: is not constant !
        //function call operator needs to be callable in host and device sides 
        // and be constant
        template<typename Particle>
        __hydra_dual__
        void operator()(Particle p){
            
            //check if particle went outside the volume
            if(hydra::get<_tc_isin>(p)<0){
                hydra::get<_tc_curr>(p)=0.;
                return;
            }
            
            //get particle charge and initial position
            double charge=hydra::get<0>(p);
            double mobility=0;
            double x=hydra::get<_tc_x>(p);
            double y=hydra::get<_tc_y>(p);
            double z=hydra::get<_tc_z>(p);
            
            double vvec[3];
            size_t ind8[8];
            
            getallindices(vvec, ind8, xvec, fNPX, yvec, fNPY, zvec, fNPZ, x, y, z);
            if(charge<0.) mobility=interpol(vvec, vecem, ind8);
            else mobility=interpol(vvec, vechm, ind8);
            
            //get electric field
            double ex = interpol(vvec, xvecef, ind8);
            double ey = interpol(vvec, yvecef, ind8);
            double ez = interpol(vvec, zvecef, ind8);
            
            //get weighting field
            double wfx = interpol(vvec, xvecwf, ind8);
            double wfy = interpol(vvec, yvecwf, ind8);
            double wfz = interpol(vvec, zvecwf, ind8);

            //first step only
            double k1x=charge*mobility*ex;
            double k1y=charge*mobility*ey;
            double k1z=charge*mobility*ez;
            
            //calculate diffusion velocity
            double s = ::sqrt(2*(fDiff*mobility)/fStep)*hydra::get<_tc_gauss_x>(p);
            
            hydra::Vector3R direction(ex,ey,ez);
            hydra::Vector3R ref(1.,0.,0.);
            double theta = hydra::get<_tc_angle_2>(p);
            double phi = hydra::get<_tc_angle_1>(p);
            if(direction.d3mag()==0){
                ref.set(s*sin(theta)*cos(phi),s*sin(theta)*sin(phi),s*cos(theta));
            }
            else{
            
                direction = direction.unit();
                ref = (direction.cross(ref)).unit()*s;
                //rotate ref
                ref = ref*cos(phi) + (direction.cross(ref))*sin(phi)+direction*(direction.dot(ref))*(1-cos(phi));
            }
            
            hydra::get<_tc_gauss_x>(p) = ref.get(0);
            hydra::get<_tc_gauss_y>(p) = ref.get(1);
            hydra::get<_tc_gauss_z>(p) = ref.get(2);
            
            double Vdriftx=k1x+ref.get(0);
            double Vdrifty=k1y+ref.get(1);
            double Vdriftz=k1z+ref.get(2);
            
//             Vdrift=sqrt((Vdriftx*Vdriftx+Vdrifty*Vdrifty+Vdriftz*Vdriftz));
            
//             MaxV=MAXVDRIFTE;
//             if(charge>0) MaxV=MAXVDRIFTH;
//             
//             if(Vdrift>MaxV){
//                 Vdriftx*=(MaxV/(Vdrift));
//                 Vdrifty*=(MaxV/(Vdrift));
//                 Vdriftz*=(MaxV/(Vdrift));
//             }
            
            //calculate current at this instant
            hydra::get<5>(p)=charge*HCHARGE*((Vdriftx)*wfx+(Vdrifty)*wfy+(Vdriftz)*wfz);
            

        }
        
        //FIXME: pointers initialized here and in the class ctors. Consider to enable the default ctor.
        //anyway you never test the pointers for null.
        size_t fN;
        double fStep;
        double fDiff;
        size_t fNPX = 0;
        size_t fNPY = 0;
        size_t fNPZ = 0;
        double *xvecwf = nullptr;
        double *yvecwf = nullptr;
        double *zvecwf = nullptr;
        double *xvecef = nullptr;
        double *yvecef = nullptr;
        double *zvecef = nullptr;
        double *vecem  = nullptr;
        double *vechm  = nullptr;
        double *xvec   = nullptr;
        double *yvec   = nullptr;
        double *zvec   = nullptr;
    };
    
    //fixme: why not using the iterators instead of raw pointer?
    // your only ctor is taking vectors ...
    //calculate currents
    struct ApplyRamo_multi
    {

        //fixme: Actually you are anyway initliazing the struct with null-pointers 
        ApplyRamo_multi()=delete;                                         // avoid creating objects with undefined data
        
        // No __hydra__dual__ here. Only allow objects to be constructed in host side. 
        ApplyRamo_multi(size_t n, double timestep, double diffusion, VecDev_t<double> &inef_x, VecDev_t<double> &inef_y, VecDev_t<double> &inef_z, VecDev_t<double> &in_xef, VecDev_t<double> &in_yef, VecDev_t<double> &in_zef, VecDev_t<double> &inem_x, VecDev_t<double> &inem_y, VecDev_t<double> &inem_z, VecDev_t<double> &in_em, VecDev_t<double> &inhm_x, VecDev_t<double> &inhm_y, VecDev_t<double> &inhm_z, VecDev_t<double> &in_hm,  VecDev_t<double> &inwf_x, VecDev_t<double> &inwf_y, VecDev_t<double> &inwf_z, VecDev_t<double> &in_xwf, VecDev_t<double> &in_ywf, VecDev_t<double> &in_zwf):
        fN(n),
        fStep(timestep),
        fDiff(diffusion),
        fNPXef(inef_x.size()),
        fNPYef(inef_y.size()),
        fNPZef(inef_z.size()),
        xefvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inef_x.data())),
        yefvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inef_y.data())),
        zefvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inef_z.data())),
        fNPXwf(inwf_x.size()),
        fNPYwf(inwf_y.size()),
        fNPZwf(inwf_z.size()),
        xwfvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inwf_x.data())),
        ywfvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inwf_y.data())),
        zwfvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inwf_z.data())),
        fNPXem(inem_x.size()),
        fNPYem(inem_y.size()),
        fNPZem(inem_z.size()),
        xemvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inem_x.data())),
        yemvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inem_y.data())),
        zemvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inem_z.data())),
        fNPXhm(inhm_x.size()),
        fNPYhm(inhm_y.size()),
        fNPZhm(inhm_z.size()),
        xhmvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inhm_x.data())),
        yhmvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inhm_y.data())),
        zhmvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inhm_z.data())),
        xvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_xef.data())),
        yvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_yef.data())),
        zvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_zef.data())),
        vecem(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_em.data())),
        vechm(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_hm.data())),
        xvecwf(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_xwf.data())),
        yvecwf(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_ywf.data())),
        zvecwf(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_zwf.data()))
        {}

        // Copy constructor needs be __hydra_dual__ to be copyable
        // on device and host sides, in the time of thread distribution.
        // hydra, as well as thrust and std pass arguments by value.
        // Obs.: RooFit passes by reference... (crazy design choice!)
        __hydra_dual__
        ApplyRamo_multi( ApplyRamo_multi const& other):
        fN(other.fN),
        fStep(other.fStep),
        fDiff(other.fDiff),
        fNPXef(other.fNPXef),
        fNPYef(other.fNPYef),
        fNPZef(other.fNPZef),
        xefvec(other.xefvec),
        yefvec(other.yefvec),
        zefvec(other.zefvec),
        fNPXwf(other.fNPXwf),
        fNPYwf(other.fNPYwf),
        fNPZwf(other.fNPZwf),
        xwfvec(other.xwfvec),
        ywfvec(other.ywfvec),
        zwfvec(other.zwfvec),
        fNPXem(other.fNPXem),
        fNPYem(other.fNPYem),
        fNPZem(other.fNPZem),
        xemvec(other.xemvec),
        yemvec(other.yemvec),
        zemvec(other.zemvec),
        fNPXhm(other.fNPXhm),
        fNPYhm(other.fNPYhm),
        fNPZhm(other.fNPZhm),
        xhmvec(other.xhmvec),
        yhmvec(other.yhmvec),
        zhmvec(other.zhmvec),
        xvecef(other.xvecef),
        yvecef(other.yvecef),
        zvecef(other.zvecef),
        vecem(other.vecem),
        vechm(other.vechm),
        xvecwf(other.xvecwf),
        yvecwf(other.yvecwf),
        zvecwf(other.zvecwf)
        { }

        __hydra_dual__
        ApplyRamo_multi operator=( ApplyRamo_multi const& other){
            if(this == &other) return *this;

            fN=other.fN;
            fStep=other.fStep;
            fDiff=other.fDiff;
            fNPXef=other.fNPXef;
            fNPYef=other.fNPYef;
            fNPZef=other.fNPZef;
            xefvec=other.xefvec;
            yefvec=other.yefvec;
            zefvec=other.zefvec;
            fNPXwf=other.fNPXwf;
            fNPYwf=other.fNPYwf;
            fNPZwf=other.fNPZwf;
            xwfvec=other.xwfvec;
            ywfvec=other.ywfvec;
            zwfvec=other.zwfvec;
            fNPXem=other.fNPXem;
            fNPYem=other.fNPYem;
            fNPZem=other.fNPZem;
            xemvec=other.xemvec;
            yemvec=other.yemvec;
            zemvec=other.zemvec;
            fNPXef=other.fNPXef;
            fNPYef=other.fNPYef;
            fNPZef=other.fNPZef;
            xefvec=other.xefvec;
            yefvec=other.yefvec;
            zefvec=other.zefvec;
            xvecef=other.xvecef;
            yvecef=other.yvecef;
            zvecef=other.zvecef;
            vecem=other.vecem;
            vechm=other.vechm;
            xvecwf=other.xvecwf;
            yvecwf=other.yvecwf;
            zvecwf=other.zvecwf;

            return *this;
        }

        //function call operator needs to be callable in host and device sides 
        // and be constant
        template<typename Particle>
        __hydra_dual__
        void operator()(Particle p){
            
            //check if particle went outside the volume
            if(hydra::get<_tc_isin>(p)<0){
                hydra::get<_tc_curr>(p)=0.;
                return;
            }
            
            //get particle charge and initial position
            double charge=hydra::get<_tc_charge>(p);
            double mobility=0;
            double x=hydra::get<_tc_x>(p);
            double y=hydra::get<_tc_y>(p);
            double z=hydra::get<_tc_z>(p);
            
            double vvec[3];
            size_t ind8[8];
            
            if(charge<0.){
                getallindices(vvec, ind8, xemvec, fNPXem, yemvec, fNPYem, zemvec, fNPZem, x, y, z);
                mobility=interpol(vvec, vecem, ind8);
            }
            else{
                getallindices(vvec, ind8, xhmvec, fNPXhm, yhmvec, fNPYhm, zhmvec, fNPZhm, x, y, z);
                mobility=interpol(vvec, vechm, ind8);
            }
            
            //get electric field
            getallindices(vvec, ind8, xefvec, fNPXef, yefvec, fNPYef, zefvec, fNPZef, x, y, z);
            double ex = interpol(vvec, xvecef, ind8);
            double ey = interpol(vvec, yvecef, ind8);
            double ez = interpol(vvec, zvecef, ind8);
            
            //get weighting field
            getallindices(vvec, ind8, xwfvec, fNPXwf, ywfvec, fNPYwf, zwfvec, fNPZwf, x, y, z);
            double wfx = interpol(vvec, xvecwf, ind8);
            double wfy = interpol(vvec, yvecwf, ind8);
            double wfz = interpol(vvec, zvecwf, ind8);

            //first step only
            double k1x=charge*mobility*ex;
            double k1y=charge*mobility*ey;
            double k1z=charge*mobility*ez;
            //calculate diffusion velocity
            double s = sqrt(2*(fDiff*mobility)/fStep)*hydra::get<_tc_gauss_x>(p);
            hydra::Vector3R direction(ex,ey,ez);
            hydra::Vector3R ref(1.,0.,0.);
            
            double theta = hydra::get<_tc_angle_2>(p);
            double phi = hydra::get<_tc_angle_1>(p);
            if(direction.d3mag()==0){
                ref.set(s*sin(theta)*cos(phi),s*sin(theta)*sin(phi),s*cos(theta));
            }
            else{
            
                direction = direction.unit();
                ref = (direction.cross(ref)).unit()*s;
                //rotate ref
                ref = ref*cos(phi) + (direction.cross(ref))*sin(phi)+direction*(direction.dot(ref))*(1-cos(phi));
            }
            
            hydra::get<_tc_gauss_x>(p) = ref.get(0);
            hydra::get<_tc_gauss_y>(p) = ref.get(1);
            hydra::get<_tc_gauss_z>(p) = ref.get(2);
            
            double Vdriftx=k1x+ref.get(0);
            double Vdrifty=k1y+ref.get(1);
            double Vdriftz=k1z+ref.get(2);
            
//             Vdrift=sqrt((Vdriftx*Vdriftx+Vdrifty*Vdrifty+Vdriftz*Vdriftz));
            
//             MaxV=MAXVDRIFTE;
//             if(charge>0) MaxV=MAXVDRIFTH;
//             
//             if(Vdrift>MaxV){
//                 Vdriftx*=(MaxV/(Vdrift));
//                 Vdrifty*=(MaxV/(Vdrift));
//                 Vdriftz*=(MaxV/(Vdrift));
//             }
            
            //calculate current at this instant
            hydra::get<5>(p)=charge*HCHARGE*((Vdriftx)*wfx+(Vdrifty)*wfy+(Vdriftz)*wfz);

        }
        
        
        size_t fN;
        double fStep;
        double fDiff;
        size_t fNPXef = 0;
        size_t fNPYef = 0;
        size_t fNPZef = 0;
        size_t fNPXwf = 0;
        size_t fNPYwf = 0;
        size_t fNPZwf = 0;
        size_t fNPXem = 0;
        size_t fNPYem = 0;
        size_t fNPZem = 0;
        size_t fNPXhm = 0;
        size_t fNPYhm = 0;
        size_t fNPZhm = 0;
        double *xvecwf = nullptr;
        double *yvecwf = nullptr;
        double *zvecwf = nullptr;
        double *xvecef = nullptr;
        double *yvecef = nullptr;
        double *zvecef = nullptr;
        double *vecem  = nullptr;
        double *vechm  = nullptr;
        double *xefvec   = nullptr;
        double *yefvec   = nullptr;
        double *zefvec   = nullptr;
        double *xwfvec   = nullptr;
        double *ywfvec   = nullptr;
        double *zwfvec   = nullptr;
        double *xemvec   = nullptr;
        double *yemvec   = nullptr;
        double *zemvec   = nullptr;
        double *xhmvec   = nullptr;
        double *yhmvec   = nullptr;
        double *zhmvec   = nullptr;
    };
    
    //fixme: test pointer validity before to use!
    struct Evolve
    {
       //fixme: Actually you are anyway initliazing the struct with null-pointers 
        Evolve()=delete;                                         // avoid creating objects with undefined data

        // No __hydra__dual__ here. Only allow objects to be constructed in host side. 
        Evolve(size_t n, double timestep, double *in_x, size_t dim_x,double *in_y, size_t dim_y,double *in_z, size_t dim_z, double *in_xef,double *in_yef, double *in_zef, double *in_em):
        fN(n),
        fStep(timestep),
        fNPX(dim_x),
        fNPY(dim_y),
        fNPZ(dim_z),
        xvec(in_x),
        yvec(in_y),
        zvec(in_z),
        xvecef(in_xef),
        yvecef(in_yef),
        zvecef(in_zef),
        vecem(in_em),
        vechm(in_hm)
        {}
        
        // No __hydra__dual__ here. Only allow objects to be constructed in host side. 
        Evolve(size_t n, double timestep, VecDev_t<double> &in_x, VecDev_t<double> &in_y, VecDev_t<double> &in_z, VecDev_t<double> &in_xef, VecDev_t<double> &in_yef, VecDev_t<double> &in_zef, VecDev_t<double> &in_em, VecDev_t<double> &in_hm):
        fN(n),
        fStep(timestep),
        fNPX(in_x.size()),
        fNPY(in_y.size()),
        fNPZ(in_z.size()),
        xvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_x.data())),
        yvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_y.data())),
        zvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_z.data())),
        xvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_xef.data())),
        yvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_yef.data())),
        zvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_zef.data())),
        vecem(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_em.data())),
        vechm(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_hm.data()))
        {}
        
        // Copy constructor needs be __hydra_dual__ to be copyable
        // on device and host sides, in the time of thread distribution.
        // hydra, as well as thrust and std pass arguments by value.
        // Obs.: RooFit passes by reference... (crazy design choice!)
        __hydra_dual__
        Evolve( Evolve const& other):
        fN(other.fN),
        fStep(other.fStep),
        fNPX(other.fNPX),
        fNPY(other.fNPY),
        fNPZ(other.fNPZ),
        xvec(other.xvec),
        yvec(other.yvec),
        zvec(other.zvec),
        xvecef(other.xvecef),
        yvecef(other.yvecef),
        zvecef(other.zvecef),
        vecem(other.vecem),
        vechm(other.vechm)
        { }

        
        __hydra_dual__
        Evolve operator=( Evolve const& other){
            if(this == &other) return *this;

            fN=other.fN;
            fStep=other.fStep;
            fNPX=other.fNPX;
            fNPY=other.fNPY;
            fNPZ=other.fNPZ;
            xvec=other.xvec;
            yvec=other.yvec;
            zvec=other.zvec;
            xvecef=other.xvecef;
            yvecef=other.yvecef;
            zvecef=other.zvecef;
            vecem=other.vecem;
            vechm=other.vechm;
            return *this;
        }
        
        //function call operator needs to be callable in host and device sides 
        // and be constant
        template<typename Particle>
        __hydra_dual__
        void operator()(Particle p){
            
            //check if particle went outside the volume
            if(hydra::get<_tc_isin>(p)<0) return;
            
            //get particle charge and initial position
            double charge=hydra::get<_tc_charge>(p);
            double mobility=0;
            double x=hydra::get<_tc_x>(p);
            double y=hydra::get<_tc_y>(p);
            double z=hydra::get<_tc_z>(p);
            
            //get diffusion
            double vdiffx = hydra::get<_tc_gauss_x>(p);
            double vdiffy = hydra::get<_tc_gauss_y>(p);
            double vdiffz = hydra::get<_tc_gauss_z>(p);
            
            double vvec[3];
            size_t ind8[8];
            
            getallindices(vvec, ind8, xvec, fNPX, yvec, fNPY, zvec, fNPZ, x, y, z);
            if(charge<0.) mobility=interpol(vvec, vecem, ind8);
            else mobility=interpol(vvec, vechm, ind8);
            
            //get electric field
            double ex = interpol(vvec, xvecef, ind8);
            double ey = interpol(vvec, yvecef, ind8);
            double ez = interpol(vvec, zvecef, ind8);
            
            //first step
            double k1x=charge*mobility*ex;
            double k1y=charge*mobility*ey;
            double k1z=charge*mobility*ez;
            
            //get third mobility
            x = x+k1x*fStep*0.25;
            y = y+k1x*fStep*0.25;
            z = z+k1x*fStep*0.25;
            
            //get indices and versor
            getallindices(vvec, ind8, xvec, fNPX, yvec, fNPY, zvec, fNPZ, x, y, z);
            if(charge<0.) mobility=interpol(vvec, vecem, ind8);
            else mobility=interpol(vvec, vechm, ind8);
            
            //get electric field
            ex = interpol(vvec, xvecef, ind8);
            ey = interpol(vvec, yvecef, ind8);
            ez = interpol(vvec, zvecef, ind8);

            //second step
            double k2x=charge*mobility*ex;
            double k2y=charge*mobility*ey;
            double k2z=charge*mobility*ez;
            
            x = x+k2x*fStep*0.25;
            y = y+k2x*fStep*0.25;
            z = z+k2x*fStep*0.25;
            
            //get indices and versor
            getallindices(vvec, ind8, xvec, fNPX, yvec, fNPY, zvec, fNPZ, x, y, z);
            if(charge<0.) mobility=interpol(vvec, vecem, ind8);
            else mobility=interpol(vvec, vechm, ind8);
            
            //get electric field
            ex = interpol(vvec, xvecef, ind8);
            ey = interpol(vvec, yvecef, ind8);
            ez = interpol(vvec, zvecef, ind8);

            //third step
            double k3x=charge*mobility*ex;
            double k3y=charge*mobility*ey;
            double k3z=charge*mobility*ez;
            
            
            x = x+k3x*fStep;
            y = y+k3x*fStep;
            z = z+k3x*fStep;
            
            //get indices and versor
            getallindices(vvec, ind8, xvec, fNPX, yvec, fNPY, zvec, fNPZ, x, y, z);
            if(charge<0.) mobility=interpol(vvec, vecem, ind8);
            else mobility=interpol(vvec, vechm, ind8);
            
            //get electric field
            ex = interpol(vvec, xvecef, ind8);
            ey = interpol(vvec, yvecef, ind8);
            ez = interpol(vvec, zvecef, ind8);
            
            //fourth step
            double k4x=charge*mobility*ex;
            double k4y=charge*mobility*ey;
            double k4z=charge*mobility*ez;
            
            
            //evolution
            double Vdriftx=0.1666667*(k1x+2*k2x+2*k3x+k4x)+vdiffx;
            double Vdrifty=0.1666667*(k1y+2*k2y+2*k3y+k4y)+vdiffy;
            double Vdriftz=0.1666667*(k1z+2*k2z+2*k3z+k4z)+vdiffz;
            
//             Vdrift=sqrt((Vdriftx*Vdriftx+Vdrifty*Vdrifty+Vdriftz*Vdriftz));
            
//             MaxV=MAXVDRIFTE;
//             if(charge>0) MaxV=MAXVDRIFTH;
//             
//             if(Vdrift>MaxV){
//                 Vdriftx*=(MaxV/(Vdrift));
//                 Vdrifty*=(MaxV/(Vdrift));
//                 Vdriftz*=(MaxV/(Vdrift));
//             }
            
            
            if(x<xvec[0] || x>xvec[fNPX-1] || y<yvec[0] || y>yvec[fNPY-1] || z<zvec[0] || z>zvec[fNPZ-1]){
                hydra::get<_tc_isin>(p)=-1;
            }

            hydra::get<_tc_x>(p)=x+Vdriftx*fStep;
            hydra::get<_tc_y>(p)=y+Vdrifty*fStep;
            hydra::get<_tc_z>(p)=z+Vdriftz*fStep;
            
        }
        
        
        
        size_t fN;
        double fStep;
        size_t fNPX = 0;
        size_t fNPY = 0;
        size_t fNPZ = 0;
        double *xvecef = nullptr;
        double *yvecef = nullptr;
        double *zvecef = nullptr;
        double *vecem = nullptr;
        double *vechm = nullptr;
        double *xvec = nullptr;
        double *yvec = nullptr;
        double *zvec = nullptr;
        double *in_em = nullptr;
        double *in_hm = nullptr;
    };
    
    struct Evolve_multi
    {

        Evolve_multi()=delete;                                         // avoid creating objects with undefined data
        
        // No __hydra__dual__ here. Only allow objects to be constructed in host side. 
        Evolve_multi(size_t n, double timestep, VecDev_t<double> &inef_x, VecDev_t<double> &inef_y, VecDev_t<double> &inef_z, VecDev_t<double> &in_xef, VecDev_t<double> &in_yef, VecDev_t<double> &in_zef, VecDev_t<double> &inem_x, VecDev_t<double> &inem_y, VecDev_t<double> &inem_z, VecDev_t<double> &in_em, VecDev_t<double> &inhm_x, VecDev_t<double> &inhm_y, VecDev_t<double> &inhm_z, VecDev_t<double> &in_hm, double _xmin, double _xmax, double _ymin, double _ymax, double _zmin, double _zmax):
        fN(n),
        fStep(timestep),
        xmin(_xmin),
        xmax(_xmax),
        ymin(_ymin),
        ymax(_ymax),
        zmin(_zmin),
        zmax(_zmax),
        fNPXef(inef_x.size()),
        fNPYef(inef_y.size()),
        fNPZef(inef_z.size()),
        xefvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inef_x.data())),
        yefvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inef_y.data())),
        zefvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inef_z.data())),
        fNPXem(inem_x.size()),
        fNPYem(inem_y.size()),
        fNPZem(inem_z.size()),
        xemvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inem_x.data())),
        yemvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inem_y.data())),
        zemvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inem_z.data())),
        fNPXhm(inhm_x.size()),
        fNPYhm(inhm_y.size()),
        fNPZhm(inhm_z.size()),
        xhmvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inhm_x.data())),
        yhmvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inhm_y.data())),
        zhmvec(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(inhm_z.data())),
        xvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_xef.data())),
        yvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_yef.data())),
        zvecef(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_zef.data())),
        vecem(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_em.data())),
        vechm(HYDRA_EXTERNAL_NS::thrust::raw_pointer_cast(in_hm.data()))
        {}
        
        // Copy constructor needs be __hydra_dual__ to be copyable
        // on device and host sides, in the time of thread distribution.
        // hydra, as well as thrust and std pass arguments by value.
        // Obs.: RooFit passes by reference... (crazy design choice!)
        __hydra_dual__
        Evolve_multi( Evolve_multi const& other):
        fN(other.fN),
        fStep(other.fStep),
        xmin(other.xmin),
        xmax(other.xmax),
        ymin(other.ymin),
        ymax(other.ymax),
        zmin(other.zmin),
        zmax(other.zmax),
        fNPXef(other.fNPXef),
        fNPYef(other.fNPYef),
        fNPZef(other.fNPZef),
        xefvec(other.xefvec),
        yefvec(other.yefvec),
        zefvec(other.zefvec),
        fNPXem(other.fNPXem),
        fNPYem(other.fNPYem),
        fNPZem(other.fNPZem),
        xemvec(other.xemvec),
        yemvec(other.yemvec),
        zemvec(other.zemvec),
        fNPXhm(other.fNPXhm),
        fNPYhm(other.fNPYhm),
        fNPZhm(other.fNPZhm),
        xhmvec(other.xhmvec),
        yhmvec(other.yhmvec),
        zhmvec(other.zhmvec),
        xvecef(other.xvecef),
        yvecef(other.yvecef),
        zvecef(other.zvecef),
        vecem(other.vecem),
        vechm(other.vechm)
        { }

        
        __hydra_dual__
        Evolve_multi operator=( Evolve_multi const& other){
            if(this == &other) return *this;

            fN=other.fN;
            fStep=other.fStep;
            xmin=other.xmin;
            xmax=other.xmax;
            ymin=other.ymin;
            ymax=other.ymax;
            zmin=other.zmin;
            zmax=other.zmax;
            fNPXef=other.fNPXef;
            fNPYef=other.fNPYef;
            fNPZef=other.fNPZef;
            xefvec=other.xefvec;
            yefvec=other.yefvec;
            zefvec=other.zefvec;
            fNPXem=other.fNPXem;
            fNPYem=other.fNPYem;
            fNPZem=other.fNPZem;
            xemvec=other.xemvec;
            yemvec=other.yemvec;
            zemvec=other.zemvec;
            fNPXef=other.fNPXef;
            fNPYef=other.fNPYef;
            fNPZef=other.fNPZef;
            xefvec=other.xefvec;
            yefvec=other.yefvec;
            zefvec=other.zefvec;
            xvecef=other.xvecef;
            yvecef=other.yvecef;
            zvecef=other.zvecef;
            vecem=other.vecem;
            vechm=other.vechm;
            return *this;
        }
        
        //function call operator needs to be callable in host and device sides 
        // and be constant
        template<typename Particle>
        __hydra_dual__
        void operator()(Particle p){
            
            //check if particle went outside the volume
            if(hydra::get<_tc_isin>(p)<0) return;
            
            //get particle charge and initial position
            double charge=hydra::get<_tc_charge>(p);
            double mobility=0;
            double x=hydra::get<_tc_x>(p);
            double y=hydra::get<_tc_y>(p);
            double z=hydra::get<_tc_z>(p);
            
//             std::cout << p << std::endl;
            
            //get diffusion
            double vdiffx = hydra::get<_tc_gauss_x>(p);
            double vdiffy = hydra::get<_tc_gauss_y>(p);
            double vdiffz = hydra::get<_tc_gauss_z>(p);
            
            double vvec[3];
            size_t ind8[8];
            
            if(charge<0.){
                getallindices(vvec, ind8, xemvec, fNPXem, yemvec, fNPYem, zemvec, fNPZem, x, y, z);
                mobility=interpol(vvec, vecem, ind8);
            }
            else{
                getallindices(vvec, ind8, xhmvec, fNPXhm, yhmvec, fNPYhm, zhmvec, fNPZhm, x, y, z);
                mobility=interpol(vvec, vechm, ind8);
            }
            
            //get electric field
            getallindices(vvec, ind8, xefvec, fNPXef, yefvec, fNPYef, zefvec, fNPZef, x, y, z);
            double ex = interpol(vvec, xvecef, ind8);
            double ey = interpol(vvec, yvecef, ind8);
            double ez = interpol(vvec, zvecef, ind8);
            
            //first step
            double k1x=charge*mobility*ex;
            double k1y=charge*mobility*ey;
            double k1z=charge*mobility*ez;
            
            //get third mobility
            x = x+k1x*fStep*0.25;
            y = y+k1x*fStep*0.25;
            z = z+k1x*fStep*0.25;
            
            //get indices and versor
            if(charge<0.){
                getallindices(vvec, ind8, xemvec, fNPXem, yemvec, fNPYem, zemvec, fNPZem, x, y, z);
                mobility=interpol(vvec, vecem, ind8);
            }
            else{
                getallindices(vvec, ind8, xhmvec, fNPXhm, yhmvec, fNPYhm, zhmvec, fNPZhm, x, y, z);
                mobility=interpol(vvec, vechm, ind8);
            }
            
            //get electric field
            getallindices(vvec, ind8, xefvec, fNPXef, yefvec, fNPYef, zefvec, fNPZef, x, y, z);
            ex = interpol(vvec, xvecef, ind8);
            ey = interpol(vvec, yvecef, ind8);
            ez = interpol(vvec, zvecef, ind8);
            
            //second step
            double k2x=charge*mobility*ex;
            double k2y=charge*mobility*ey;
            double k2z=charge*mobility*ez;
            
            x = x+k2x*fStep*0.25;
            y = y+k2x*fStep*0.25;
            z = z+k2x*fStep*0.25;
            
            //get indices and versor
            if(charge<0.){
                getallindices(vvec, ind8, xemvec, fNPXem, yemvec, fNPYem, zemvec, fNPZem, x, y, z);
                mobility=interpol(vvec, vecem, ind8);
            }
            else{
                getallindices(vvec, ind8, xhmvec, fNPXhm, yhmvec, fNPYhm, zhmvec, fNPZhm, x, y, z);
                mobility=interpol(vvec, vechm, ind8);
            }
            
            //get electric field
            getallindices(vvec, ind8, xefvec, fNPXef, yefvec, fNPYef, zefvec, fNPZef, x, y, z);
            ex = interpol(vvec, xvecef, ind8);
            ey = interpol(vvec, yvecef, ind8);
            ez = interpol(vvec, zvecef, ind8);
            
            //third step
            double k3x=charge*mobility*ex;
            double k3y=charge*mobility*ey;
            double k3z=charge*mobility*ez;
            
            
            x = x+k3x*fStep;
            y = y+k3x*fStep;
            z = z+k3x*fStep;
            
            //get indices and versor
            if(charge<0.){
                getallindices(vvec, ind8, xemvec, fNPXem, yemvec, fNPYem, zemvec, fNPZem, x, y, z);
                mobility=interpol(vvec, vecem, ind8);
            }
            else{
                getallindices(vvec, ind8, xhmvec, fNPXhm, yhmvec, fNPYhm, zhmvec, fNPZhm, x, y, z);
                mobility=interpol(vvec, vechm, ind8);
            }
            
            //get electric field
            getallindices(vvec, ind8, xefvec, fNPXef, yefvec, fNPYef, zefvec, fNPZef, x, y, z);
            ex = interpol(vvec, xvecef, ind8);
            ey = interpol(vvec, yvecef, ind8);
            ez = interpol(vvec, zvecef, ind8);
            
            //fourth step
            double k4x=charge*mobility*ex;
            double k4y=charge*mobility*ey;
            double k4z=charge*mobility*ez;
            
            
            //evolution
            double Vdriftx=0.1666667*(k1x+2*k2x+2*k3x+k4x)+vdiffx;
            double Vdrifty=0.1666667*(k1y+2*k2y+2*k3y+k4y)+vdiffy;
            double Vdriftz=0.1666667*(k1z+2*k2z+2*k3z+k4z)+vdiffz;
            
//             Vdrift=sqrt((Vdriftx*Vdriftx+Vdrifty*Vdrifty+Vdriftz*Vdriftz));
            
//             MaxV=MAXVDRIFTE;
//             if(charge>0) MaxV=MAXVDRIFTH;
//             
//             if(Vdrift>MaxV){
//                 Vdriftx*=(MaxV/(Vdrift));
//                 Vdrifty*=(MaxV/(Vdrift));
//                 Vdriftz*=(MaxV/(Vdrift));
//             }
            
            
            if(x<xmin || x>xmax || y<ymin || y>ymax || z<zmin || z>zmax){
                hydra::get<_tc_isin>(p)=-1;
            }
//             std::cout << Vdriftx << "\t" << Vdrifty << "\t" << Vdriftz << std::endl;
            hydra::get<_tc_x>(p)=x+Vdriftx*fStep;
            hydra::get<_tc_y>(p)=y+Vdrifty*fStep;
            hydra::get<_tc_z>(p)=z+Vdriftz*fStep;
            
//             std::cout << fStep << std::endl;
//             std::cout << Vdriftx << "\t" << Vdrifty << "\t" << Vdriftz << std::endl;
//             std::cout << p << std::endl;
//             std::cout << "gen"<<std::endl;
        }
        
        
        
        size_t fN;
        double fStep;
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double zmin;
        double zmax;
        size_t fNPXef = 0;
        size_t fNPYef = 0;
        size_t fNPZef = 0;
        size_t fNPXem = 0;
        size_t fNPYem = 0;
        size_t fNPZem = 0;
        size_t fNPXhm = 0;
        size_t fNPYhm = 0;
        size_t fNPZhm = 0;
        double *xvecef = nullptr;
        double *yvecef = nullptr;
        double *zvecef = nullptr;
        double *vecem  = nullptr;
        double *vechm  = nullptr;
        double *xefvec   = nullptr;
        double *yefvec   = nullptr;
        double *zefvec   = nullptr;
        double *xemvec   = nullptr;
        double *yemvec   = nullptr;
        double *zemvec   = nullptr;
        double *xhmvec   = nullptr;
        double *yhmvec   = nullptr;
        double *zhmvec   = nullptr;
    };
}
#endif
