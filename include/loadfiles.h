/*----------------------------------------------------------------------------
 * 
 *   Copyright (C) 2018-2020 Andrea Contu e Angelo Loi
 *
 *   This file is part of TCode software.
 *
 *   TCode is free software: you can redistribute it and/or dify
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
 * loadfiles.h
 *
 *  Created on: 12/11/2018
 *      Author: Andrea Contu
 */


//Functions to load files

using namespace hydra::placeholders;

namespace loaddata{

    
    RunningStateHost_t getDummy(size_t nparticles, double length, size_t group=DEFAULT_GROUP){
        INFO_LINE("Generating deposit")
        hydra::Random<> Generator( std::chrono::system_clock::now().time_since_epoch().count());
        size_t np=(nparticles<=0) ? MAXPARTICLES : nparticles;

        if(np==0) ERROR_LINE("Number of particles is zero! Check your configuration.")
        if(np==1){
            if(group >1 && np/group==0) ERROR_LINE("Number of particles (or bunches of particles) is zero! Check your configuration.")
            auto data_d1=RunningStateHost_t(np,RunningTuple_t(std::copysign(1.,length)*group,0.,0.,0., 1.,0.,0.,0.,0.,0.,0.,0.));
            return data_d1;
        }
        size_t halfsize=np/2;
        if(group>1) halfsize/=group;
        RunningStateHost_t data_d(halfsize*2);
        
        
//------------------------------------------------------------------------------------------
        int save_group = group;
            //need the float in order to compute the lost charge in the calculus
        double idep_f = deposit/group;
            //lost charge will be divided by number of cycles and multiplied by group. total charge is conv
        double rest_charge = ((idep_f - idep)/idep)*group;
            //if division is less than 1, correct in order to not loose the charge!
        if(idep_f<1){
            idep=1;
            rest_charge = idep_f*group;
            group=0;
        }
//------------------------------------------------------------------------------------------
        
        
        auto data_de=RunningStateHost_t(halfsize,RunningTuple_t(-1.*group,0.,0.,0., 1.,0.,0.,0.,0.,0.,0.,0.));
        auto data_dh=RunningStateHost_t(halfsize,RunningTuple_t(1.*group,0.,0.,0., 1.,0.,0.,0.,0.,0.,0.,0.));
        
        
        //distribute them randomly in a line along Z
        if(length>0.){
            int seed=0;
            Generator.SetSeed(seed); // IMPORTANT, otherwise the seed stays the same
            Generator.Uniform(0., length, data_de.begin(_tc_z), data_de.end(_tc_z));
            Generator.SetSeed(seed);
            Generator.Uniform(0., length, data_dh.begin(_tc_z), data_dh.end(_tc_z));
        }
        hydra::copy(data_de,hydra::make_range(data_d.begin(),data_d.begin()+halfsize));
        hydra::copy(data_dh,hydra::make_range(data_d.begin()+halfsize,data_d.end()));
        
        return data_d;
    }

C@mpid@n@1
    
    /**********************************************************************************
    Load energy deposit
    ***********************************************************************************/
//     RunningStateHost_t loadfile(std::string path, size_t group=DEFAULT_GROUP){
//         std::ifstream infile;
//         
//         RunningStateHost_t data_d;
// //         if(path=="dummy"){
// //             return getDummy(MAXPARTICLES);
// //         }
//         
//         infile.open(path.c_str());// file containing numbers in 3 columns 
//         
//         
//         if(infile.fail()){ 
//             ERROR_LINE("Error loading file " << path << ". Please check!")
//             exit(1);
//         }
//     
//         double pix=0;
//         double piy=0;
//         double piz=0;
//         double pfx=0;
//         double pfy=0;
//         double pfz=0;
//         double typ=0;
//         
//         std::string line;
//         
//         
//         double deposit=0;
//         int ct=0;
//         
// //         //new code
//         while(std::getline(infile,line)){ 
//             typ=0;
//             std::stringstream lineStream;
//             lineStream << line;
//             lineStream >> pix >> piy >> piz >> pfx >> pfy >> pfz >> deposit >> typ;
//     //         std::cout << line << std::endl;
//             if(deposit<1) continue;
//             
//             unsigned int idep=deposit/group;
//             
//             /***************************************************************************
//                                         ***NEW FEATURE START-dode***
//             **************************************************************************/
//             int save_group = group;
//             //need the float in order to compute the lost charge in the calculus
//             double idep_f = deposit/group;
//             //lost charge will be divided by number of cycles and multiplied by group. total charge is conv
//             double rest_charge = ((idep_f - idep)/idep)*group;
//             //if division is less than 1, correct in order to not loose the charge!
//             std::cout<<"Before:"<<idep<<" "<<idep_f<<" "<<rest_charge<<" "<<deposit<<"\n";
//             if(idep_f<1){
//                     idep=1;
//                     rest_charge = idep_f*group;
//                     group=0;
//             }
//             std::cout<<"After:"<<idep<<" "<<idep_f<<" "<<rest_charge<<"\n";
//             
//             /***************************************************************************
//                                         ***NEW FEATURE END-dode***
//             **************************************************************************/
//             
//             if(idep<1) continue;
//     //         std::cout << pix << "\t" << piy << "\t" << piz << "\t" << pfx << "\t" << pfy << "\t" << pfz << "\t" << deposit << std::endl;
//             double d=sqrt((pfx-pix)*(pfx-pix)+(pfy-piy)*(pfy-piy)+(pfz-piz)*(pfz-piz));
//             
//             if(d==0){
//                 for(unsigned int i=0;i<idep;i++){
//                     data_d.push_back(RunningTuple_t(-1.*(group+rest_charge),pix, piy, piz, 1.,0.,0.,0.,0.,0.,0.,typ));
//                     data_d.push_back(RunningTuple_t(1.*(group+rest_charge),pix, piy, piz, 1.,0.,0.,0.,0.,0.,0.,typ));
//                 }
//             }
//             else{
//                 double a=(pfx-pix)/d;
//                 double b=(pfy-piy)/d;
//                 double c=(pfz-piz)/d;
//                 
//                 std::random_device rd;
//                 std::mt19937 mt(rd());
//                 std::uniform_real_distribution<double> dist(0.,d);
//                 for(unsigned int i=0;i<idep;i++){
//                     double vec=dist(mt);
//                     data_d.push_back(RunningTuple_t(-1.*(group+rest_charge),a*vec+pix, b*vec+piy, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));
//                     data_d.push_back(RunningTuple_t(1.*(group+rest_charge),a*vec+pix, b*vec+piy, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));
//                 }
// //                 double vec=dist(mt);
// //                 data_d.push_back(RunningTuple_t(-1.*(rest_charge),a*vec+pix, b*vec+piy, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));
// //                 data_d.push_back(RunningTuple_t(1.*(rest_charge),a*vec+pix, b*vec+piy, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));
//                 
//             }
// 
//             ct++;
//             if(group == 0){
//                 group=save_group;
//             }
//         } 
//         infile.close(); 
//     //     std::cout << "Size: "<< data_d.size()<< std::endl;
//         
//         
//         return data_d;
//     }
    
// // //       RunningStateHost_t loadfile(std::string path, size_t group=DEFAULT_GROUP){
// // //         std::ifstream infile;
// // //         
// // //         RunningStateHost_t data_d;
// // // //         if(path=="dummy"){
// // // //             return getDummy(MAXPARTICLES);
// // // //         }
// // //         
// // //         infile.open(path.c_str());// file containing numbers in 3 columns 
// // //         
// // //         
// // //         if(infile.fail()){ 
// // //             ERROR_LINE("Error loading file " << path << ". Please check!")
// // //             exit(1);
// // //         }
// // //     
// // //         double pix=0;
// // //         double piy=0;
// // //         double piz=0;
// // //         double pfx=0;
// // //         double pfy=0;
// // //         double pfz=0;
// // //         double typ=0;
// // //        
// // //         double sigm_f=0;
// // //         double sigm_i=0;
// // //         
// // //         std::string line;
// // //         
// // //         
// // //         double deposit=0;
// // //         int ct=0;
// // //         
// // // //         //new code
// // //         while(std::getline(infile,line)){ 
// // //             typ=0;
// // //             std::default_random_engine generator;
// // //             std::stringstream lineStream;
// // //             lineStream << line;
// // //             lineStream >> pix >> piy >> piz >> pfx >> pfy >> pfz >> deposit >> typ >> sigm_i >> sigm_f;
// // //     //         std::cout << line << std::endl;
// // //             if(deposit<1) continue;
// // //             
// // //             unsigned int idep=deposit/group;
// // //             
// // //             /***************************************************************************
// // //                                         ***NEW FEATURE START-dode***
// // //             **************************************************************************/
// // //             int save_group = group;
// // //             //need the float in order to compute the lost charge in the calculus
// // //             double idep_f = deposit/group;
// // //             //lost charge will be divided by number of cycles and multiplied by group. total charge is conv
// // //             double rest_charge = ((idep_f - idep)/idep)*group;
// // //             //if division is less than 1, correct in order to not loose the charge!
// // //             std::cout<<"Before:"<<idep<<" "<<idep_f<<" "<<rest_charge<<" "<<deposit<<"\n";
// // //             if(idep_f<1){
// // //                     idep=1;
// // //                     rest_charge = idep_f*group;
// // //                     group=0;
// // //             }
// // //             std::cout<<"After:"<<idep<<" "<<idep_f<<" "<<rest_charge<<"\n";
// // //             
// // //             /***************************************************************************
// // //                                         ***NEW FEATURE END-dode***
// // //             **************************************************************************/
// // //             
// // //             if(idep<1) continue;
// // //     //         std::cout << pix << "\t" << piy << "\t" << piz << "\t" << pfx << "\t" << pfy << "\t" << pfz << "\t" << deposit << std::endl;
// // //             double d=sqrt((pfx-pix)*(pfx-pix)+(pfy-piy)*(pfy-piy)+(pfz-piz)*(pfz-piz));
// // //             
// // //             if(d==0){
// // //                 for(unsigned int i=0;i<idep;i++){
// // //                     data_d.push_back(RunningTuple_t(-1.*(group+rest_charge),pix, piy, piz, 1.,0.,0.,0.,0.,0.,0.,typ));
// // //                     data_d.push_back(RunningTuple_t(1.*(group+rest_charge),pix, piy, piz, 1.,0.,0.,0.,0.,0.,0.,typ));
// // //                 }
// // //             }
// // //             else{
// // //                 double a=(pfx-pix)/d;
// // //                 double b=(pfy-piy)/d;
// // //                 double c=(pfz-piz)/d;
// // //                 
// // //                 double sigma_a=(sigm_f-sigm_i)/d;
// // //                 
// // //                 std::random_device rd;
// // //                 std::mt19937 mt(rd());
// // //                 std::uniform_real_distribution<double> dist(0.,d);
// // //                 for(unsigned int i=0;i<idep;i++){
// // //                     double vec=dist(mt);
// // //                     std::normal_distribution<double> Xdistribution(a*vec+pix,sigma_a*vec+sigm_i);
// // //                     std::normal_distribution<double> Ydistribution(b*vec+piy,sigma_a*vec+sigm_i);
// // //                     double pos_X = Xdistribution(mt);
// // //                     double pos_Y = Ydistribution(mt);
// // //                     data_d.push_back(RunningTuple_t(-1.*(group+rest_charge),pos_X, pos_Y, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));
// // //                     data_d.push_back(RunningTuple_t(1.*(group+rest_charge),pos_X, pos_Y, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));                    
// // //                     
// // // //                     data_d.push_back(RunningTuple_t(-1.*(group+rest_charge),a*vec+pix, b*vec+piy, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));
// // // //                     data_d.push_back(RunningTuple_t(1.*(group+rest_charge),a*vec+pix, b*vec+piy, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));
// // //                 }
// // // //                 double vec=dist(mt);
// // // //                 data_d.push_back(RunningTuple_t(-1.*(rest_charge),a*vec+pix, b*vec+piy, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));
// // // //                 data_d.push_back(RunningTuple_t(1.*(rest_charge),a*vec+pix, b*vec+piy, c*vec+piz, 1.,0.,0.,0.,0.,0.,0.,typ));
// // //                 
// // //             }
// // // 
// // //             ct++;
// // //             if(group == 0){
// // //                 group=save_group;
// // //             }
// // //         } 
// // //         infile.close(); 
// // //     //     std::cout << "Size: "<< data_d.size()<< std::endl;
// // //         
// // //         
// // //         return data_d;
// // //     }
    
    //prime difiche per LASER  ----------------- START --------------------------------
    
    RunningStateHost_t loadfile(std::string path, size_t group=DEFAULT_GROUP){
        std::ifstream infile;
        
        RunningStateHost_t data_d;  //data vector containing carrier positions
        
        infile.open(path.c_str());  // file containing numbers in 3 columns 
        
        
        if(infile.fail()){ 
            ERROR_LINE("Error loading file " << path << ". Please check!")
            exit(1);
        }
    
        double pix=0;                       //initialising coordinate
        double piy=0;                       //initialising coordinate
        double piz=0;                       //initialising coordinate
        double pfx=0;                       //initialising coordinate
        double pfy=0;                       //initialising coordinate
        double pfz=0;                       //initialising coordinate
        double typ=0;                       //initialising particle type
       
        double sigm_f=0;                    //initialising dispersion
        double sigm_i=0;                    //initialising dispersion
        char dispersion[20];                //initialising dispersion TYPE
        double absorption = 0;              //initialising absorption 
        double deposit=0;                   //initialising deposit
        int ct=0;                           //initialising couple number
        
        std::string line;                   //initialising line
        while(std::getline(infile,line)){ 
            //reading line
            std::stringstream lineStream;
            lineStream << line;
            lineStream >> pix >> piy >> piz >> pfx >> pfy >> pfz >> deposit >> typ >> sigm_i >> sigm_f >> dispersion >> absorption;
    //         std::cout << line << std::endl;
            if(deposit<1) continue;
            
            unsigned int idep=deposit/group;
            
            /***************************************************************************
                                        ***NEW FEATURE START-code***
            **************************************************************************/
            int save_group = group;
            //need the float in order to compute the lost charge in the calculus
            double idep_f = deposit/group;
            //lost charge will be divided by number of cycles and multiplied by group. total charge is conv
            double rest_charge = ((idep_f - idep)/idep)*group;
            //if division is less than 1, correct in order to not loose the charge!
            if(idep_f<1){
                    idep=1;
                    rest_charge = idep_f*group;
                    group=0;
            }
            
            /***************************************************************************
                                        ***NEW FEATURE END-code***
            **************************************************************************/
            
            if(idep<1) continue;
            double d=sqrt((pfx-pix)*(pfx-pix)+(pfy-piy)*(pfy-piy)+(pfz-piz)*(pfz-piz));
            
            if(d==0){
                for(unsigned int i=0;i<idep;i++){
                    data_d.push_back(RunningTuple_t(-1.*(group+rest_charge),pix, piy, piz, 1.,0.,0.,0.,0.,0.,0.,typ));
                    data_d.push_back(RunningTuple_t(1.*(group+rest_charge),pix, piy, piz, 1.,0.,0.,0.,0.,0.,0.,typ));
                }
            }
            else{
                double a=(pfx-pix)/d;                                               // parameter for linear variation of coordiante
                double b=(pfy-piy)/d;                                               // parameter for linear variation of coordiante
                double c=(pfz-piz)/d;                                               // parameter for linear variation of coordiante
                double sigma_a=(sigm_f-sigm_i)/d;                                   // parameter for linear variation of dispersion
                double vec=0;                                                       // position along line
                double disp_radius =0;                                              // position along dispersion
                //random device for track-----------|
                std::random_device rd;                                              // Random device object
                std::mt19937 mt(rd());                                              // Random device number generator
                //random device for angle-----------|
                std::random_device ang_d;
                std::mt19937 ang_mt(ang_d());
                double phi, theta, gamma;
                
                for(unsigned int i=0; i<idep ; i++){
                    //particle distribution along track
                    //###############################################################################################################################
                    if(absorption == 0){
                        std::uniform_real_distribution<double> dist(0.,d);          //uniform real number distribution
                        vec=dist(mt);                                               //random position along line for point generation
                    }
                    else if(absorption != 0){                                       //If light is used, absorption becomes important. 
                        std::exponential_distribution<double> dist(absorption);     //exponential distribution uses lambert-beer-like equation
                        vec=dist(mt);                                               //to redistribute points along tracks
                        while( vec>d){                                              //if exponential is larger than module, a while loop tries again
                            vec=dist(mt);                                           //until vec<d
                        }
                    }
                    //particle distribution radial to track
                    //###############################################################################################################################
                    
                    std::string str1 ("UNIFORM");
                    std::string str2 ("GAUSSIAN");
                    if(str1.compare(dispersion) == 0){
                        std::uniform_real_distribution<double> Radius_distribution(0,sigma_a*vec+sigm_i);    //Introduces UNIFORM Distribution for radial dispersion of carriers
                        disp_radius = Radius_distribution(mt);
                    }
                    else if(str2.compare(dispersion) == 0){
                        std::normal_distribution<double> Radius_distribution(0,sigma_a*vec+sigm_i);          //Introduces Gaussian Distribution for radial dispersion of carriers
                        disp_radius = Radius_distribution(mt);
                    }
                    else{                                                                                    //If you set an unexisting dispersion shape
                            ERROR_LINE("Dispersion shape does not exist! Set Dispersion at 0!");
                            disp_radius = 0;
                    }
                    //RANDOM DISTRIBUTION ARROUND TRACK
                    //###############################################################################################################################
                        
                    std::uniform_real_distribution<double> angle(0.,2*PI);
                    phi = acos(((pfz-piz)/d));                                      //determines angulation of track in space
                    if(pfy == piy && pfx == pix){                                   //in spherical coordinates
                        theta = 0;                                                  //if theta is undetermined, set by default to 0 ...
                    }
                    else{
                        theta = atan(abs(pfy-piy)/(pfx-pix));                       //... or direct computed                                    
                    }
                    
                    double x_Y = disp_radius*sin(theta + PI/2)*cos(phi);            //defines perpendicular dispersion vector to line
                    double y_Y = disp_radius*sin(theta + PI/2)*sin(phi);    
                    double z_Y = disp_radius*cos(theta + PI/2);
                    
                    gamma = angle(ang_d);                                           // random rotation arround track
                    //final point obtained by rotating the dispersion vector arround the track by gamma
                    double dispX = (((pfx - pix)*(pfx - pix)*(1/(d*d)) + (1-(pfx - pix)*(pfx - pix)*(1/(d*d)))*cos(gamma))*x_Y + ((pfx - pix)*(pfy - piy)*(1/(d*d))*(1-cos(gamma)) - (pfz - piz)*(1/(d))*sin(gamma))*y_Y + ((pfx - pix)*(pfz - piz)*(1/(d*d))*(1 - cos(gamma)) + (pfy - piy)*(1/(d))*sin(gamma))*z_Y);
                    
                    double dispY = ((pfx - pix)*(1/(d*d))*(pfy - piy)*(1 - cos(gamma))*x_Y + (pfz - piz)*(1/(d))*sin(gamma)*x_Y + (pfy - piy)*(1/(d*d))*(pfy - piy)*y_Y + (1-(pfy - piy)*(1/(d*d))*(pfy - piy))*cos(gamma)*y_Y + (pfy - piy)*(1/(d*d))*(pfz - piz)*(1 - cos(gamma))*z_Y - (pfx - pix)*(1/(d))*sin(gamma)*z_Y);
                    
                    double dispZ = ((pfx - pix)*(pfz - piz)*(1/(d*d))*(1-cos(gamma))*x_Y - (pfy - piy)*(1/(d))*sin(gamma)*x_Y + (pfy - piy)*(1/(d*d))*(pfz - piz)*(1-cos(gamma))*y_Y + (pfx - pix)*(1/(d*d))*sin(gamma)*y_Y + (pfz - piz)*(pfz - piz)*(1/(d*d))*z_Y + (1-(pfz - piz)*(1/(d*d))*(pfz - piz))*cos(gamma)*z_Y);
                 
                    //adds poins
                    data_d.push_back(RunningTuple_t(-1.*(group+rest_charge),a*vec+pix+dispX, b*vec+piy+dispY, c*vec+piz+dispZ, 1.,0.,0.,0.,0.,0.,0.,typ));
                    data_d.push_back(RunningTuple_t(1.*(group+rest_charge) ,a*vec+pix+dispX, b*vec+piy+dispY, c*vec+piz+dispZ, 1.,0.,0.,0.,0.,0.,0.,typ));
                }
            }

            ct++;
            if(group == 0){
                group=save_group;
            }
        } 
        infile.close(); 
        
        return data_d;
    }

    
    
    
    
    // functor to transform input deposit
    struct Translate
    {

        Translate()=delete;                                         // avoid creating objects with undefined data

        // No __hydra__dual__ here. Only allow objects to be constructed in host side. 
        Translate(double x, double y, double z, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax):
        fX(x),
        fY(y),
        fZ(z),
        fxmin(xmin),
        fxmax(xmax),
        fymin(ymin),
        fymax(ymax),
        fzmin(zmin),
        fzmax(zmax)
        {}

        // Copy constructor needs be __hydra_dual__ to be copyable
        // on device and host sides, in the time of thread distribution.
        // hydra, as well as thrust and std pass arguments by value.
        // Obs.: RooFit passes by reference... (crazy design choice!)
        __hydra_dual__
        Translate( Translate const& other):
        fX(other.fX),
        fY(other.fY),
        fZ(other.fZ),
        fxmin(other.fxmin),
        fxmax(other.fxmax),
        fymin(other.fymin),
        fymax(other.fymax),
        fzmin(other.fzmin),
        fzmax(other.fzmax)
        { }


        //function call operator needs to be callable in host and device sides 
        // and be constant
        template<typename Particle>
        __hydra_dual__
        void operator()(Particle p){
            
            hydra::get<_tc_x>(p)=hydra::get<_tc_x>(p)+fX;
            hydra::get<_tc_y>(p)=hydra::get<_tc_y>(p)+fY;
            hydra::get<_tc_z>(p)=hydra::get<_tc_z>(p)+fZ;
            
            if(hydra::get<_tc_x>(p)<fxmin || hydra::get<_tc_x>(p)>fxmax || hydra::get<_tc_y>(p)<fymin || hydra::get<_tc_y>(p)>fymax || hydra::get<_tc_z>(p)<fzmin || hydra::get<_tc_z>(p)>fzmax) hydra::get<_tc_isin>(p)=-1;
            
        }
        
        
    // private:
        double fX;
        double fY;
        double fZ;
        double fxmin;
        double fxmax;
        double fymin;
        double fymax;
        double fzmin;
        double fzmax;

    };
}
