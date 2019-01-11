/*----------------------------------------------------------------------------
 * 
 *   Copyright (C) 2018 Andrea Contu e Angelo Loi
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
 * analysis.h
 *
 *  Created on: 12/11/2018
 *      Author: Andrea Contu
 */

//ROOT
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

using namespace hydra::placeholders;
namespace analysis{
    template <typename Data>
    TH3D *getHist3DDraw(Data mapx, Data mapy, Data mapz, std::set<double> _vx, std::set<double> _vy, std::set<double> _vz){
        
        std::vector<double> vx(_vx.begin(),_vx.end());
        std::vector<double> vy(_vy.begin(),_vy.end());
        std::vector<double> vz(_vz.begin(),_vz.end());
        std::vector<double> vex(vx.size()+1);
        std::vector<double> vey(vy.size()+1);
        std::vector<double> vez(vz.size()+1);
        
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
        
        TH3D *h3=new TH3D("physmap","physmap;x [#mum];y [#mum];z [#mum]",vx.size(),vex.data(),vy.size(),vey.data(),vz.size(),vez.data());
        

        size_t count=0;
        for(size_t iz=0;iz<vz.size();iz++){
            for(size_t ix=0;ix<vx.size();ix++){
                for(size_t iy=0;iy<vy.size();iy++){
                    h3->SetBinContent(ix+1,iy+1,iz+1, sqrt(mapx[count]*mapx[count]+mapy[count]*mapy[count]+mapz[count]*mapz[count]));
                    count++;
                }
            }
        }
        return h3;
    }

    TGraph *setGraph(size_t points, std::string nm, int color, int mstyle=8, int msize=1, int lwidth=2){
        TGraph *g = new TGraph(points);
        g->SetName(nm.c_str());
        g->SetMarkerStyle(mstyle);
        g->SetMarkerSize(msize);
        g->SetMarkerColor(color);
        g->SetLineColor(color);
        g->SetLineWidth(lwidth);
        return g;
    }

    struct SelectChargeAndSec
    {
        
        template<typename Particle>
        __hydra_dual__
        ReducedTuple_t operator()(Particle p){
            return ReducedTuple_t(hydra::get<2>(p)*(hydra::get<0>(p)>0)*(hydra::get<1>(p)>0), //holes current
                                hydra::get<2>(p)*(hydra::get<0>(p)<0)*(hydra::get<1>(p)>0), //electrons current
                                hydra::get<2>(p)*(hydra::get<0>(p)>0)*(hydra::get<3>(p)>0)*(hydra::get<1>(p)>0), //holes current secondaries
                                hydra::get<2>(p)*(hydra::get<0>(p)<0)*(hydra::get<3>(p)>0)*(hydra::get<1>(p)>0), //electrons current secondaries
                                1.0*(hydra::get<0>(p)>0)*(hydra::get<1>(p)<0), // count lost holes
                                1.0*(hydra::get<0>(p)<0)*(hydra::get<1>(p)<0), // count lost electrons
                                1.0*(hydra::get<0>(p)>0)*(hydra::get<1>(p)<0)*(hydra::get<3>(p)>0), // count lost holes secondaries
                                1.0*(hydra::get<0>(p)<0)*(hydra::get<1>(p)<0)*(hydra::get<3>(p)>0)  // count lost electrons secondaries
                                                                                        );
        }
        
    };



    struct SumTuples{
    __hydra_dual__ 
    ReducedTuple_t operator()(const ReducedTuple_t &a,
                                                const ReducedTuple_t &b) {
        return ReducedTuple_t(hydra::get<0>(a)+hydra::get<0>(b),
                    hydra::get<1>(a)+hydra::get<1>(b),
                    hydra::get<2>(a)+hydra::get<2>(b),
                    hydra::get<3>(a)+hydra::get<3>(b),
                    hydra::get<4>(a)+hydra::get<4>(b),
                    hydra::get<5>(a)+hydra::get<5>(b),
                    hydra::get<6>(a)+hydra::get<6>(b),
                    hydra::get<7>(a)+hydra::get<7>(b)
            );
        }
    };


    void AnalyseSim(std::vector<ReducedTuple_t> &tp_currs, double timestep, std::map<std::string,std::string> settings){
        
        INFO_LINE("Analysing simulation...")
        
        TFile *outf=new TFile(Form("%s/%s.root",settings["outputdir"].c_str(),settings["name"].c_str()),"RECREATE");
        TTree *tree=new TTree("data","data");
        double t_ec=0, t_hc=0, t_tc=0, t_time=0, t_einduced=0, t_hinduced=0, t_tinduced=0, t_elost=0, t_hlost=0, t_tlost=0;
        double t_ec_sec=0, t_hc_sec=0, t_tc_sec=0, t_einduced_sec=0, t_hinduced_sec=0, t_tinduced_sec=0, t_elost_sec=0, t_hlost_sec=0, t_tlost_sec=0;
        tree->Branch("I_h",&t_hc,"I_h/D");
        tree->Branch("I_e",&t_ec,"I_e/D");
        tree->Branch("I_t",&t_tc,"I_t/D");
        tree->Branch("time",&t_time,"time/D");
        tree->Branch("Qe_ind",&t_einduced,"Qe_ind/D");
        tree->Branch("Qe_lost",&t_elost,"Qe_lost/D");
        tree->Branch("Qh_ind",&t_hinduced,"Qh_ind/D");
        tree->Branch("Qh_lost",&t_hlost,"Qh_lost/D");
        tree->Branch("Qt_ind",&t_tinduced,"Qt_ind/D");
        tree->Branch("Qt_lost",&t_tlost,"Qt_lost/D");
        
        tree->Branch("I_h_sec",&t_hc_sec,"I_h_sec/D");
        tree->Branch("I_e_sec",&t_ec_sec,"I_e_sec/D");
        tree->Branch("I_t_sec",&t_tc_sec,"I_t_sec/D");
        tree->Branch("Qe_ind_sec",&t_einduced_sec,"Qe_ind_sec/D");
        tree->Branch("Qe_lost_sec",&t_elost_sec,"Qe_lost_sec/D");
        tree->Branch("Qh_ind_sec",&t_hinduced_sec,"Qh_ind_sec/D");
        tree->Branch("Qh_lost_sec",&t_hlost_sec,"Qh_lost_sec/D");
        tree->Branch("Qt_ind_sec",&t_tinduced_sec,"Qt_ind_sec/D");
        tree->Branch("Qt_lost_sec",&t_tlost_sec,"Qt_lost_sec/D");
        
        TGraph *ghc = setGraph(tp_currs.size(),"holescur",COLORHOLE);
        TGraph *gec = setGraph(tp_currs.size(),"electronscur",COLORELEC);
        TGraph *gtc = setGraph(tp_currs.size(),"totcur",1);
        
        TGraph *ghind=(TGraph*)ghc->Clone("holeind");
        ghind->SetLineStyle(7);
        TGraph *geind=(TGraph*)gec->Clone("electronind");
        geind->SetLineStyle(7);
        TGraph *gtind=(TGraph*)gtc->Clone("totind");
        gtind->SetLineStyle(7);
        
        TGraph *ghlost=(TGraph*)ghc->Clone("holelost");
        ghlost->SetLineStyle(9);
        TGraph *gelost=(TGraph*)gec->Clone("electronlost");
        gelost->SetLineStyle(9);
        TGraph *gtlost=(TGraph*)gtc->Clone("totlost");
        gtlost->SetLineStyle(9);
        
        TGraph *ghc_sec = setGraph(tp_currs.size(),"holescur_sec",COLORHOLESEC);
        TGraph *gec_sec = setGraph(tp_currs.size(),"electronscur_sec",COLORELECSEC);
        TGraph *gtc_sec = setGraph(tp_currs.size(),"totcur_sec",17);
        
        TGraph *ghind_sec=(TGraph*)ghc_sec->Clone("holeind_sec");
        ghind_sec->SetLineStyle(7);
        TGraph *geind_sec=(TGraph*)gec_sec->Clone("electronind_sec");
        geind_sec->SetLineStyle(7);
        TGraph *gtind_sec=(TGraph*)gtc_sec->Clone("totind_sec");
        gtind_sec->SetLineStyle(7);
        
        TGraph *ghlost_sec=(TGraph*)ghc_sec->Clone("holelost_sec");
        ghlost_sec->SetLineStyle(9);
        TGraph *gelost_sec=(TGraph*)gec_sec->Clone("electronlost_sec");
        gelost_sec->SetLineStyle(9);
        TGraph *gtlost_sec=(TGraph*)gtc_sec->Clone("totlost_sec");
        gtlost_sec->SetLineStyle(9);
        
        double currmax=-999999;
        double chargemax=-999999;
        double currmin=-0;
        double chargemin=-0;
        unsigned pl=0;
        
        for(auto& curr:tp_currs){
            
            double ce=hydra::get<1>(curr);
            double ch=hydra::get<0>(curr);
            double ct=ce+ch;
            double ce_sec=hydra::get<3>(curr);
            double ch_sec=hydra::get<2>(curr);
            double ct_sec=ce_sec+ch_sec;
            
//             std::cout << hydra::get<0>(curr) << "\t" << hydra::get<1>(curr) << std::endl;
            
            t_ec=ce;
            t_hc=ch;
            t_tc=ct;
            t_ec_sec=ce_sec;
            t_hc_sec=ch_sec;
            t_tc_sec=ct_sec;
            
            t_elost=hydra::get<5>(curr);
            t_hlost=hydra::get<4>(curr);
            t_tlost=t_elost+t_hlost;
            t_elost_sec=hydra::get<7>(curr);
            t_hlost_sec=hydra::get<6>(curr);
            t_tlost_sec=t_elost_sec+t_hlost_sec;
            
            ghc->SetPoint(pl,pl*(timestep/UTIME),ch);
            gec->SetPoint(pl,pl*(timestep/UTIME),ce);
            gtc->SetPoint(pl,pl*(timestep/UTIME),ct);
            ghc_sec->SetPoint(pl,pl*(timestep/UTIME),ch_sec);
            gec_sec->SetPoint(pl,pl*(timestep/UTIME),ce_sec);
            gtc_sec->SetPoint(pl,pl*(timestep/UTIME),ct_sec);
            if(ct>currmax) currmax=ct;
            if(ct<currmin) currmin=ct;
            
            t_time=pl*timestep;
            t_einduced+=t_ec*timestep;
            t_hinduced+=t_hc*timestep;
            t_tinduced+=t_tc*timestep;
            
            t_einduced_sec+=t_ec_sec*timestep;
            t_hinduced_sec+=t_hc_sec*timestep;
            t_tinduced_sec+=t_tc_sec*timestep;
            
            if(t_tinduced>chargemax) chargemax=t_tinduced;
            if(t_tinduced<chargemin) chargemin=t_tinduced;
            
            ghind->SetPoint(pl,pl*(timestep/UTIME),t_hinduced);
            geind->SetPoint(pl,pl*(timestep/UTIME),t_einduced);
            gtind->SetPoint(pl,pl*(timestep/UTIME),t_tinduced);
            ghlost->SetPoint(pl,pl*(timestep/UTIME),t_hlost);
            gelost->SetPoint(pl,pl*(timestep/UTIME),t_elost);
            gtlost->SetPoint(pl,pl*(timestep/UTIME),t_tlost);
            
            ghind_sec->SetPoint(pl,pl*(timestep/UTIME),t_hinduced_sec);
            geind_sec->SetPoint(pl,pl*(timestep/UTIME),t_einduced_sec);
            gtind_sec->SetPoint(pl,pl*(timestep/UTIME),t_tinduced_sec);
            ghlost_sec->SetPoint(pl,pl*(timestep/UTIME),t_hlost_sec);
            gelost_sec->SetPoint(pl,pl*(timestep/UTIME),t_elost_sec);
            gtlost_sec->SetPoint(pl,pl*(timestep/UTIME),t_tlost_sec);
            
            pl++;
            tree->Fill();
        }
        
        INFO_LINE("Max Current: "<<currmax<<" A")
        INFO_LINE("Int Charge: "<<t_tinduced<<" C")
        INFO_LINE("Int Charge e: "<<t_einduced<<" C")
        INFO_LINE("Int Charge h: "<<t_hinduced<<" C")
        TH1D *hcurr=new TH1D("h1curr",";time [ps]; Current [A]",1,0,(timestep)/UTIME * tp_currs.size());
        hcurr->GetYaxis()->SetTitleOffset(0.9);
        hcurr->SetMaximum(currmax*1.1);
        hcurr->SetMinimum(currmin);
        
        std::cout << std::endl;
        
        TCanvas *cancurr= new TCanvas("current","Current vs time",1000,600);
        cancurr->SetTopMargin(0.1);
        cancurr->SetLeftMargin(0.15);
        cancurr->SetRightMargin(0.15);
        gPad->SetTicks(1,0);
        hcurr->Draw();
        gtc->Draw("L SAME");
        gec->Draw("L SAME");
        ghc->Draw("L SAME");
        gtc_sec->Draw("L SAME");
        gec_sec->Draw("L SAME");
        ghc_sec->Draw("L SAME");
        TLegend *legc=new TLegend(0.68,0.55,0.83,0.9);
        legc->SetNColumns(2);
        legc->AddEntry(gec,"I_{e}","l");
        legc->AddEntry(geind,"Q^{ind}_{e}","l");
        legc->AddEntry(ghc,"I_{h}","l");
        legc->AddEntry(ghind,"Q^{ind}_{h}","l");
        legc->AddEntry(gtc,"I_{e+h}","l");
        legc->AddEntry(gtind,"Q^{ind}_{e+h}","l");
        legc->AddEntry(gec_sec,"I_{e,sec}","l");
        legc->AddEntry(geind_sec,"Q^{ind}_{e}","l");
        legc->AddEntry(ghc_sec,"I_{h,sec}","l");
        legc->AddEntry(ghind_sec,"Q^{ind}_{h}","l");
        legc->AddEntry(gtc_sec,"I_{e+h,sec}","l");
        legc->AddEntry(gtind_sec,"Q^{ind}_{e+h}","l");
        legc->Draw();
        
        cancurr->cd();
        TPad *overlay = new TPad("overlay","",0.,0.,1.,1.);
        overlay->SetTicks(0,0);
        overlay->SetTopMargin(0.1);
        overlay->SetLeftMargin(0.15);
        overlay->SetRightMargin(0.15);
        overlay->SetFillStyle(4000);
        overlay->SetFillColor(0);
        overlay->SetFrameFillStyle(4000);
        
        overlay->Draw();
        overlay->cd();
        TH1D *hcurroverlay=new TH1D("h1curroverlay"," ; ; Integrated charge [C]",1,0,(timestep)/UTIME * tp_currs.size());
        hcurroverlay->SetMaximum(chargemax*2.1);
        hcurroverlay->SetMinimum(chargemin);
    //     hcurroverlay->GetXaxis()->SetAxisColor(0);
        hcurroverlay->GetYaxis()->SetAxisColor(2);
        hcurroverlay->GetYaxis()->SetTitleColor(2);
        hcurroverlay->GetYaxis()->SetTitleOffset(1.1);
    //     hcurroverlay->GetYaxis()->SetTitle("induced charge [C]");
        hcurroverlay->GetYaxis()->SetLabelColor(2);
        hcurroverlay->Draw("Y+");
        gtind->Draw("L SAME");
        geind->Draw("L SAME");
        ghind->Draw("L SAME");
        gtind_sec->Draw("L SAME");
        geind_sec->Draw("L SAME");
        ghind_sec->Draw("L SAME");
        
        TPaveText* settingstext = new TPaveText(0.55, 0.6, 0.68,0.87, "BRNDC");
        for(auto s:settings){
            if(s.first=="name" || s.first=="outputdir") continue;
            settingstext->SetFillColor(0);
            settingstext->SetTextAlign(12);
            settingstext->SetBorderSize(0);
    //         settingstext->SetTextSize(timetext->GetTextSize()*0.6);
            settingstext->AddText(Form("%s = %s",s.first.c_str(),s.second.c_str()));
            
        }
        settingstext->Draw();
    //     cancurr->SaveAs(Form("%s/%s_current.pdf",outputdir.Data(),filename.Data()),"pdf");
        
        outf->WriteTObject(gtc,gtc->GetName(),"Overwrite");
        outf->WriteTObject(gec,gec->GetName(),"Overwrite");
        outf->WriteTObject(ghc,ghc->GetName(),"Overwrite");
        outf->WriteTObject(gtind,gtind->GetName(),"Overwrite");
        outf->WriteTObject(geind,geind->GetName(),"Overwrite");
        outf->WriteTObject(ghind,ghind->GetName(),"Overwrite");
        outf->WriteTObject(gtlost,gtlost->GetName(),"Overwrite");
        outf->WriteTObject(gelost,gelost->GetName(),"Overwrite");
        outf->WriteTObject(ghlost,ghlost->GetName(),"Overwrite");
        outf->WriteTObject(gtc_sec,gtc_sec->GetName(),"Overwrite");
        outf->WriteTObject(gec_sec,gec_sec->GetName(),"Overwrite");
        outf->WriteTObject(ghc_sec,ghc_sec->GetName(),"Overwrite");
        outf->WriteTObject(gtind_sec,gtind_sec->GetName(),"Overwrite");
        outf->WriteTObject(geind_sec,geind_sec->GetName(),"Overwrite");
        outf->WriteTObject(ghind_sec,ghind_sec->GetName(),"Overwrite");
        outf->WriteTObject(gtlost_sec,gtlost_sec->GetName(),"Overwrite");
        outf->WriteTObject(gelost_sec,gelost_sec->GetName(),"Overwrite");
        outf->WriteTObject(ghlost_sec,ghlost_sec->GetName(),"Overwrite");
        outf->WriteTObject(hcurr,"currframe","Overwrite");
        outf->WriteTObject(cancurr,"current_plots","Overwrite");
        outf->WriteTObject(tree,tree->GetName(),"Overwrite");
        
        //write settings
        for(auto s:settings){
            TNamed *n=new TNamed(s.first.c_str(),s.second.c_str());
            outf->WriteTObject(n,n->GetName(),"Overwrite");
        }
        
        outf->Close();
        
        delete gtc, gec, ghc;
        delete gtind, geind, ghind;
        delete gtlost, gelost, ghlost;
        delete legc;
        delete cancurr;
    }
    
    void ExtraPlots(UniverseHost_t states, std::vector<ReducedTuple_t> &tp_currs, TH3D* hdraw, double timestep, std::map<std::string,std::string> settings, bool storeextra=true, bool drawgif=false){
        
        if(drawgif){
            WARNING_LINE("Animated GIF requested, it may take a long time...")
            std::remove(Form("%s/%s.gif",settings["outputdir"].c_str(),settings["name"].c_str()));
        }
        
        if(storeextra){
            WARNING_LINE("Full info is being saved, the output file may be large!")
            std::remove(Form("%s/%s.gif",settings["outputdir"].c_str(),settings["name"].c_str()));
        }
        
        TFile *outf=new TFile(Form("%s/%s.root",settings["outputdir"].c_str(),settings["name"].c_str()),"UPDATE");
        double t_ec=0, t_hc=0, t_tc=0, t_einduced=0, t_hinduced=0, t_tinduced=0;
        double t_ec_sec=0, t_hc_sec=0, t_tc_sec=0, t_einduced_sec=0, t_hinduced_sec=0, t_tinduced_sec=0;
        
        double currmax=-999999;
        double chargemax=-999999;
        double currmin=-0;
        double chargemin=-0;
        unsigned pl=0;
        
        for(auto& curr:tp_currs){
            
            double ce=hydra::get<1>(curr);
            double ch=hydra::get<0>(curr);
            double ct=ce+ch;
            double ce_sec=hydra::get<3>(curr);
            double ch_sec=hydra::get<2>(curr);
            double ct_sec=ce_sec+ch_sec;
            
            t_ec=ce;
            t_hc=ch;
            t_tc=ct;
            t_ec_sec=ce_sec;
            t_hc_sec=ch_sec;
            t_tc_sec=ct_sec;
            
            if(ct>currmax) currmax=ct;
            if(ct<currmin) currmin=ct;
            
            t_einduced+=t_ec*timestep;
            t_hinduced+=t_hc*timestep;
            t_tinduced+=t_tc*timestep;
            
            t_einduced_sec+=t_ec_sec*timestep;
            t_hinduced_sec+=t_hc_sec*timestep;
            t_tinduced_sec+=t_tc_sec*timestep;
            
            if(t_tinduced>chargemax) chargemax=t_tinduced;
            if(t_tinduced<chargemin) chargemin=t_tinduced;
            
            pl++;
        }
        
        INFO_LINE("Max Current: "<<currmax<<" A")
        TH1D *hcurr=new TH1D("h1curr",";time [ps]; Current [A]",1,0,(timestep)/UTIME * states.size());
        hcurr->GetYaxis()->SetTitleOffset(0.9);
        hcurr->SetMaximum(currmax*1.1);
        hcurr->SetMinimum(currmin);
        
        if(storeextra || drawgif){
            double maxhist=hdraw->GetMaximum();
            double minhist=hdraw->GetMinimum();
            INFO_LINE("Max E: "<<maxhist<<" V/micron")
            INFO_LINE("Min E: "<<minhist<<" V/micron")
            std::string extrafolder="extra_info";
            TDirectory *extradir=outf->mkdir(extrafolder.c_str());
            double sizepadcurrent=0.4;
            
            TH2D *pxy=(TH2D*)hdraw->Project3D("yx");    
            TH2D *pyz=(TH2D*)hdraw->Project3D("zy");
            TH2D *pxz=(TH2D*)hdraw->Project3D("zx");
            
            pxy->Scale(1./hdraw->GetZaxis()->GetNbins());
            pxy->SetMinimum(minhist);
            pyz->Scale(1./hdraw->GetXaxis()->GetNbins());
            pyz->SetMinimum(minhist);
            pxz->Scale(1./hdraw->GetYaxis()->GetNbins());
            pxz->SetMinimum(minhist);
            
            
            double xdim=pxy->GetXaxis()->GetXmax()-pxy->GetXaxis()->GetXmin();
            double ydim=pxy->GetYaxis()->GetXmax()-pxy->GetYaxis()->GetXmin();
            double zdim=pxz->GetYaxis()->GetXmax()-pxz->GetYaxis()->GetXmin();
            
            TCanvas *can= new TCanvas("can","can",2*(xdim+ydim),2*(ydim+zdim)*(1./(1.-sizepadcurrent)));
            

            double vfz=zdim/(ydim+zdim);
            double ofx=xdim/(xdim+ydim);

            TPad *c1=new TPad("c1","xy",0., sizepadcurrent + vfz*(1.-sizepadcurrent) , ofx , 1.);
            c1->SetBorderSize(0);
            c1->SetLeftMargin(0.15);
            c1->SetRightMargin(0);
            c1->SetTopMargin(0.1);
            c1->SetBottomMargin(0);
            c1->Draw();
            can->cd();
            TPad *c2=new TPad("c2","current",ofx, sizepadcurrent + vfz*(1.-sizepadcurrent), 1., 1.);
            c2->SetBorderSize(0);
            c2->SetLeftMargin(0);
            c2->SetRightMargin(0.05);
            c2->SetTopMargin(0.1);
            c2->SetBottomMargin(0);
            c2->Draw();
            can->cd();
            TPad *c3=new TPad("c3","xz",0., sizepadcurrent, ofx, sizepadcurrent+vfz*(1-sizepadcurrent));
            c3->SetBorderSize(0);
            c3->SetLeftMargin(0.15);
            c3->SetRightMargin(0);
            c3->SetTopMargin(0);
            c3->SetBottomMargin(0.15);
            c3->Draw();
            can->cd();
            TPad *c4=new TPad("c4","yz",ofx, sizepadcurrent, 1., sizepadcurrent+vfz*(1-sizepadcurrent));
            c4->SetBorderSize(0);
            c4->SetLeftMargin(0);
            c4->SetRightMargin(0.05);
            c4->SetTopMargin(0);
            c4->SetBottomMargin(0.15);
            c4->Draw();
            can->cd();
            
            TPad *c5=new TPad("c5","current",0., 0., 1., sizepadcurrent);
            c5->SetBorderSize(0);
            c5->SetLeftMargin(0.15);
            c5->SetRightMargin(0.05);
            c5->SetTopMargin(0.1);
            c5->SetBottomMargin(0.2);
            c5->Draw();
            can->cd();
            
            
            TMarker *emark=new TMarker(0,0,1);
            emark->SetMarkerStyle(2);
            emark->SetMarkerSize(1);
            emark->SetMarkerColor(COLORELEC);
            TMarker *hmark=new TMarker(0,0,1);
            hmark->SetMarkerStyle(5);
            hmark->SetMarkerSize(1);
            hmark->SetMarkerColor(COLORHOLE);
            
            TLegend *leg=new TLegend(0.,0.5,1.,1.);
            leg->AddEntry(emark,"electrons","p");
            leg->AddEntry(hmark,"holes","p");
            
            
            pxy->GetYaxis()->SetTitleOffset(0.35);
            pxz->GetYaxis()->SetTitleOffset(0.95);
            
            pxy->GetYaxis()->SetLabelSize(pxy->GetYaxis()->GetLabelSize()*(c3->GetWNDC()*c3->GetHNDC())/(c1->GetWNDC()*c1->GetHNDC()));
            pyz->GetXaxis()->SetLabelSize(pyz->GetXaxis()->GetLabelSize()*(c3->GetWNDC()*c3->GetHNDC())/(c4->GetWNDC()*c4->GetHNDC())*0.9);
            pyz->GetXaxis()->SetTitleSize(pyz->GetXaxis()->GetTitleSize()*(c3->GetWNDC()*c3->GetHNDC())/(c4->GetWNDC()*c4->GetHNDC())*0.9);
            pxy->GetYaxis()->SetTitleSize(pxy->GetYaxis()->GetTitleSize()*(c3->GetWNDC()*c3->GetHNDC())/(c1->GetWNDC()*c1->GetHNDC()));
            pyz->GetXaxis()->SetTitleOffset(0.36);
            pyz->GetXaxis()->SetLabelOffset(-0.070);
            
            
            c1->cd();
            pxy->Draw("col");
            c2->cd();
            leg->Draw();
            c3->cd();
            pxz->Draw("col");
            c4->cd();
            pyz->Draw("col");
            
            
            c5->cd();
            hcurr->Draw("");
            
            pl=0;
            for(auto& final_state:states){
                
                StateHost_t holes;
                StateHost_t electrons;
                for(auto t:final_state){
                    if(hydra::get<0>(t)>0) holes.push_back(t);
                    else electrons.push_back(t);
                }

                
                c1->cd();
                TGraph *_gxy = new TGraph(electrons.size(),const_cast<double*>(electrons.column(1).data()),const_cast<double*>(electrons.column(2).data()));
                _gxy->SetName(Form("e_xyproj_step%i",pl));
                _gxy->SetMarkerStyle(1);
                _gxy->SetMarkerSize(1);
                _gxy->SetMarkerColor(COLORELEC);
                if(_gxy->GetN()>0)_gxy->Draw("P SAME");
                
                TGraph *_ghxy = new TGraph(holes.size(),const_cast<double*>(holes.column(1).data()),const_cast<double*>(holes.column(2).data()));
                _ghxy->SetName(Form("h_xyproj_step%i",pl));
                _ghxy->SetMarkerStyle(1);
                _ghxy->SetMarkerSize(1);
                _ghxy->SetMarkerColor(COLORHOLE);
                if(_ghxy->GetN()>0)_ghxy->Draw("P SAME");
                
                c4->cd();
                TGraph *_gyz = new TGraph(electrons.size(),const_cast<double*>(electrons.column(2).data()),const_cast<double*>(electrons.column(3).data()));
                _gyz->SetName(Form("e_yzproj_step%i",pl));
                _gyz->SetMarkerStyle(1);
                _gyz->SetMarkerSize(1);
                _gyz->SetMarkerColor(COLORELEC);
                if(_gyz->GetN()>0) _gyz->Draw("P SAME");
                
                TGraph *_ghyz = new TGraph(holes.size(),const_cast<double*>(holes.column(2).data()),const_cast<double*>(holes.column(3).data()));
                _ghyz->SetName(Form("h_yzproj_step%i",pl));
                _ghyz->SetMarkerStyle(1);
                _ghyz->SetMarkerSize(1);
                _ghyz->SetMarkerColor(COLORHOLE);
                if(_ghyz->GetN()>0)_ghyz->Draw("P SAME");
                
                c3->cd();
                TGraph *_gxz = new TGraph(electrons.size(),const_cast<double*>(electrons.column(1).data()),const_cast<double*>(electrons.column(3).data()));
                _gxz->SetName(Form("e_xzproj_step%i",pl));
                _gxz->SetMarkerStyle(1);
                _gxz->SetMarkerSize(1);
                _gxz->SetMarkerColor(COLORELEC);
                if(_gxz->GetN()>0)_gxz->Draw("P SAME");
                
                
                TGraph *_ghxz = new TGraph(holes.size(),const_cast<double*>(holes.column(1).data()),const_cast<double*>(holes.column(3).data()));
                _ghyz->SetName(Form("h_yzproj_step%i",pl));
                _ghxz->SetMarkerStyle(1);
                _ghxz->SetMarkerSize(1);
                _ghxz->SetMarkerColor(COLORHOLE);
                if(_ghxz->GetN()>0)_ghxz->Draw("P SAME");
                
                
                c5->cd();
                
                double ce=hydra::get<1>(tp_currs[pl]);
                double ch=hydra::get<0>(tp_currs[pl]);
                double ct=ce+ch;
                
                TMarker *_emark=new TMarker(timestep/UTIME*pl,ce,1);
                _emark->SetMarkerStyle(2);
                _emark->SetMarkerSize(0.5);
                _emark->SetMarkerColor(COLORELEC);
                _emark->Draw();
                
                TMarker *_hmark=new TMarker(timestep/UTIME*pl,ch,1);
                _hmark->SetMarkerStyle(5);
                _hmark->SetMarkerSize(0.5);
                _hmark->SetMarkerColor(COLORHOLE);
                _hmark->Draw();
                
                TMarker *_tmark=new TMarker(timestep/UTIME*pl,ct,1);
                _tmark->SetMarkerStyle(31);
                _tmark->SetMarkerSize(0.5);
                _tmark->SetMarkerColor(1);
                _tmark->Draw();
                
                
                c2->cd();
                TPaveText* timetext = new TPaveText(0.1, 0.1, 0.9,0.4, "BRNDC");
                timetext->SetFillColor(0);
                timetext->SetTextAlign(12);
                timetext->SetBorderSize(0);
                timetext->SetTextSize(timetext->GetTextSize()*0.6);
                timetext->AddText(Form("Time [ps]: %3.0f",timestep/UTIME*pl));
                timetext->Draw();
                
                c5->cd();
                TPaveText* settingstext = new TPaveText(0.65, 0.6, 0.85,0.87, "BRNDC");
                for(auto s:settings){
                    if(s.first=="name" || s.first=="outputdir") continue;
                    settingstext->SetFillColor(0);
                    settingstext->SetTextAlign(12);
                    settingstext->SetBorderSize(0);
            //         settingstext->SetTextSize(timetext->GetTextSize()*0.6);
                    settingstext->AddText(Form("%s = %s",s.first.c_str(),s.second.c_str()));
                    
                }

                settingstext->Draw();
                
                if(drawgif) can->Print(Form("%s/%s.gif+1",settings["outputdir"].c_str(),settings["name"].c_str()));
                
                if(storeextra){
                    if(_gxy->GetN()>0)  extradir->WriteTObject(_gxy,_gxy->GetName(),"Overwrite");
                    if(_ghxy->GetN()>0) extradir->WriteTObject(_ghxy,_ghxy->GetName(),"Overwrite");
                    if(_gyz->GetN()>0)  extradir->WriteTObject(_gyz,_gyz->GetName(),"Overwrite");
                    if(_ghyz->GetN()>0) extradir->WriteTObject(_ghyz,_ghyz->GetName(),"Overwrite");
                    if(_gxz->GetN()>0)  extradir->WriteTObject(_gxz,_gxz->GetName(),"Overwrite");
                    if(_ghxz->GetN()>0) extradir->WriteTObject(_ghxz,_ghxz->GetName(),"Overwrite");
                    extradir->WriteTObject(can,Form("Projections_step%i",pl),"Overwrite");
                }
                
                delete settingstext;
                delete timetext;
                delete _gxy;
                delete _gyz;
                delete _gxz;
                delete _ghxy;
                delete _ghyz;
                delete _ghxz;
//                 delete _emark;
//                 delete _hmark;
//                 delete _tmark;
                
                pl++;
                
                printProgress((double)pl/(double)states.size(),"Producing extra plots:");
            }
            
            
            
            
            can->Print(Form("%s/%s.gif++",settings["outputdir"].c_str(),settings["name"].c_str()));
            delete can;
            outf->cd();
        }
        std::cout << std::endl;
        
        
        outf->Close();
        
    }
}
