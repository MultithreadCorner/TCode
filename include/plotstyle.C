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
 * plotstyle.C
 *
 *  Created on: 12/11/2018
 *      Author: Andrea Contu
 */

#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPaveText.h"
void setstyle(){// all users - please change the name of this file to lhcbStyle.C

  // Use times new roman, precision 2 
  Int_t plotFont        = 132;  // Old TIMESPOT style: 62;
  // Line thickness
  Double_t plotWidth    = 2.00; // Old TIMESPOT style: 3.00;
  // Text size
  Double_t plotTSize    = 0.06; 
  
  // use plain black on white colors
  gROOT->SetStyle("Plain"); 
  TStyle *plotStyle= new TStyle("plotStyle","TIMESPOT plots style");
  
  //plotStyle->SetErrorX(0); //  don't suppress the error bar along X

  plotStyle->SetFillColor(1);
  plotStyle->SetFillStyle(1001);   // solid
  plotStyle->SetFrameFillColor(0);
  plotStyle->SetFrameBorderMode(0);
  plotStyle->SetPadBorderMode(0);
  plotStyle->SetPadColor(0);
  plotStyle->SetCanvasBorderMode(0);
  plotStyle->SetCanvasColor(0);
  plotStyle->SetStatColor(0);
  plotStyle->SetLegendBorderSize(0);
  plotStyle->SetLegendFont(132);

  // If you want the usual gradient palette (blue -> red)
  plotStyle->SetPalette(1);
  // If you want colors that correspond to gray scale in black and white:
  int colors[8] = {0,5,7,3,6,2,4,1};
  plotStyle->SetPalette(8,colors);

  // set the paper & margin sizes
  plotStyle->SetPaperSize(20,26);
  plotStyle->SetPadTopMargin(0.05);
  plotStyle->SetPadRightMargin(0.05); // increase for colz plots
  plotStyle->SetPadBottomMargin(0.17);
  plotStyle->SetPadLeftMargin(0.19);
  
  // use large fonts
  plotStyle->SetTextFont(plotFont);
  plotStyle->SetTextSize(plotTSize);
  plotStyle->SetLabelFont(plotFont,"x");
  plotStyle->SetLabelFont(plotFont,"y");
  plotStyle->SetLabelFont(plotFont,"z");
  plotStyle->SetLabelSize(plotTSize,"x");
  plotStyle->SetLabelSize(plotTSize,"y");
  plotStyle->SetLabelSize(plotTSize,"z");
  plotStyle->SetTitleFont(plotFont);
  plotStyle->SetTitleFont(plotFont,"x");
  plotStyle->SetTitleFont(plotFont,"y");
  plotStyle->SetTitleFont(plotFont,"z");
  plotStyle->SetTitleSize(1.2*plotTSize,"x");
  plotStyle->SetTitleSize(1.2*plotTSize,"y");
  plotStyle->SetTitleSize(1.2*plotTSize,"z");

  // use medium bold lines and thick markers
  plotStyle->SetLineWidth(plotWidth);
  plotStyle->SetFrameLineWidth(plotWidth);
  plotStyle->SetHistLineWidth(plotWidth);
  plotStyle->SetFuncWidth(plotWidth);
  plotStyle->SetGridWidth(plotWidth);
  plotStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  plotStyle->SetMarkerStyle(20);
  plotStyle->SetMarkerSize(1.0);

  // label offsets
  plotStyle->SetLabelOffset(0.012,"X");
  plotStyle->SetLabelOffset(0.012,"Y");

  // by default, do not display histogram decorations:
  plotStyle->SetOptStat(0);  
  //plotStyle->SetOptStat("emr");  // show only nent -e , mean - m , rms -r
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  plotStyle->SetStatFormat("6.3g"); // specified as c printf options
  plotStyle->SetOptTitle(0);
  plotStyle->SetOptFit(0);
  //plotStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
  //titles
  plotStyle->SetTitleOffset(0.95,"X");
  plotStyle->SetTitleOffset(0.95,"Y");
  plotStyle->SetTitleOffset(1.2,"Z");
  plotStyle->SetTitleFillColor(0);
  plotStyle->SetTitleStyle(0);
  plotStyle->SetTitleBorderSize(0);
  plotStyle->SetTitleFont(plotFont,"title");
  plotStyle->SetTitleX(0.0);
  plotStyle->SetTitleY(1.0); 
  plotStyle->SetTitleW(1.0);
  plotStyle->SetTitleH(0.05);
  
  // look of the statistics box:
  plotStyle->SetStatBorderSize(0);
  plotStyle->SetStatFont(plotFont);
  plotStyle->SetStatFontSize(0.05);
  plotStyle->SetStatX(0.9);
  plotStyle->SetStatY(0.9);
  plotStyle->SetStatW(0.25);
  plotStyle->SetStatH(0.15);

  // put tick marks on top and RHS of plots
  plotStyle->SetPadTickX(1);
  plotStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  plotStyle->SetNdivisions(505,"x");
  plotStyle->SetNdivisions(510,"y");

    const Int_t NRGBs = 3;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = {0.00, 0.5, 1.00};
    Double_t red[NRGBs]   = {0.99, 0.4, 0.81};
    Double_t green[NRGBs] = {0.99, 0.7, 0.30};
    Double_t blue[NRGBs]  = {0.99, 0.8, 0.30};
    Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    Int_t MyPalette[NCont];
    for (int i=0;i<NCont;i++) MyPalette[i] = FI+i;
    plotStyle->SetPalette(NCont, MyPalette);
    plotStyle->SetTitleOffset(1.5,"Y");
    plotStyle->SetTitleOffset(1.0,"X");
    plotStyle->SetNdivisions(505,"Y");
  
  
  gROOT->SetStyle("plotStyle");
  gROOT->ForceStyle();

  // add project label
  TPaveText* plotName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                                      0.87 - gStyle->GetPadTopMargin(),
                                      gStyle->GetPadLeftMargin() + 0.20,
                                      0.95 - gStyle->GetPadTopMargin(),
                                      "BRNDC");
  plotName->AddText(PROJECT_NAME);
  plotName->SetFillColor(0);
  plotName->SetTextAlign(12);
  plotName->SetBorderSize(0);

  TText *plotLabel = new TText();
  plotLabel->SetTextFont(plotFont);
  plotLabel->SetTextColor(1);
  plotLabel->SetTextSize(plotTSize);
  plotLabel->SetTextAlign(12);

  TLatex *plotLatex = new TLatex();
  plotLatex->SetTextFont(plotFont);
  plotLatex->SetTextColor(1);
  plotLatex->SetTextSize(plotTSize);
  plotLatex->SetTextAlign(12);
  
  

//   std::cout << "-------------------------" << std::endl;  
//   std::cout << "Set LHCb Style - Feb 2012" << std::endl;
//   std::cout << "-------------------------" << std::endl;  
  
}


