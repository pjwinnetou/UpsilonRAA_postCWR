#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"
#include "../commonUtility.h"

void strickland_RAA_cent_newglobal_all_asym_postCWRConstrain(bool isArrow =true)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 100; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nState = 3; // Y(1S), Y(2S), and Y(3S)
  double xmax = 420.0;
  double xmin_int = 0.6;
  double xmax_int = 1.8;
  double boxw = 6.5; // for syst. box (vs cent)
  double boxw_int = 0.09;
//  double relsys = 0.1;


  ///////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn[nState];
	TGraphErrors* gRAA[nState]; // vs centrality
  TGraphAsymmErrors* gRAA_sys[nState];
	TGraphErrors* gRAA_int[nState]; // centrality-integrated
  TGraphAsymmErrors* gRAA_int_sys[nState];
  for (int is=0; is<nState; is++){
  	fIn[is] = new TFile(Form("Ups_%d_RAA.root",is+1),"READ");
    gRAA[is]=(TGraphErrors*)fIn[is]->Get("gRAA_cent");
    //gRAA_sys[is]=(TGraphErrors*)fIn[is]->Get("gRAA_cent");
    gRAA_sys[is]= new TGraphAsymmErrors();
    gRAA_int[is]=(TGraphErrors*)fIn[is]->Get("gRAA_int");
    //gRAA_int_sys[is]=(TGraphAsymmErrors*)fIn[is]->Get("gRAA_int");
    gRAA_int_sys[is]= new TGraphAsymmErrors();
    //cout << "gRAA["<<is<<"] = " <<gRAA[is] << endl;
  }
  //// read input file : syst.
  TFile* fInSys_Hi[nState];
  TFile* fInSys_Lo[nState];
  TH1D* hSys_Hi[nState];
  TH1D* hSys_Lo[nState];
  TH1D* hSys_int_Hi[nState];
  TH1D* hSys_int_Lo[nState];
  int npoint_Hi[nState];
  int npoint_Lo[nState];
  int npoint_int_Hi[nState];
  int npoint_int_Lo[nState];
  for (int is=0; is<nState; is++){
    fInSys_Hi[is] = new TFile(Form("../Systematic_new/postCWR_mergedSys_constrain_ups%ds_asymHi.root",is+1),"READ");
    fInSys_Lo[is] = new TFile(Form("../Systematic_new/postCWR_mergedSys_constrain_ups%ds_asymLo.root",is+1),"READ");
    hSys_Hi[is]=(TH1D*)fInSys_Hi[is]->Get("hcentRAA_merged");
    hSys_Lo[is]=(TH1D*)fInSys_Lo[is]->Get("hcentRAA_merged");
    npoint_Hi[is] = hSys_Hi[is]->GetSize()-2;
    npoint_Lo[is] = hSys_Lo[is]->GetSize()-2;
    cout << "*** CENT Hi *** Y("<<is+1<<") : # of point = " << npoint_Hi[is] << endl;
    cout << "*** CENT Lo *** Y("<<is+1<<") : # of point = " << npoint_Lo[is] << endl;
    hSys_int_Hi[is]=(TH1D*)fInSys_Hi[is]->Get("hintRAA_merged");
    hSys_int_Lo[is]=(TH1D*)fInSys_Lo[is]->Get("hintRAA_merged");
    npoint_int_Hi[is] = hSys_int_Hi[is]->GetSize()-2;
    npoint_int_Lo[is] = hSys_int_Lo[is]->GetSize()-2;
    cout << "*** INT Hi *** Y("<<is+1<<") : # of point = " << npoint_int_Hi[is] << endl;
    cout << "*** INT Lo *** Y("<<is+1<<") : # of point = " << npoint_int_Lo[is] << endl;
  }   
  
  //// set bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys_Lo, relsys_Hi;
  double xshift = 0.2;
  //// --- vs centrality
  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
    if (npoint_Lo[is] != gRAA[is]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
    for (int ipt=0; ipt< npoint_Lo[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys_Lo=0; relsys_Hi=0;
      gRAA[is]->GetPoint(ipt, pxtmp, pytmp);
      extmp=gRAA[is]->GetErrorX(ipt);
      eytmp=gRAA[is]->GetErrorY(ipt);
      relsys_Hi=hSys_Hi[is]->GetBinContent(npoint_Hi[is]-ipt);
      relsys_Lo=hSys_Lo[is]->GetBinContent(npoint_Lo[is]-ipt);
      cout << ipt <<"th bin RAA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      //cout << ipt <<"th bin rel. syst. = " << relsys << endl;
      cout << ipt <<"th bin Hi syst. = " << pytmp*relsys_Hi << endl; 
      cout << ipt <<"th bin Lo syst. = " << pytmp*relsys_Lo << endl; 
      //// 1) remove ex from gRAA
      gRAA[is]->SetPointError(ipt, 0, eytmp);
      //// 2) set ey for gRAA_sys 
      //gRAA_sys[is]->SetPointError(ipt, extmp, pytmp*relsys);
      //cout << "asdasd" << endl;
      gRAA_sys[is]->SetPoint(ipt, pxtmp,pytmp);
      cout << "asdasd" << endl;
      gRAA_sys[is]->SetPointError(ipt, boxw, boxw, pytmp*relsys_Lo,pytmp*relsys_Hi);
    }
  }

  //// --- centrality-integrated
  cout << " " << endl;
  cout << " INTEGRATED" << endl;
  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
  cout << "npoint_int_Lo " << npoint_int_Lo[is] << endl;
    if (npoint_int_Lo[is] != gRAA_int[is]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }    
    for (int ipt=0; ipt< npoint_int_Lo[is] ; ipt++) {
      pxtmp=0; pytmp=0; extmp=0; eytmp=0;
      gRAA_int[is]->GetPoint(ipt, pxtmp, pytmp);
      extmp=gRAA_int[is]->GetErrorX(ipt);
      eytmp=gRAA_int[is]->GetErrorY(ipt);
      relsys_Hi = hSys_int_Hi[is]->GetBinContent(ipt+1);
      relsys_Hi = TMath::Sqrt(relsys_Hi*relsys_Hi+lumi_unc_pp*lumi_unc_pp + nMB_unc*nMB_unc);
      relsys_Lo = hSys_int_Lo[is]->GetBinContent(ipt+1);
      relsys_Lo = TMath::Sqrt(relsys_Lo*relsys_Lo+lumi_unc_pp*lumi_unc_pp + nMB_unc*nMB_unc);
      cout << ipt <<"th bin RAA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      //cout << ipt <<"th bin rel. syst. = " << relsys << endl;
      cout << ipt <<"th bin Hi syst. = " << pytmp*relsys_Hi << endl; 
      cout << ipt <<"th bin Lo syst. = " << pytmp*relsys_Lo << endl; 
      //// 1) remove ex from gRAA
      gRAA_int[is]->SetPoint(ipt, pxtmp+xshift*is, pytmp);
      gRAA_int[is]->SetPointError(ipt, 0, eytmp);
      //// 2) set ey for gRAA_sys
      gRAA_int_sys[is]->SetPoint(ipt, pxtmp+xshift*is, pytmp);
      //gRAA_int_sys[is]->SetPointError(ipt, extmp, pytmp*relsys);
      gRAA_int_sys[is]->SetPointError(ipt, boxw_int,boxw_int, pytmp*relsys_Lo,pytmp*relsys_Hi); //extemp fixed
    }
  }

  ////////////////////////////////////////////////////////////////
  //// 3S upper limit (arrow)
  int ulstate = 2; //3S
  static const int n3s = 2;
  double lower68[n3s] = {lower68_c1,lower68_c2};
  double upper68[n3s] = {upper68_c1,upper68_c2};
  double lower95[n3s] = {lower95_c1,lower95_c2};
  double upper95[n3s] = {upper95_c1,upper95_c2};
  static const int n3s_int = 1;
  double lower68_int[n3s_int] = {lower68_cint};
  double upper68_int[n3s_int] = {upper68_cint};
  double lower95_int[n3s_int] = {lower95_cint};
  double upper95_int[n3s_int] = {upper95_cint};

  if (n3s != npoint_Lo[ulstate]) {cout<<"ERROR!! # of bins for UL is wrong!!"<<endl;return;} 
  if (n3s_int != npoint_int_Lo[ulstate]) {cout<<"ERROR!! # of bins for UL (int) is wrong!!"<<endl;return;} 

  //// --- vs centrality
  TBox *box68per[n3s];
  TArrow *arr95per[n3s];
  for (int ipt=0; ipt< n3s ; ipt++) { //bin by bin
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; 
    //lower68=0; upper68=0; lower95=0; upper95=0; 
    gRAA[ulstate]->GetPoint(ipt, pxtmp, pytmp);
    box68per[ipt] = new TBox(pxtmp-boxw,lower68[ipt],pxtmp+boxw,upper68[ipt]);
    arr95per[ipt] = new TArrow(pxtmp,lower95[ipt],pxtmp,upper95[ipt],0.027,"<-|"); //95%    box68per[ipt]->SetLineColor(kGreen+2);
    box68per[ipt]->SetFillColorAlpha(kGreen-10,0.5);
    box68per[ipt]->SetLineWidth(1);
    arr95per[ipt]->SetLineColor(kGreen+2);
    arr95per[ipt]->SetLineWidth(2);
  }
  //// --- centrality-integrated
  TBox *box68per_int[n3s_int];
  TArrow *arr95per_int[n3s_int];
  for (int ipt=0; ipt< n3s_int ; ipt++) { //bin by bin
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; 
    gRAA_int[ulstate]->GetPoint(ipt, pxtmp, pytmp);
    box68per_int[ipt] = new TBox(pxtmp-boxw_int,lower68_int[ipt],pxtmp+boxw_int,upper68_int[ipt]);
    arr95per_int[ipt] = new TArrow(pxtmp,lower95_int[ipt],pxtmp,upper95_int[ipt],0.02,"<-|"); //95%
    box68per_int[ipt]->SetLineColor(kGreen+2);
    box68per_int[ipt]->SetFillColorAlpha(kGreen-10,0.5);
    box68per_int[ipt]->SetLineWidth(1);
    arr95per_int[ipt]->SetLineColor(kGreen+2);
    arr95per_int[ipt]->SetLineWidth(2);
  }
  
  double lower68_2s = 0;
  double lower95_2s = 0;
  double upper68_2s = upper68_2s_c1;
  double upper95_2s = upper95_2s_c1;

  TBox *box68per_2s;
  TArrow *arr95per_2s;
  pxtmp=0; pytmp=0; extmp=0; eytmp=0; 
  gRAA[1]->GetPoint(8,pxtmp,pytmp);
  box68per_2s = new TBox(pxtmp-boxw,lower68_2s,pxtmp+boxw,upper68_2s);
  arr95per_2s = new TArrow(pxtmp,lower95_2s,pxtmp,upper95_2s,0.027,"<-|");

  gRAA[1]->SetPoint(8,-10,-10);
  gRAA[1]->SetPointError(8,0,0);
  gRAA_sys[1]->SetPoint(8,-10,-10);
  gRAA_sys[1]->SetPointError(8,0,0,0,0);

  ////////////////////////////////////////////////////////////////
  
  //// graph style 
  for (int is=0; is<nState; is++){
    SetGraphStyle(gRAA[is], is, is); 
    SetGraphStyleSys(gRAA_sys[is], is); 
    SetGraphStyle(gRAA_int[is], is, is); 
    SetGraphStyleSys(gRAA_int_sys[is], is); 
	}
  
  double xlonger = 120;
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.0387);
  
  //// axis et. al
  gRAA_sys[0]->GetXaxis()->SetTitle("<N_{part}>");
  gRAA_sys[0]->GetXaxis()->CenterTitle();
  gRAA_sys[0]->GetYaxis()->SetTitle("R_{AA}");
  gRAA_sys[0]->GetYaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetLimits(0.,xmax);
  gRAA_sys[0]->SetMinimum(0.0);
  gRAA_sys[0]->SetMaximum(1.3);
  //// for cent
  gRAA_sys[0]->GetXaxis()->SetTitleSize(0.06*1.0);
  gRAA_sys[0]->GetYaxis()->SetTitleSize(0.06*1.0);
  gRAA_sys[0]->GetXaxis()->SetLabelSize(0.05*1.0);
  gRAA_sys[0]->GetYaxis()->SetLabelSize(0.05*1.0);
  
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600+xlonger,600);
  TPad* pad_diff = new TPad("pad_diff", "",0, 0, 600/(600.+xlonger), 1.0); // vs centrality
  pad_diff->SetRightMargin(0);
  pad_diff->SetBottomMargin(0.14);
  pad_diff->SetTopMargin(0.067);
  TPad* pad_int = new TPad("pad_int", "",600/(600.+xlonger), 0, 1.0, 1.0); // centrality-integrated
  pad_int->SetBottomMargin(0.14);
  pad_int->SetTopMargin(0.067);
  pad_int->SetLeftMargin(0);
  pad_int->SetRightMargin(0.032*600/xlonger);

  //// --- 1st pad!!!   
  c1->cd();
  pad_diff->Draw(); 
  pad_diff->cd(); 
  //// syst
  for (int is=0; is<nState; is++){
    if ( is==0) { gRAA_sys[is]->Draw("A5"); }
    else if (is==ulstate && isArrow==true) { 
      for (int ipt=0; ipt< n3s ; ipt++) { //bin by bin
        box68per[ipt]->Draw("l"); 
      }
    }
    else { gRAA_sys[is]->Draw("5"); }
	}
  //// point
  for (int is=0; is<nState; is++){
    if (is==ulstate && isArrow==true) {
      for (int ipt=0; ipt< n3s ; ipt++) { //bin by bin
        arr95per[ipt]->Draw();
      }
    }
    else { gRAA[is]->Draw("P"); }
	}
  box68per_2s->SetLineColor(kBlue-3);
//  box68per_2s->SetFillColorAlpha(kBlue-10,0.5);
  box68per_2s->SetFillColor(kBlue-3);
  box68per_2s->SetFillStyle(3005);
  box68per_2s->SetLineWidth(1);
  arr95per_2s->SetLineColor(kBlue-3);
  arr95per_2s->SetLineWidth(2);
  box68per_2s->Draw("lf");
  arr95per_2s->Draw();
  dashedLine(0.,1.,xmax,1.,1,1);
  
  //// legend
  //TLegend *leg= new TLegend(0.75, 0.50, 0.95, 0.70);
  TLegend *leg= new TLegend(0.39, 0.59, 0.59, 0.7);
  SetLegendStyle(leg);
  leg->SetTextSize(0.0385);
  
  TLegend *leg_upp_68= new TLegend(0.57, 0.59, 0.74, 0.76);
  SetLegendStyle(leg_upp_68);
  leg_upp_68->SetTextSize(0.0386);
  TLegend *leg_upp_95= new TLegend(0.74, 0.59, 0.91, 0.76);
  SetLegendStyle(leg_upp_95);
  leg_upp_95->SetTextSize(0.0386);

  TArrow *arrLeg_2s = new TArrow(300.,0.844,300.,0.898,0.022,"<-|");
  arrLeg_2s->SetLineColor(kBlue-3);
  arrLeg_2s->SetLineWidth(2);
  TArrow *arrLeg_3s = new TArrow(300.,0.749,300.,0.802,0.022,"<-|");
  arrLeg_3s->SetLineColor(kGreen+2);
  arrLeg_3s->SetLineWidth(2);

  leg_upp_68->SetHeader("68\% CL","");
  leg_upp_95->SetHeader("95\% CL","");
  
  if (isArrow==false) { 
    for (int is=0; is<nState; is++){
      leg -> AddEntry(gRAA[is],Form(" #varUpsilon(%dS)",is+1),"lp");
    }
  }
  else {
    leg -> AddEntry(gRAA[0]," #varUpsilon(1S)","lp");
    leg -> AddEntry(gRAA[1]," #varUpsilon(2S)","lp");
    TLegendEntry *ent2s=leg_upp_68->AddEntry("ent"," #varUpsilon(2S)","f");
    ent2s->SetLineColor(kBlue-3);
    ent2s->SetFillColor(kBlue-3);
//    ent->SetFillColorAlpha(kGreen-10,0.5);
    ent2s->SetFillStyle(3005);
    
    TLegendEntry *ent=leg_upp_68->AddEntry("ent"," #varUpsilon(3S)","f");
    ent->SetLineColor(kGreen+2);
    ent->SetFillColor(kGreen+2);
 //   ent->SetFillColorAlpha(kGreen-10,0.5);
    ent->SetFillStyle(3005);
    ent2s=leg_upp_95->AddEntry("ent"," #varUpsilon(2S)","");
    ent=leg_upp_95->AddEntry("ent"," #varUpsilon(3S)","");
  }
  leg->Draw("same");
  leg_upp_68->Draw("same");
  leg_upp_95->Draw("same");
  if(isArrow==true) {arrLeg_3s->Draw();arrLeg_2s->Draw();}
  
  //// draw text
  double sz_init = 0.874; double sz_step = 0.0558;
//  globtex->DrawLatex(0.22+0.04, sz_init, "p_{T}^{#mu} > 4 GeV/c");
  globtex->DrawLatex(0.22+0.1, sz_init, "p_{T} < 30 GeV");
//  globtex->DrawLatex(0.46+0.04, sz_init+0.002, "|#eta|^{#mu} < 2.4");
  globtex->DrawLatex(0.22+0.1, sz_init-sz_step, "|y| < 2.4");
/*
  TLatex* centtex = new TLatex();
  centtex->SetNDC();
  centtex->SetTextAlign(12); //left-center
  centtex->SetTextFont(42);
  centtex->SetTextSize(0.029);

  centtex->DrawLatex(0.908,0.37,"0-5%");
  centtex->DrawLatex(0.802,0.37,"5-10%");
  centtex->DrawLatex(0.666,0.37,"10-20%");
  centtex->DrawLatex(0.518,0.421,"20-30%");
  centtex->DrawLatex(0.403,0.465,"30-40%");
  centtex->DrawLatex(0.318,0.495,"40-50%");
  centtex->DrawLatex(0.258,0.555,"50-60%");
  centtex->DrawLatex(0.242,0.698,"60-70%");
  centtex->DrawLatex(0.181,0.781,"70-100%");
*/

//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Centrality 0-100%");

  TFile *fstrickland = new TFile("TheoryCurve/StrickLand_RAA_5023.root","READ");
  
  TGraphErrors *gRAA_1S_strickland[3]; 
  TGraphErrors *gRAA_2S_strickland[3]; 
  TGraphErrors *gRAA_3S_strickland[3]; 
  
  for(int i=0;i<3;i++)
  {
    gRAA_1S_strickland[i] = (TGraphErrors*) fstrickland-> Get(Form("RAA_strick_nPart_1S_%d",i));
    gRAA_2S_strickland[i] = (TGraphErrors*) fstrickland-> Get(Form("RAA_strick_nPart_2S_%d",i));
    gRAA_3S_strickland[i] = (TGraphErrors*) fstrickland-> Get(Form("RAA_strick_nPart_3S_%d",i));
    gRAA_1S_strickland[i] -> SetLineWidth(3.);
    gRAA_2S_strickland[i] -> SetLineWidth(3.0);
    gRAA_3S_strickland[i] -> SetLineWidth(3);
  }
  gRAA_1S_strickland[0]->SetLineColor(kRed+3);
  gRAA_1S_strickland[1]->SetLineColor(kRed+3);
  gRAA_1S_strickland[2]->SetLineColor(kRed+3);
  gRAA_1S_strickland[0]->SetLineStyle(3);
  gRAA_1S_strickland[1]->SetLineStyle(1);
  gRAA_1S_strickland[2]->SetLineStyle(8);
  
  gRAA_2S_strickland[0]->SetLineColor(kBlue+3);
  gRAA_2S_strickland[1]->SetLineColor(kBlue+3);
  gRAA_2S_strickland[2]->SetLineColor(kBlue+3);
  gRAA_2S_strickland[0]->SetLineStyle(3);
  gRAA_2S_strickland[1]->SetLineStyle(1);
  gRAA_2S_strickland[2]->SetLineStyle(8);
  
  gRAA_3S_strickland[0]->SetLineColor(kGreen+2);
  gRAA_3S_strickland[1]->SetLineColor(kGreen+2);
  gRAA_3S_strickland[2]->SetLineColor(kGreen+2);
  gRAA_3S_strickland[0]->SetLineStyle(3);
  gRAA_3S_strickland[1]->SetLineStyle(1);
  gRAA_3S_strickland[2]->SetLineStyle(8);
  

  for(int i=0;i<3;i++){
    gRAA_1S_strickland[i]->Draw("same");
    gRAA_2S_strickland[i]->Draw("same");
    gRAA_3S_strickland[i]->Draw("same");
  }
   
  TLegend *leg_strick= new TLegend(0.29, 0.586, 0.46, 0.716);
  SetLegendStyle(leg_strick);
  leg_strick->SetTextSize(0.040);
  leg_strick->AddEntry(gRAA_1S_strickland[2],"Y(1S)","l");
  leg_strick->AddEntry(gRAA_2S_strickland[2],"Y(2S)","l");
//  leg_strick->Draw("same");

  double line_y = 0.64;
  double line_y_diff = 0.07;
  double line_y_diff_in = 0.02;
  double line_x_end = 277;//122
  double line_x_start = 258;//97
  TLine* t1 = new TLine(line_x_start,line_y,line_x_end,line_y);
  t1->SetLineStyle(3);
  t1->SetLineWidth(2);
  t1->SetLineColor(kRed+3);
  t1->Draw("same");

  TLine* t11 = new TLine(line_x_start,line_y-line_y_diff_in,line_x_end,line_y-line_y_diff_in);
  t11->SetLineStyle(3);
  t11->SetLineWidth(2);
  t11->SetLineColor(kBlue-3);
  t11->Draw("same");

  TLine* t111 = new TLine(line_x_start,line_y-line_y_diff_in*2,line_x_end,line_y-line_y_diff_in*2);
  t111->SetLineStyle(3);
  t111->SetLineWidth(2);
  t111->SetLineColor(kGreen+2);
  t111->Draw("same");

  TLine* t2 = new TLine(line_x_start,line_y-line_y_diff-line_y_diff_in,line_x_end,line_y-line_y_diff-line_y_diff_in);
  t2->SetLineStyle(1);
  t2->SetLineWidth(2);
  t2->SetLineColor(kRed+3);
  t2->Draw("same");

  TLine* t22 = new TLine(line_x_start,line_y-line_y_diff-line_y_diff_in*2,line_x_end,line_y-line_y_diff-line_y_diff_in*2);
  t22->SetLineStyle(1);
  t22->SetLineWidth(2);
  t22->SetLineColor(kBlue-3);
  t22->Draw("same");

  TLine* t222 = new TLine(line_x_start,line_y-line_y_diff-line_y_diff_in*3,line_x_end,line_y-line_y_diff-line_y_diff_in*3);
  t222->SetLineStyle(1);
  t222->SetLineWidth(2);
  t222->SetLineColor(kGreen+2);
  t222->Draw("same");

  TLine* t3 = new TLine(line_x_start,line_y-line_y_diff*2-line_y_diff_in*2,line_x_end,line_y-line_y_diff*2-line_y_diff_in*2);
  t3->SetLineStyle(8);
  t3->SetLineWidth(2);
  t3->SetLineColor(kRed+3);
  t3->Draw("same");

  TLine* t33 = new TLine(line_x_start,line_y-line_y_diff*2-line_y_diff_in*3,line_x_end,line_y-line_y_diff*2-line_y_diff_in*3);
  t33->SetLineStyle(8);
  t33->SetLineWidth(2);
  t33->SetLineColor(kBlue-3);
  t33->Draw("same");

  TLine* t333 = new TLine(line_x_start,line_y-line_y_diff*2-line_y_diff_in*4,line_x_end,line_y-line_y_diff*2-line_y_diff_in*4);
  t333->SetLineStyle(8);
  t333->SetLineWidth(2);
  t333->SetLineColor(kGreen+2);
  t333->Draw("same");


  double ThTextSize = 21;
  drawText2("4#pi #eta/s=1", line_x_end+7, line_y-0.038, ThTextSize);
  drawText2("4#pi #eta/s=2", line_x_end+7, line_y-line_y_diff*1-0.038 - line_y_diff_in, ThTextSize);
  drawText2("4#pi #eta/s=3", line_x_end+7, line_y-line_y_diff*2-0.038 - line_y_diff_in*2, ThTextSize);

  drawText2("Krouppa, Strickland",line_x_start,line_y+0.034,ThTextSize);

  //Global Unc.
  TFile* fppInt = new TFile("../fitResults/Constrain/fitresults_upsilon_fixParm1_seed2_DoubleCB_PP_DATA_pt0.0-30.0_y0.0-2.4_muPt4.0.root");
  TH1D* hppInt = (TH1D*) fppInt -> Get("fitResults");
  double relsys_ppInt1S = hppInt->GetBinError(1)/hppInt->GetBinContent(1);
  double relsys_ppInt2S = hppInt->GetBinError(2)/hppInt->GetBinContent(2);

  TH1D* hSys_glb[nState];
  double sys_global_pp[nState];
  double sys_global_val;
  double accept_sys;
  TFile* f_acc[nState]; 
  TH1D* hSys_glb_acc[nState];
  for(int is=0; is<nState; is++){
    hSys_glb[is] = (TH1D*) fInSys_Lo[is]->Get("hintPP_merged");
    f_acc[is] = new TFile(Form("../acceptance_new/sys_acceptance_ups%dS_1804.root",is+1));
    hSys_glb_acc[is] = (TH1D*) f_acc[is]->Get("hcentSysPP");
    accept_sys = hSys_glb_acc[is]->GetBinContent(1);
    sys_global_pp[is] = TMath::Sqrt(hSys_glb[is]->GetBinContent(1)*hSys_glb[is]->GetBinContent(1)+accept_sys*accept_sys);
  } 
  
  sys_global_val = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp + nMB_unc*nMB_unc);
  double sys_global_y = sys_global_val; 
  double sys_global_x = 15;
  double sys_pp_1S = sys_global_pp[0];
  double sys_pp_2S = sys_global_pp[1];

  sys_pp_1S = TMath::Sqrt(sys_pp_1S*sys_pp_1S+relsys_ppInt1S*relsys_ppInt1S);
  sys_pp_2S = TMath::Sqrt(sys_pp_2S*sys_pp_2S+relsys_ppInt2S*relsys_ppInt2S);

  TBox *globalUncBox = new TBox(xmax-sys_global_x*3,1-sys_global_y,xmax-sys_global_x*2,1+sys_global_y);
  globalUncBox -> SetFillColorAlpha(kGray+2,0);
  globalUncBox -> SetLineColor(kBlack);
  globalUncBox -> SetLineWidth(2);
  globalUncBox -> Draw("l same");

  TBox *ppRefUncBox1S = new TBox(xmax-sys_global_x*2,1-sys_pp_1S,xmax-sys_global_x+1,1+sys_pp_1S);
  ppRefUncBox1S -> SetFillColor(kPink-6);
  ppRefUncBox1S -> Draw("same");

  TBox *ppRefUncBox2S = new TBox(xmax-sys_global_x,1-sys_pp_2S,xmax,1+sys_pp_2S);
  ppRefUncBox2S -> SetFillColor(kBlue-3);
  ppRefUncBox2S -> Draw("same");

  pad_diff->Update();
//  CMS_lumi( c1, iPeriod, iPos );
  CMS_lumi_raaCent( pad_diff, iPeriod, iPos );
  pad_diff->Update();
  //// --- 2nd pad!!!   
  c1->cd();
  pad_int->Draw(); 
  pad_int->cd(); 
  
  //// for int
  gRAA_int_sys[0]->GetXaxis()->SetLimits(xmin_int,xmax_int);
  gRAA_int_sys[0]->SetMinimum(0.0);
  gRAA_int_sys[0]->SetMaximum(1.3);
  gRAA_int_sys[0]->GetXaxis()->SetNdivisions(101);
  gRAA_int_sys[0]->GetXaxis()->SetLabelSize(0);
  gRAA_int_sys[0]->GetYaxis()->SetTickLength(0.03*600/xlonger);
  gRAA_int_sys[0]->GetYaxis()->SetLabelSize(0);
  
  //// syst
  for (int is=0; is<nState; is++){
    if ( is==0) { gRAA_int_sys[is]->Draw("A5"); }
    else if (is==ulstate && isArrow==true) { 
      for (int ipt=0; ipt< n3s_int ; ipt++) { //bin by bin
        box68per_int[ipt]->Draw("l"); 
      }
    }
    else { gRAA_int_sys[is]->Draw("5"); }
	}
  //// point
  for (int is=0; is<nState; is++){
    if (is==ulstate && isArrow==true) {
      for (int ipt=0; ipt< n3s_int ; ipt++) { //bin by bin
        arr95per_int[ipt]->Draw();
      }
    }
    else { gRAA_int[is]->Draw("P"); }
	}
  dashedLine(0.,1.,xmax_int,1.,1,1);
  
  pad_int->Update();
  //// draw text
  double sz_allign = 0.034;
  globtex->SetTextAlign(22); //center-center
  globtex->SetTextSize(0.038*600./xlonger);
  globtex->DrawLatex(0.5*(1-0.032*600/xlonger), sz_init-sz_step-sz_allign, "Cent.");
  globtex->DrawLatex(0.5*(1-0.032*600/xlonger), sz_init-sz_step*2-sz_allign, "0-100%");

	c1->Update();
  c1->SaveAs(Form("plots/Strickland_RAA_vs_cent_isArrow%d_all_newglobal_asym_postCWRConstrain.png",(int)isArrow));
  c1->SaveAs(Form("plots/Strickland_RAA_vs_cent_isArrow%d_all_newglobal_asym_postCWRConstrain.pdf",(int)isArrow));
/*
	///////////////////////////////////////////////////////////////////
	//// save as a root file
	TFile *outFile = new TFile("RAA_vs_cent.root", "RECREATE");
	outFile->cd();
	for (int is=0; is<nState; is++){
		gRAA_sys[is]->Write();	
		gRAA[is]->Write();	
	}
	outFile->Close();
*/	
	return;

} // end of main func.

