#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "../cutsAndBin.h"
#include "../commonUtility.h"

void strickland_RAA_3S_cent_postCWRConstrain(bool isArrow =true)
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

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn[nState];
	TGraphErrors* gRAA[nState]; // vs centrality
  TGraphErrors* gRAA_sys[nState];
	TGraphErrors* gRAA_int[nState]; // centrality-integrated
  TGraphErrors* gRAA_int_sys[nState];
  for (int is=0; is<nState; is++){
  	fIn[is] = new TFile(Form("Ups_%d_RAA.root",is+1),"READ");
    gRAA[is]=(TGraphErrors*)fIn[is]->Get("gRAA_cent");
    gRAA_sys[is]=(TGraphErrors*)fIn[is]->Get("gRAA_cent");
    gRAA_int[is]=(TGraphErrors*)fIn[is]->Get("gRAA_int");
    gRAA_int_sys[is]=(TGraphErrors*)fIn[is]->Get("gRAA_int");
    //cout << "gRAA["<<is<<"] = " <<gRAA[is] << endl;
  }
  //// read input file : syst.
  TFile* fInSys[nState];
  TH1D* hSys[nState];
  TH1D* hSys_int[nState];
  int npoint[nState];
  int npoint_int[nState];
  for (int is=0; is<nState; is++){
  	fInSys[is] = new TFile(Form("../Systematic/postCWR_mergedSys_constrain_ups%ds_asymHi.root",is+1),"READ");
    hSys[is]=(TH1D*)fInSys[is]->Get("hcentRAA_merged");
    npoint[is] = hSys[is]->GetSize()-2;
    cout << "*** CENT *** Y("<<is+1<<") : # of point = " << npoint[is] << endl;
    hSys_int[is]=(TH1D*)fInSys[is]->Get("hintRAA_merged");
    npoint_int[is] = hSys_int[is]->GetSize()-2;
    cout << "*** INT *** Y("<<is+1<<") : # of point = " << npoint_int[is] << endl;
  }   
  
  //// set bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;
  double xshift = 0.2;
  //// --- vs centrality
  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
    if (npoint[is] != gRAA[is]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
    for (int ipt=0; ipt< npoint[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRAA[is]->GetPoint(ipt, pxtmp, pytmp);
      extmp=gRAA[is]->GetErrorX(ipt);
      eytmp=gRAA[is]->GetErrorY(ipt);
      relsys=hSys[is]->GetBinContent(npoint[is]-ipt);
      cout << ipt <<"th bin RAA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      //cout << ipt <<"th bin rel. syst. = " << relsys << endl;
      cout << ipt <<"th bin syst. = " << pytmp*relsys << endl; 
      //// 1) remove ex from gRAA
      gRAA[is]->SetPointError(ipt, 0, eytmp);
      //// 2) set ey for gRAA_sys 
      //gRAA_sys[is]->SetPointError(ipt, extmp, pytmp*relsys);
      gRAA_sys[is]->SetPointError(ipt, boxw, pytmp*relsys);
    }
  }
  //// --- centrality-integrated
  cout << " " << endl;
  cout << " INTEGRATED" << endl;
  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
    if (npoint_int[is] != gRAA_int[is]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }    
    for (int ipt=0; ipt< npoint_int[is] ; ipt++) {
      pxtmp=0; pytmp=0; extmp=0; eytmp=0;
      gRAA_int[is]->GetPoint(ipt, pxtmp, pytmp);
      extmp=gRAA_int[is]->GetErrorX(ipt);
      eytmp=gRAA_int[is]->GetErrorY(ipt);
      relsys=hSys_int[is]->GetBinContent(ipt+1);
      cout << ipt <<"th bin RAA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      //cout << ipt <<"th bin rel. syst. = " << relsys << endl;
      cout << ipt <<"th bin syst. = " << pytmp*relsys << endl; 
      //// 1) remove ex from gRAA
      gRAA_int[is]->SetPoint(ipt, pxtmp+xshift*is, pytmp);
      gRAA_int[is]->SetPointError(ipt, 0, eytmp);
      //// 2) set ey for gRAA_sys
      gRAA_int_sys[is]->SetPoint(ipt, pxtmp+xshift*is, pytmp);
      //gRAA_int_sys[is]->SetPointError(ipt, extmp, pytmp*relsys);
      gRAA_int_sys[is]->SetPointError(ipt, boxw_int, pytmp*relsys); //extemp fixed
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
  /*
  double lower68[n3s] = {0., 0., 0.0183, 0.};
  double upper68[n3s] = {0.154250585, 0.01689862101, 0.0943 , 0.03990571614};
  double lower95[n3s] = {0., 0., 0., 0.};
  double upper95[n3s] = {0.22041872, 0.07136456236, 0.1333, 0.06203828951};
  static const int n3s_int = 1;
  double lower68_int[n3s_int] = {0.};
  double upper68_int[n3s_int] = {0.1536360214};
  double lower95_int[n3s_int] = {0.};
  double upper95_int[n3s_int] = {0.2477551813};
 */
  if (n3s != npoint[ulstate]) {cout<<"ERROR!! # of bins for UL is wrong!!"<<endl;return;} 
  if (n3s_int != npoint_int[ulstate]) {cout<<"ERROR!! # of bins for UL (int) is wrong!!"<<endl;return;} 

  //// --- vs centrality
  TBox *box68per[n3s];
  TArrow *arr95per[n3s];
  for (int ipt=0; ipt< n3s ; ipt++) { //bin by bin
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; 
    //lower68=0; upper68=0; lower95=0; upper95=0; 
    gRAA[ulstate]->GetPoint(ipt, pxtmp, pytmp);
    box68per[ipt] = new TBox(pxtmp-boxw,lower68[ipt],pxtmp+boxw,upper68[ipt]);
    arr95per[ipt] = new TArrow(pxtmp,lower95[ipt],pxtmp,upper95[ipt],0.027,"<-|"); //95%
    box68per[ipt]->SetLineColor(kGreen+2);
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
  
  ////////////////////////////////////////////////////////////////
  
  //// graph style 
  for (int is=0; is<nState; is++){
    SetGraphStyle(gRAA[is], is, is); 
    SetGraphStyleSys(gRAA_sys[is], is); 
    SetGraphStyle(gRAA_int[is], is, is); 
    SetGraphStyleSys(gRAA_int_sys[is], is); 
	}
  
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);


  //dashedLine(0.,1.,xmax,1.,1,1);
  //function from SONGKYO still need to draw it t1
  TF1 *tx = new TF1("f1","1",0,xmax);
  //TLine* tx = new TLine(0,1,xmax,1);
      tx->SetLineWidth(1);
         tx->SetLineStyle(7);
            tx->SetLineColor(1);

  //// axis et. al
  tx->GetXaxis()->SetTitle("N_{part}");
  tx->GetXaxis()->CenterTitle();
  tx->GetYaxis()->SetTitle("R_{AA}");
  tx->GetYaxis()->CenterTitle();
  tx->GetXaxis()->SetLimits(0.,xmax);
  tx->SetMinimum(0.0);
  tx->SetMaximum(1.);
  //for cent
  tx->GetXaxis()->SetTitleSize(0.06*1.0);
  tx->GetYaxis()->SetTitleSize(0.06*1.0);
  tx->GetXaxis()->SetLabelSize(0.05*1.0);
  tx->GetYaxis()->SetLabelSize(0.05*1.0);
  
  //// draw  
  double xlonger = 120;
  TCanvas* c1 = new TCanvas("c1","c1",600+xlonger,600);
  TPad* pad_diff = new TPad("pad_diff", "",0, 0, 600/(600.+xlonger), 1.0); // vs centrality
  pad_diff->SetRightMargin(0);
  pad_diff->Draw(); 
  
  TPad* pad_int = new TPad("pad_int", "",600/(600.+xlonger), 0, 1.0, 1.0); // centrality-integrated
  pad_int->SetLeftMargin(0);
  pad_int->SetRightMargin(0.032*600/xlonger);
  pad_int->Draw(); 
  //// --- 1st pad!!!   
  //c1->cd();
  pad_diff->cd(); 
  //// syst
  tx->Draw();
  for (int is=2; is<nState; is++){
    if ( is==0) { gRAA_sys[is]->Draw("A5"); }
    else if (is==ulstate && isArrow==true) { 
      for (int ipt=0; ipt< n3s ; ipt++) { //bin by bin
        box68per[ipt]->Draw("l"); 
      }
    }
    else { gRAA_sys[is]->Draw("5"); }
	}
  //// point
  for (int is=2; is<nState; is++){
    if (is==ulstate && isArrow==true) {
      for (int ipt=0; ipt< n3s ; ipt++) { //bin by bin
        arr95per[ipt]->Draw();
      }
    }
    else { gRAA[is]->Draw("P"); }
	}
  
  //// legend
  //TLegend *leg= new TLegend(0.75, 0.50, 0.95, 0.70);
  //TLegend *leg= new TLegend(0.65, 0.51, 0.85, 0.76);
  TLegend *leg= new TLegend(0.68, 0.60, 0.88, 0.716);
  SetLegendStyle(leg);
  leg->SetTextSize(0.040);
  TArrow *arrLeg = new TArrow(270.,0.58,270.,0.63,0.025,"<-|");
  arrLeg->SetLineColor(kGreen+2);
  arrLeg->SetLineWidth(2);
  
  if (isArrow==false) { 
    for (int is=0; is<nState; is++){
      leg -> AddEntry(gRAA[is],Form(" #Upsilon(%dS)",is+1),"lp");
    }
  }
  else {
    //leg -> AddEntry(gRAA[0]," #Upsilon(1S)","lp");
    //leg -> AddEntry(gRAA[1]," #Upsilon(2S)","lp");
    TLegendEntry *ent=leg->AddEntry("ent"," #Upsilon(3S) 68\% CL","f");
    ent->SetLineColor(kGreen+2);
    ent->SetFillColorAlpha(kGreen-10,0.5);
    ent->SetFillStyle(1001);
    ent=leg->AddEntry("ent"," #Upsilon(3S) 95\% CL","f");
    ent->SetLineColor(kWhite);
  }
  leg->Draw("same");
  arrLeg->Draw();
  
  //// draw text
  double sz_init = 0.889; double sz_step = 0.0558;
//  globtex->DrawLatex(0.22+0.04, sz_init, "p_{T}^{#mu} > 4 GeV/c");
  globtex->DrawLatex(0.22+0.04, sz_init, "p_{T} < 30 GeV");
//  globtex->DrawLatex(0.46+0.04, sz_init+0.002, "|#eta|^{#mu} < 2.4");
  globtex->DrawLatex(0.22+0.04, sz_init-sz_step, "|y| < 2.4");

/*  TLatex* centtex = new TLatex();
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

  globtex->DrawLatex(0.24, sz_init-sz_step*10+0.04, "30-100%");
  globtex->DrawLatex(0.69, sz_init-sz_step*12+0.01, "0-30%");

  TFile *fstrickland = new TFile("TheoryCurve/StrickLand_RAA_5023.root","READ");
  
  TGraphErrors *gRAA_1S_strickland[3]; 
  TGraphErrors *gRAA_2S_strickland[3]; 
  
  for(int i=0;i<3;i++)
  {
    gRAA_1S_strickland[i] = (TGraphErrors*) fstrickland-> Get(Form("RAA_strick_nPart_1S_%d",i));
    gRAA_2S_strickland[i] = (TGraphErrors*) fstrickland-> Get(Form("RAA_strick_nPart_2S_%d",i));
    gRAA_1S_strickland[i] -> SetLineWidth(3.);
    gRAA_2S_strickland[i] -> SetLineWidth(3.0);
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
  

  for(int i=0;i<3;i++){
    //gRAA_1S_strickland[i]->Draw("same");
    //gRAA_2S_strickland[i]->Draw("same");
  }
   
//----3S strickland
const int npt = 407;

TGraphErrors* g3S = new TGraphErrors("TheoryCurve/Y3SNpart5023Xi0.tsv","%lg %lg %lg %lg","\t");

int nPts = g3S->GetN();
cout<<"Number of Points: "<<nPts<<endl;
double y3S_1[npt];
double y3S_2[npt];
double y3S_3[npt];


for(int i=0;i< nPts; i++){
    y3S_1[i]=g3S->GetY()[i];
    y3S_2[i]=g3S->GetErrorX(i);
    y3S_3[i]=g3S->GetErrorY(i);
}

TGraph* g1t = new TGraph(g3S->GetN(),g3S->GetX(),y3S_1);
TGraph* g2t = new TGraph(g3S->GetN(),g3S->GetX(),y3S_2);
TGraph* g3t = new TGraph(g3S->GetN(),g3S->GetX(),y3S_3);
g1t->Draw("lsame");
g1t->SetLineWidth(3);
g1t->SetLineStyle(3);
g1t->SetLineColor(kGreen+2);
g2t->Draw("SAMEl");
g2t->SetLineWidth(3);
g2t->SetLineStyle(1);
g2t->SetLineColor(kGreen+2);
g3t->Draw("SAMEl");
g3t->SetLineWidth(3);
g3t->SetLineColor(kGreen+2);
g3t->SetLineStyle(8);



  TLegend *leg_strick= new TLegend(0.29, 0.586, 0.46, 0.716);
  SetLegendStyle(leg_strick);
  leg_strick->SetTextSize(0.04);
  leg_strick->AddEntry(gRAA_1S_strickland[2],"Y(1S)","l");
  leg_strick->AddEntry(gRAA_2S_strickland[2],"Y(2S)","l");
//  leg_strick->Draw("same");

  double line_y = 0.62;
  double line_y_diff = 0.07;
  double line_x_end = 75;
  double line_x_start = 50;
  TLine* t1 = new TLine(line_x_start,line_y,line_x_end,line_y);
  t1->SetLineStyle(3);
  t1->SetLineWidth(2);
  t1->SetLineColor(kGreen+2);
  t1->Draw("same");

  TLine* t2 = new TLine(line_x_start,line_y-line_y_diff,line_x_end,line_y-line_y_diff);
  t2->SetLineStyle(1);
  t2->SetLineWidth(2);
  t2->SetLineColor(kGreen+2);
  t2->Draw("same");

  TLine* t3 = new TLine(line_x_start,line_y-line_y_diff*2,line_x_end,line_y-line_y_diff*2);
  t3->SetLineStyle(8);
  t3->SetLineWidth(2);
  t3->SetLineColor(kGreen+2);
  t3->Draw("same");

  drawText2("4#pi #eta/s=1", line_x_end+7, line_y-0.015, 22);
  drawText2("4#pi #eta/s=2", line_x_end+7, line_y-line_y_diff*1-0.015, 22);
  drawText2("4#pi #eta/s=3", line_x_end+7, line_y-line_y_diff*2-0.015, 22);
  drawText2("Krouppa, Strickland",line_x_start,line_y+0.05,22);

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
    hSys_glb[is] = (TH1D*) fInSys[is]->Get("hintPP_merged");
    f_acc[is] = new TFile(Form("../acceptance/sys_acceptance_ups%dS_1804.root",is+1));
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
  //globalUncBox -> Draw("l same");

  TBox *ppRefUncBox1S = new TBox(xmax-sys_global_x*2,1-sys_pp_1S,xmax-sys_global_x+1,1+sys_pp_1S);
  ppRefUncBox1S -> SetFillColor(kPink-6);
 // ppRefUncBox1S -> Draw("same");

  TBox *ppRefUncBox2S = new TBox(xmax-sys_global_x,1-sys_pp_2S,xmax,1+sys_pp_2S);
  ppRefUncBox2S -> SetFillColor(kBlue-3);
 // ppRefUncBox2S -> Draw("same");

//  CMS_lumi( c1, iPeriod, iPos );
  CMS_lumi( pad_diff, iPeriod, iPos);
 
  //// --- 2nd pad!!!   
  pad_diff->Update();
  c1->Update();
  pad_int->cd(); 
  
  //// for int
  TF1 *tI = new TF1("f1","1",0,xmax_int);
  //TLine* tx = new TLine(0,1,xmax,1);
  tI->SetLineWidth(1);
  tI->SetLineStyle(7);
  tI->SetLineColor(1);

  tI->GetXaxis()->SetLimits(xmin_int,xmax_int);
  tI->SetMinimum(0.0);
  tI->SetMaximum(1.);
  tI->GetXaxis()->SetNdivisions(101);
  tI->GetXaxis()->SetLabelSize(0);
  tI->GetYaxis()->SetTickLength(0.03*600/xlonger);
  tI->GetYaxis()->SetLabelSize(0);
  tI->Draw(); 
  //// syst
  for (int is=2; is<nState; is++){
    if ( is==0) { gRAA_int_sys[is]->Draw("A5"); }
    else if (is==ulstate && isArrow==true) { 
      for (int ipt=0; ipt< n3s_int ; ipt++) { //bin by bin
        box68per_int[ipt]->Draw("l"); 
      }
    }
    else { gRAA_int_sys[is]->Draw("5"); }
	}
  //// point
  for (int is=2; is<nState; is++){
    if (is==ulstate && isArrow==true) {
      for (int ipt=0; ipt< n3s_int ; ipt++) { //bin by bin
        arr95per_int[ipt]->Draw();
      }
    }
    else { gRAA_int[is]->Draw("P"); }
	}
  
  double sz_d = 0.1391;
  //// draw text
  globtex->SetTextAlign(22); //center-center
  globtex->SetTextSize(0.038*600./xlonger);
  globtex->DrawLatex(0.5*(1-0.032*600/xlonger), sz_init-sz_step-sz_d, "Cent.");
  globtex->DrawLatex(0.5*(1-0.032*600/xlonger), sz_init-sz_step*2-sz_d, "0-100%");
 pad_int->Update();
 c1->Update();
 

 c1->SaveAs("plots/Strickland_RAA_vs_cent_3S_postCWRConstrain.png");
 c1->SaveAs("plots/Strickland_RAA_vs_cent_3S_postCWRConstrain.pdf");
// c1->SaveAs("Strickland_RAA_vs_cent_3S.C");

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

