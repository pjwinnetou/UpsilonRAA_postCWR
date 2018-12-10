#include "Jaebeom.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent_projection.C"
#include "cutsAndBin.h"
#include "commonUtility.h"

void projection_RAA_pt_wo3S(bool isArrow=false)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 10001; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nState = 3; // Y(1S), Y(2S), and Y(3S)
  double xmax = 30.0;
//  double relsys = 0.1;

  double exsys_1s[6] =  {1., 1., 1., 1.5, 1.5, 9.};
  double exsys_2s[3] =  {2., 2.5, 10.5};
  double exsys_3s[2] =  {3., 12.};

  double theory_1s[6] = {(3.109810e-01 + 3.582435e-01)/2, (3.219182e-01 + 3.707651e-01)/2, (3.416386e-01 + 3.934574e-01)/2, (3.623386e-01 + 4.193548e-01)/2, (3.655781e-01 + 4.227464e-01)/2, (3.079565e-01 + 3.571198e-01)/2};
  double theory_2s[3] = {(8.354941e-02 + 1.054606e-01)/2, (9.655152e-02 + 1.148232e-01)/2, (9.326153e-02 + 1.165367e-01)/2};
  double theory_3s[2] = {(4.046849e-02 + 4.529941e-02)/2, (6.016030e-02 + 6.758260e-02)/2};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn[nState];
	TGraphErrors* gRAA[nState];
	TGraphErrors* gRAA_sys[nState];
  for (int is=0; is<nState; is++){
  	fIn[is] = new TFile(Form("../Ups_%d_RAA.root",is+1),"READ");
    gRAA[is]=(TGraphErrors*)fIn[is]->Get("gRAA_pt");
    gRAA_sys[is]=(TGraphErrors*)fIn[is]->Get("gRAA_pt");
    //cout << "gRAA["<<is<<"] = " <<gRAA[is] << endl;
  }
  //// read input file : syst.
  TFile* fInSys[nState];
  TFile* fInSys_woBkg[nState];
  TH1D* hSys[nState];
  TH1D* hSys_woBkg[nState];
  int npoint[nState];
  for (int is=0; is<nState; is++){
  	fInSys[is] = new TFile(Form("postCWR_mergedSys_constrain_ups%ds_asymHi.root",is+1),"READ");
    hSys[is]=(TH1D*)fInSys[is]->Get("hptRAA_merged");
    npoint[is] = hSys[is]->GetSize()-2;
  	fInSys_woBkg[is] = new TFile(Form("WOBkg_postCWR_mergedSys_constrain_ups%ds_asymHi.root",is+1),"READ");
    hSys_woBkg[is]=(TH1D*)fInSys_woBkg[is]->Get("hptRAA_merged");
    cout << "*** Y("<<is+1<<") : # of point = " << npoint[is] << endl;
  } 
  
  //// set bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;
  double sys_red = 3;
  double lumi_red = TMath::Sqrt(27);

  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
    if (npoint[is] != gRAA[is]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
    for (int ipt=0; ipt< npoint[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRAA[is]->GetPoint(ipt, pxtmp, pytmp); 
      extmp=gRAA[is]->GetErrorX(ipt);
      eytmp=gRAA[is]->GetErrorY(ipt);
      if( (is==0 && ipt<2) || (is==1 && ipt==0) ) relsys = hSys_woBkg[is]->GetBinContent(ipt+1);
      else relsys=hSys[is]->GetBinContent(ipt+1);
      cout << ipt <<"th bin RAA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      cout << ipt <<"th bin syst. = " << pytmp*relsys << endl; 
      
      if (is==0){
        gRAA[is]->SetPoint(ipt,pxtmp,theory_1s[ipt]);
        gRAA[is]->SetPointError(ipt,0,(theory_1s[ipt]*eytmp/pytmp)/lumi_red);
        gRAA_sys[is]->SetPoint(ipt, pxtmp, theory_1s[ipt]);
        gRAA_sys[is]->SetPointError(ipt, exsys_1s[ipt], theory_1s[ipt]*relsys/sys_red);
      }
      else if (is==1){
        gRAA[is]->SetPoint(ipt,pxtmp,theory_2s[ipt]);
        gRAA[is]->SetPointError(ipt,0,(theory_2s[ipt]*eytmp/pytmp)/lumi_red);
        gRAA_sys[is]->SetPoint(ipt, pxtmp, theory_2s[ipt]);
        gRAA_sys[is]->SetPointError(ipt, exsys_2s[ipt], theory_2s[ipt]*relsys/sys_red);
      }
      else if (is==2){
        gRAA[is]->SetPoint(ipt,pxtmp,theory_3s[ipt]);
        gRAA[is]->SetPointError(ipt,0,(theory_3s[ipt]*eytmp/pytmp)/lumi_red);
        gRAA_sys[is]->SetPoint(ipt, pxtmp, theory_3s[ipt]);
        gRAA_sys[is]->SetPointError(ipt, exsys_3s[ipt], theory_3s[ipt]*relsys/sys_red);
      }
    }
  }
 
  ////////////////////////////////////////////////////////////////
  int ulstate = 2; //3S
  static const int n3s = 2;
  double boxw = 0.6; // for syst. box (vs cent)
  double lower68[n3s] = {lower68_pt1,lower68_pt2};
  double upper68[n3s] = {upper68_pt1,upper68_pt2};
  double lower95[n3s] = {lower95_pt1,lower95_pt2};
  double upper95[n3s] = {upper95_pt1,upper95_pt2};
  if (n3s != npoint[ulstate]) {cout<<"ERROR!! # of bins for UL is wrong!!"<<endl;return;} 

  //// --- vs centrality
  TBox *box68per[n3s];
  TArrow *arr95per[n3s];
  for (int ipt=0; ipt< n3s ; ipt++) { //bin by bin
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; 
    //lower68=0; upper68=0; lower95=0; upper95=0; 
    gRAA[ulstate]->GetPoint(ipt, pxtmp, pytmp);
    box68per[ipt] = new TBox(pxtmp-boxw,lower68[ipt],pxtmp+boxw,upper68[ipt]);
    arr95per[ipt] = new TArrow(pxtmp,lower95[ipt],pxtmp,upper95[ipt],0.027,"<-|"); //95%
    box68per[ipt]->SetLineColor(kGreen+3);
    box68per[ipt]->SetFillColorAlpha(kGreen-6,0.5);
    box68per[ipt]->SetLineWidth(1);
    arr95per[ipt]->SetLineColor(kGreen+2);
    arr95per[ipt]->SetLineWidth(2);
  }

  //// graph style 
  for (int is=0; is<nState; is++){
    SetGraphStyle(gRAA[is], is, is); 
    SetGraphStyleSys(gRAA_sys[is], is); 
    gRAA_sys[is]-> SetFillColorAlpha(kWhite,0.);
  }
  
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  
  //// legend
  //// axis et. al
  gRAA_sys[0]->GetXaxis()->SetTitle("p_{T} (GeV)");
  gRAA_sys[0]->GetXaxis()->CenterTitle();
  gRAA_sys[0]->GetYaxis()->SetTitle("R_{AA}");
  gRAA_sys[0]->GetYaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetLimits(0.,xmax);
  gRAA_sys[0]->SetMinimum(0.0);
  gRAA_sys[0]->SetMaximum(1.14);


  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  gRAA_sys[0]->Draw("A5");
  
  dashedLine(0.,1.,xmax,1.,1,1);
  TLegend *leg= new TLegend(0.64, 0.545, 0.855, 0.705);
  SetLegendStyle(leg);
  TLegend *leg_up= new TLegend(0.64, 0.50, 0.85, 0.62);
  SetLegendStyle(leg_up);

  TArrow *arrLeg = new TArrow(18.59,0.532,18.59,0.582,0.02,"<-|");
  arrLeg->SetLineColor(kGreen+2);
  arrLeg->SetLineWidth(2);

  if (isArrow==false) { 
    for (int is=0; is<nState-1; is++){
      leg -> AddEntry(gRAA[is],Form(" #varUpsilon(%dS)",is+1),"lp");
    }
  }
  else {
    leg -> AddEntry(gRAA[0]," #varUpsilon(1S)","lp");
    leg -> AddEntry(gRAA[1]," #varUpsilon(2S)","lp");
//    leg -> AddEntry(gRAA[2]," #Upsilon(3S)","lp");
    TLegendEntry *ent=leg_up->AddEntry("ent"," #varUpsilon(3S) 68\% CL","f");
    ent->SetLineColor(kGreen+3);
    ent->SetFillColorAlpha(kGreen-6,0.5);
    ent->SetFillStyle(1001);
    ent=leg_up->AddEntry("ent"," #Upsilon(3S) 95\% CL","f");
    ent->SetLineColor(kWhite);
//    leg_up->SetTextSize(0.03);
//    leg_up->Draw("same");
//    arrLeg->Draw();
  }

    leg->Draw("same");

  //// draw text
  double sz_init = 0.925; double sz_step = 0.0535;
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu#mu} < 30 GeV/c");
  globtex->DrawLatex(0.22, sz_init-sz_step, "|y| < 2.4");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 2.4");
  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 0-100%");


  
  TFile *RappTot = new TFile("Rapp_RAA_5023_pt_tot.root","READ");
  TFile *RappReg = new TFile("Rapp_RAA_5023_pt_reg.root","READ");
  
  TGraphErrors *gRAA_tot_max[nState]; 
  TGraphErrors *gRAA_tot_min[nState]; 
  TGraphErrors *gRAA_tot_shade[nState]; 
  TGraphErrors *gRAA_reg_max[nState]; 
  TGraphErrors *gRAA_reg_min[nState]; 
  TGraphErrors *gRAA_reg_shade[nState]; 
  
  for(int i=0;i<nState;i++)
  {
    gRAA_tot_max[i] = (TGraphErrors*) RappTot-> Get(Form("RAA_rapp_pT_tot_%dS_max",i+1));
    gRAA_tot_min[i] = (TGraphErrors*) RappTot-> Get(Form("RAA_rapp_pT_tot_%dS_min",i+1));
    gRAA_tot_shade[i] = (TGraphErrors*) RappTot-> Get(Form("RAA_%ds_pt_tot_shade",i+1));
    gRAA_tot_max[i] -> SetLineWidth(2.);
    gRAA_tot_min[i] -> SetLineWidth(2.0);
    gRAA_tot_shade[i] -> SetLineWidth(2);
  }
  for(int i=0;i<nState;i++)
  {
    gRAA_reg_max[i] = (TGraphErrors*) RappReg-> Get(Form("RAA_rapp_pT_reg_%dS_max",i+1));
    gRAA_reg_min[i] = (TGraphErrors*) RappReg-> Get(Form("RAA_rapp_pT_reg_%dS_min",i+1));
    gRAA_reg_shade[i] = (TGraphErrors*) RappReg-> Get(Form("RAA_%ds_pt_reg_shade",i+1));
    gRAA_reg_max[i] -> SetLineWidth(2.);
    gRAA_reg_min[i] -> SetLineWidth(2.0);
    gRAA_reg_shade[i] -> SetLineWidth(2);
  }

  gRAA_tot_max[0]->SetLineColor(kRed+1);
  gRAA_tot_max[1]->SetLineColor(kAzure-3);
  
  gRAA_tot_min[0]->SetLineColor(kRed+1);
  gRAA_tot_min[1]->SetLineColor(kAzure-3);
   
  gRAA_tot_shade[0]->SetFillStyle(3004);
  gRAA_tot_shade[1]->SetFillStyle(3005);
  gRAA_tot_shade[0]->SetFillColor(kRed+1);
  gRAA_tot_shade[1]->SetFillColor(kAzure-3);

  gRAA_reg_max[0]->SetLineColor(kOrange-3);
  gRAA_reg_max[1]->SetLineColor(kTeal-5);
  
  gRAA_reg_min[0]->SetLineColor(kOrange-3);
  gRAA_reg_min[1]->SetLineColor(kTeal-5);
   
  gRAA_reg_shade[0]->SetFillStyle(3004);
  gRAA_reg_shade[1]->SetFillStyle(3004);
  gRAA_reg_shade[0]->SetFillColor(kOrange-3);
  gRAA_reg_shade[1]->SetFillColor(kTeal-5);
  
  gRAA_tot_shade[0]->SetLineColor(kRed+1);
  gRAA_tot_shade[1]->SetLineColor(kBlue);
  gRAA_reg_shade[0]->SetLineColor(kOrange-3);
  gRAA_reg_shade[1]->SetLineColor(kTeal-5);

  gRAA_tot_shade[0]->SetFillColorAlpha(kRed+1,0.8);
  gRAA_tot_shade[1]->SetFillColorAlpha(kAzure-3,0.8);

  for(int is=0;is<nState-1;is++){
    gRAA_tot_shade[is]->Draw("f");
    gRAA_tot_max[is]->Draw("l");
    gRAA_tot_min[is]->Draw("l");
    gRAA_sys[is]->Draw("5");
    gRAA[is]->Draw("P");
  }
  
  double line_y = 0.827;
  double line_y_diff = 0.07;
  double line_y_diff_in = 0.02;
  double line_x_end = 4.4;
  double line_x_start = 2.4;
  double leg_y_text_diff = 0.095;
  double leg_y_diff = 0.23;

  TLegend *leg_strick= new TLegend(0.22, line_y-leg_y_text_diff-leg_y_diff+0.05, 0.42, line_y-leg_y_text_diff-0.03);
  SetLegendStyle(leg_strick);
  leg_strick->SetTextSize(0.038);
  leg_strick->AddEntry(gRAA_tot_shade[0]," #varUpsilon(1S) Total","f");
  leg_strick->AddEntry(gRAA_tot_shade[1]," #varUpsilon(2S) Total","f");
  leg_strick->Draw("same");
 
  drawText2("X. Du, M. He, R. Rapp",line_x_start,line_y+0.034,20);


  //Global Unc.
  double TAA_unc_Global_Hi = 0.028;
  double TAA_unc_Global_Lo = 0.034;

  double sys_global_val_Hi = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+TAA_unc_Global_Hi*TAA_unc_Global_Hi+nMB_unc*nMB_unc);
  double sys_global_val_Lo = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+TAA_unc_Global_Lo*TAA_unc_Global_Lo+nMB_unc*nMB_unc);
  double sys_global_y_Hi = sys_global_val_Hi;
  double sys_global_y_Lo = sys_global_val_Lo;
  double sys_global_x = 0.8;
  TBox *globalUncBox = new TBox(xmax-sys_global_x,1-sys_global_y_Lo,xmax,1+sys_global_y_Hi);
  globalUncBox -> SetLineColor(kBlack);
  globalUncBox -> SetFillColorAlpha(kGray+2,0.6);
  globalUncBox -> SetLineWidth(1);
  globalUncBox -> Draw("l same");
  
  CMS_lumi_raaCent_projection( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs("Projection_RAA_vs_pt_wo3S.pdf");
  c1->SaveAs("Projection_RAA_vs_pt_wo3S.png");

/*
	///////////////////////////////////////////////////////////////////
	//// save as a root file
	TFile *outFile = new TFile("RAA_vs_pt.root", "RECREATE");
	outFile->cd();
	for (int is=0; is<nState; is++){
		gRAA_sys[is]->Write();	
		gRAA[is]->Write();	
	}
	outFile->Close();
*/	
	return;

} // end of main func.

