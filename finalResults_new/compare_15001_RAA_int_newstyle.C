#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"

void compare_15001_RAA_int_newstyle() //1 or 2 (1S or 2S)
{
  setTDRStyle();
  writeExtraText = false;       // if extra text
  int iPeriod = 101; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nfile = 5; // 0: 15001, 1: ours
  double xmax = 2.85;
//  double relsys = 0.1;
  
  double exsys_1s[6] =  {1., 1., 1., 1.5, 1.5, 9.};
  double exsys_2s[3] =  {2., 2.5, 10.5};
  double exsys_3s[2] =  {3.,12.};

  double exsys = 0.05;
  double exsys_align = 0.075;

  //// 15001 values
  const int cn_1s =  3;
  double cpx_1s[cn_1s] =  {0.51-exsys_align, 1.425-exsys_align, 2.3-exsys_align};
  double cpx_1s_exsys[cn_1s] = {0.51+exsys_align, 1.425+exsys_align, 2.3+exsys_align};
  double cpy_1s[cn_1s] =  {0.453, 0.119, 0.145}; 
  double cex_1s[cn_1s] =  {0., 0., 0};
  double cey_1s[cn_1s] =  {0.014, 0.028, 0.031};
  double cexsys_1s[cn_1s] =  {exsys, exsys, exsys};
  double ceysys_1s[cn_1s] =  {0.046, 0.015, 0.058};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
	TGraphErrors* gRAA[nfile];
	TGraphAsymmErrors* gRAA_sys[nfile];
	TGraphErrors* gRAA_new[nfile];
	TGraphAsymmErrors* gRAA_sys_new[nfile];
  //// 1) 15001
  /*
  gRAA[0] = new TGraphErrors(cn_1s, cpx_1s, cpy_1s, cex_1s, cey_1s); 
  gRAA_sys[0] = new TGraphAsymmErrors(cn_1s, cpx_1s, cpy_1s, cexsys_1s, cexsys_1s, ceysys_1s, ceysys_1s); 
  */
  for(int i=0;i<3;i++){
    gRAA[i] = new TGraphErrors();
    gRAA_sys[i] = new TGraphAsymmErrors();
    gRAA[i]->SetPoint(0,cpx_1s[i],cpy_1s[i]);
    gRAA[i]->SetPointError(0,cex_1s[i],cey_1s[i]);
    gRAA_sys[i]->SetPoint(0,cpx_1s[i],cpy_1s[i]);
    gRAA_sys[i]->SetPointError(0,cexsys_1s[i],cexsys_1s[i],ceysys_1s[i],ceysys_1s[i]);
  }


  //// 2) ours
  TFile* fIn_1S = new TFile("Ups_1_RAA.root","READ");
  TFile* fIn_2S = new TFile("Ups_2_RAA.root","READ");
  TFile* fIn_3S = new TFile("Ups_3_RAA.root","READ");
  gRAA_new[1]=(TGraphErrors*)fIn_1S->Get("gRAA_int");
  gRAA_sys_new[1]= new TGraphAsymmErrors();
  gRAA_new[2]=(TGraphErrors*)fIn_2S->Get("gRAA_int");
  gRAA_sys_new[2]= new TGraphAsymmErrors();
  gRAA_new[3]=(TGraphErrors*)fIn_3S->Get("gRAA_int");
  gRAA_sys_new[3]= new TGraphAsymmErrors();
  gRAA_new[4]= new TGraphErrors(cn_1s, cpx_1s_exsys, cpy_1s, cex_1s, cey_1s);
  gRAA_sys_new[4]= new TGraphAsymmErrors(cn_1s, cpx_1s_exsys, cpy_1s, cex_1s, cex_1s, cey_1s, cey_1s);
  //// read input file : syst.
  TFile* fInSys_1S_Hi = new TFile("../Systematic_new/postCWR_mergedSys_constrain_ups1s_asymHi.root","READ");
  TFile* fInSys_2S_Hi = new TFile("../Systematic_new/postCWR_mergedSys_constrain_ups2s_asymHi.root","READ");
  TFile* fInSys_3S_Hi = new TFile("../Systematic_new/postCWR_mergedSys_constrain_ups3s_asymHi.root","READ");
  TFile* fInSys_1S_Lo = new TFile("../Systematic_new/postCWR_mergedSys_constrain_ups1s_asymLo.root","READ");
  TFile* fInSys_2S_Lo = new TFile("../Systematic_new/postCWR_mergedSys_constrain_ups2s_asymLo.root","READ");
  TFile* fInSys_3S_Lo = new TFile("../Systematic_new/postCWR_mergedSys_constrain_ups3s_asymLo.root","READ");
  TH1D* hSys_1S_Hi = (TH1D*)fInSys_1S_Hi->Get("hintRAA_merged");
  TH1D* hSys_2S_Hi = (TH1D*)fInSys_2S_Hi->Get("hintRAA_merged");
  TH1D* hSys_3S_Hi = (TH1D*)fInSys_3S_Hi->Get("hintRAA_merged");
  TH1D* hSys_1S_Lo = (TH1D*)fInSys_1S_Lo->Get("hintRAA_merged");
  TH1D* hSys_2S_Lo = (TH1D*)fInSys_2S_Lo->Get("hintRAA_merged");
  TH1D* hSys_3S_Lo = (TH1D*)fInSys_3S_Lo->Get("hintRAA_merged");
  
  //// set bin width and calculate systematic uncertainties 
  double pxtmp, pytmp, extmp, eytmp;
  double relsys_Hi, relsys_Lo;
  int npoint = 3;//gRAA[0]->GetN();
 // if (npoint !=cn_1s) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
  for (int ipt=0; ipt< npoint; ipt++) 
  {
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys_Hi=0; relsys_Lo;
    gRAA_new[ipt+1]->GetPoint(0, pxtmp, pytmp);
    extmp=gRAA_new[ipt+1]->GetErrorX(0);
    eytmp=gRAA_new[ipt+1]->GetErrorY(0);
    if(ipt==0){
      relsys_Hi=hSys_1S_Hi->GetBinContent(1);
      relsys_Lo=hSys_1S_Lo->GetBinContent(1);
    }
    else if(ipt==1){
      relsys_Hi=hSys_2S_Hi->GetBinContent(1);
      relsys_Lo=hSys_2S_Lo->GetBinContent(1);
    }
    else if(ipt==2){
      relsys_Hi=hSys_3S_Hi->GetBinContent(1);
      relsys_Lo=hSys_3S_Lo->GetBinContent(1);
    }
    relsys_Hi=TMath::Sqrt(relsys_Hi*relsys_Hi + lumi_unc_pp*lumi_unc_pp + nMB_unc*nMB_unc);
    relsys_Lo=TMath::Sqrt(relsys_Lo*relsys_Lo + lumi_unc_pp*lumi_unc_pp + nMB_unc*nMB_unc);
    gRAA_new[ipt+1]->SetPoint(0, cpx_1s_exsys[ipt], pytmp);
    gRAA_sys_new[ipt+1]->SetPoint(0, cpx_1s_exsys[ipt], pytmp);
    gRAA_new[ipt+1]->SetPointError(0, 0, eytmp);
    gRAA_sys_new[ipt+1]->SetPointError(0, exsys, exsys, pytmp*relsys_Lo, pytmp*relsys_Hi);
  
    cout << "ipt+1 : " << ipt +1 << endl;
    cout << "RAA of Y : " << eytmp << endl;
  }
 

  gRAA[2]->SetPoint(0,-100,-100); 
  gRAA_sys[2]->SetPoint(0,-100,-100); 
  gRAA_new[3]->SetPoint(0,-100,-100); 
  gRAA_sys_new[3]->SetPoint(0,-100,-100); 
  ////////////////////////////////////////////////////////////////


  //******************************************
  // Upper Limit
  //******************************************
  
  const int numComp = 2;
  double lower95_int[numComp] = {0,lower95_cint}; 
  double upper95_int[numComp] = {0.145,upper95_cint};
  
  double align_upper = 0.558;
  double exsys_upper = 0.15;

   TArrow *arr95per_int[numComp];
   for(int icomp = 0; icomp<numComp; icomp++)
   {
     arr95per_int[icomp] = new TArrow((cpx_1s[cn_1s-1]+cpx_1s[cn_1s-2])/2-exsys/2+exsys_upper*icomp+align_upper,lower95_int[icomp],(cpx_1s[cn_1s-1]+cpx_1s[cn_1s-2])/2-exsys/2+exsys_upper*icomp+align_upper,upper95_int[icomp],0.027,"<-|");
     arr95per_int[icomp]->SetLineWidth(2);
     if(icomp==0) arr95per_int[icomp]->SetLineColor(kBlack);
     else if(icomp!=0) arr95per_int[icomp]->SetLineColor(kGreen+2);
   }


  //// graph style 
  SetGraphStyleOpen(gRAA[0], 4, 4, 0); 
  SetGraphStyleSys(gRAA_sys[0], 4); 
  SetGraphStyleOpen(gRAA[1], 4, 4, 1); 
  SetGraphStyleSys(gRAA_sys[1], 4); 
  SetGraphStyleOpen(gRAA[2], 4, 4, 2); 
  SetGraphStyleSys(gRAA_sys[2], 4); 
  SetGraphStyle(gRAA_new[1], 0, 0); 
  SetGraphStyleSys(gRAA_sys_new[1], 0); 
  SetGraphStyle(gRAA_new[2], 1, 1); 
  SetGraphStyleSys(gRAA_sys_new[2], 1); 
  SetGraphStyle(gRAA_new[3], 2, 2); 
  SetGraphStyleSys(gRAA_sys_new[3], 2); 
  

  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  
  TLatex* globtex_label = new TLatex();
  globtex_label->SetNDC();
  globtex_label->SetTextAlign(12); //left-center
  globtex_label->SetTextFont(42);
  globtex_label->SetTextSize(0.047);
  
  //// axis et. al
  gRAA_sys[0]->GetXaxis()->SetTitle("");
  gRAA_sys[0]->GetXaxis()->CenterTitle();
  gRAA_sys[0]->GetYaxis()->SetTitle("R_{AA}");
  gRAA_sys[0]->GetYaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetLimits(0.,xmax);
  gRAA_sys[0]->SetMinimum(0.0);
  gRAA_sys[0]->SetMaximum(1.);
 
  TArrow *arrLeg_l = new TArrow(17.4,0.0093,17.4,0.0145,0.021,"<-|");
  arrLeg_l->SetLineColor(kGreen+2);
  arrLeg_l->SetLineWidth(2);

  for(int i=0;i<3;i++){
    gRAA[i]->GetXaxis()->SetBinLabel(10,"");
    gRAA_sys[i]->GetXaxis()->SetBinLabel(10,"");
    gRAA_new[i+1]->GetXaxis()->SetBinLabel(10,"");
    gRAA_sys_new[i+1]->GetXaxis()->SetBinLabel(10,"");
  }
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.1);
  gPad->SetTopMargin(0.067);
  for(int i=0;i<npoint;i++){
    if(i==0) gRAA_sys[i]->Draw("A5");
    else if(i!=0) gRAA_sys[i]->Draw("5");
    gRAA[i]->Draw("P");
    gRAA_sys_new[i+1]->Draw("5");
    gRAA_new[i+1]->Draw("P");
  }
  arr95per_int[0]->Draw();
  arr95per_int[1]->Draw();
  dashedLine(0.,1.,xmax,1.,1,1);
 
  double legend_ypos1 = 0.6; 
  double legend_ypos2 = 0.87; 
  //// legend
  TLegend *leg= new TLegend(0.474, legend_ypos1, 0.904, legend_ypos2);
  SetLegendStyle(leg);
  leg->SetNColumns(2);
  leg -> SetHeader("");
  //leg -> SetHeader("#Upsilon's");
  leg -> AddEntry(gRAA_new[1]," ","lp");
  leg -> AddEntry(gRAA_new[2],"      #sqrt{s_{NN}} = 5.02 TeV","lp");
//  leg -> AddEntry(gRAA[0],"#sqrt{s_{NN}} = 2.76 TeV","lp");

  double arrLegpos_x = 1.55;
  double arrLegpoydiff1 = 0.055;
  double arrLegpoydiff2 = 0.097;
  TArrow *arrLeg = new TArrow(arrLegpos_x ,legend_ypos1+arrLegpoydiff1, arrLegpos_x,legend_ypos1+arrLegpoydiff2,0.02,"<-|");
  arrLeg->SetLineColor(kGreen+2);
  arrLeg->SetLineWidth(2);
//    leg_up->SetTextSize(0.03);

  TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
  header->SetTextSize(0.046);
  header->SetTextFont(62);
  leg->Draw("same");
  arrLeg->Draw();
  
  double legend_yposdiff1 = 0.055;
  double legend_yposdiff2 = 0.135;
  TLegend *leg276= new TLegend(0.474, legend_ypos1-legend_yposdiff1, 0.904, legend_ypos2-legend_yposdiff2);
  SetLegendStyle(leg276);
  leg276->SetNColumns(2);
  leg276 -> SetHeader("");
  //leg -> SetHeader("#Upsilon's");
  leg276 -> AddEntry(gRAA[0]," ","lp");
  leg276 -> AddEntry(gRAA[1],"      #sqrt{s_{NN}} = 2.76 TeV","lp");
//  leg -> AddEntry(gRAA[0],"#sqrt{s_{NN}} = 2.76 TeV","lp");

  double arrLegpos_x276 = 1.55;
  double arrLegpo276ydiff1 = 0.025;
  double arrLegpo276ydiff2 = 0.065;
  TArrow *arrLeg276 = new TArrow(arrLegpos_x276 ,legend_ypos1-legend_yposdiff1+arrLegpo276ydiff1, arrLegpos_x276, legend_ypos1-legend_yposdiff1+arrLegpo276ydiff2,0.02,"<-|");
  arrLeg276->SetLineColor(kBlack);
  arrLeg276->SetLineWidth(2);
//    leg_up->SetTextSize(0.03);

  TLegendEntry *header276 = (TLegendEntry*)leg276->GetListOfPrimitives()->First();
  header276->SetTextSize(0.046);
  header276->SetTextFont(62);
  leg->Draw("same");
  arrLeg->Draw();
  leg276->Draw("same");
  arrLeg276->Draw();
  //// draw text
  double sz_init = 0.87; double sz_step = 0.0535;
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
  globtex->DrawLatex(0.79,sz_init-sz_step*11.5,"95% CL");
  globtex->DrawLatex(0.22, sz_init, "p_{T} < 30 GeV");
  globtex->DrawLatex(0.22, sz_init-sz_step-0.007, "|y| < 2.4");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 2.4");
  globtex->DrawLatex(0.22, sz_init-sz_step*2.2, "Cent. 0-100%");
  globtex_label->DrawLatex(0.243, sz_init-sz_step*15.24, "#varUpsilon(1S)");
  globtex_label->DrawLatex(0.505, sz_init-sz_step*15.24, "#varUpsilon(2S)");
  globtex_label->DrawLatex(0.782, sz_init-sz_step*15.24, "#varUpsilon(3S)");
  
  double TAA_unc_Global_Hi = 0.028;
  double TAA_unc_Global_Lo = 0.034;

  double sys_global_val_Hi = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+TAA_unc_Global_Hi*TAA_unc_Global_Hi+nMB_unc*nMB_unc);
  double sys_global_val_Lo = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+TAA_unc_Global_Lo*TAA_unc_Global_Lo+nMB_unc*nMB_unc);
  double sys_global_y_Hi = sys_global_val_Hi;
  double sys_global_y_Lo = sys_global_val_Lo;
  double sys_global_x = 0.8;
  //double sys_global_val = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_aa*lumi_unc_aa);
  double sys_global_y_15001 = 0.079;
  TBox *globalUncBox = new TBox(xmax-sys_global_x*2,1-sys_global_y_Lo,xmax-sys_global_x-0.05,1+sys_global_y_Hi);
  globalUncBox -> SetLineColor(kRed-2);
  globalUncBox -> SetFillColorAlpha(kPink-6,0.6);
  globalUncBox -> SetLineWidth(1);
  //globalUncBox -> Draw("l same");
  
  TBox *globalUncBox_15001 = new TBox(xmax-sys_global_x,1-sys_global_y_15001,xmax,1+sys_global_y_15001);
  globalUncBox_15001 -> SetLineColor(kBlue-3);
  globalUncBox_15001 -> SetFillColorAlpha(kBlue-3,0.6);
  globalUncBox_15001 -> SetLineWidth(1);
  //globalUncBox_15001 -> Draw("l same");
  
  
  CMS_lumi_raaCent( c1, iPeriod, iPos );


	c1->Update();
  c1->SaveAs("plots/comp15001_RAA_int_asym_newstyle_Constrain.pdf");
  c1->SaveAs("plots/comp15001_RAA_int_asym_newstyle_Constrain.png");
  c1->SaveAs("Figure7.C");

/*
	///////////////////////////////////////////////////////////////////
	//// save as a root file
	TFile *outFile = new TFile("RAA_vs_pt.root", "RECREATE");
	outFile->cd();
	for (int is=0; is<nfile; is++){
		gRAA_sys[is]->Write();	
		gRAA[is]->Write();	
	}
	outFile->Close();
*/	
	return;

} // end of main func.

