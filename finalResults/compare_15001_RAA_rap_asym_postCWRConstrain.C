#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"

void compare_15001_RAA_rap_asym_postCWRConstrain(int istate=1) //1 or 2 (1S or 2S)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 101; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nfile = 2; // 0: 15001, 1: ours
  double xmax = 2.4;
//  double relsys = 0.1;
  
  double exsys_1s[6] =  {0.2, 0.2, 0.2, 0.2, 0.2, 0.2};
  double exsys_2s[3] =  {0.4, 0.4, 0.4};
  double exsys_3s[2] =  {0.6, 0.6};

  //// 15001 values
  const int cn_1s =  6;
  double cpx_1s[cn_1s] =  {0.2, 0.6, 1.0, 1.4, 1.8, 2.2};
  double cpy_1s[cn_1s] =  {0.446, 0.415, 0.491, 0.490, 0.479, 0.399};
  double cex_1s[cn_1s] =  {0., 0., 0., 0., 0., 0.};
  double cey_1s[cn_1s] =  {0.028, 0.028, 0.033, 0.036, 0.038, 0.060};
  double cexsys_1s[cn_1s] =  {0.2, 0.2, 0.2, 0.2, 0.2, 0.2};
  double ceysys_1s_1[cn_1s] =  {0.043, 0.069, 0.065, 0.072, 0.072, 0.084};
  double ceysys_1s_2[cn_1s] =  {0.035, 0.033, 0.039, 0.038, 0.038, 0.031};
  double ceysys_1s[cn_1s] = {0.043, 0.069, 0.065, 0.072, 0.072, 0.084};
//  for (int it=0; it < cn_1s ; it++) {
//    ceysys_1s[it] = TMath::Sqrt(ceysys_1s_1[it]*ceysys_1s_1[it]+ceysys_1s_2[it]*ceysys_1s_2[it]);
//  }

  const int cn_2s =  2;
  double cpx_2s[cn_2s] =  {0.6, 1.8};
  double cpy_2s[cn_2s] =  {0.116, 0.138};
  double cex_2s[cn_2s] =  {0., 0.};
  double cey_2s[cn_2s] =  {0.034, 0.049};
  double cexsys_2s[cn_2s] =  {0.6, 0.6};
  double ceysys_2s_1[cn_2s] =  {0.036, 0.073};
  double ceysys_2s_2[cn_2s] =  {0.002, 0.002};
  double ceysys_2s[cn_2s] = {0.039, 0.064};
//  for (int it=0; it < cn_2s ; it++) {
//    ceysys_2s[it] = TMath::Sqrt(ceysys_2s_1[it]*ceysys_2s_1[it]+ceysys_2s_2[it]*ceysys_2s_2[it]);
//  }

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
	TGraphErrors* gRAA[nfile];
	TGraphErrors* gRAA_sys[nfile];
  //// 1) 15001
  if (istate==1) {
    gRAA[0] = new TGraphErrors(cn_1s, cpx_1s, cpy_1s, cex_1s, cey_1s); 
    gRAA_sys[0] = new TGraphErrors(cn_1s, cpx_1s, cpy_1s, cexsys_1s, ceysys_1s); 
  }
  else {
    gRAA[0] = new TGraphErrors(cn_2s, cpx_2s, cpy_2s, cex_2s, cey_2s); 
    gRAA_sys[0] = new TGraphErrors(cn_2s, cpx_2s, cpy_2s, cexsys_2s, ceysys_2s); 
  } 
  //// 2) ours
  TFile* fIn = new TFile(Form("Ups_%d_RAA.root",istate),"READ");
  gRAA[1]=(TGraphErrors*)fIn->Get("gRAA_rap");
  gRAA_sys[1]=(TGraphErrors*)fIn->Get("gRAA_rap");
  //// read input file : syst. 
  TFile* fInSys = new TFile(Form("../Systematic/postCWR_mergedSys_constrain_ups%ds_asymLo.root",istate),"READ");
  TH1D* hSys = (TH1D*)fInSys->Get("hrapRAA_merged");
  int npoint = hSys->GetSize()-2;
  cout << "*** Y("<<istate<<") : # of point = " << npoint << endl;
  
  //// set bin width and calculate systematic uncertainties 
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;
  if (npoint != gRAA[1]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
  for (int ipt=0; ipt< npoint; ipt++) {
    pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
    gRAA[1]->GetPoint(ipt, pxtmp, pytmp);
    extmp=gRAA[1]->GetErrorX(ipt);
    eytmp=gRAA[1]->GetErrorY(ipt);
    relsys=hSys->GetBinContent(ipt+1);
    // 1) remove ex from gRAA
    gRAA[1]->SetPointError(ipt, 0, eytmp);
    // 2) set ey for gRAA_sys (assign 10% temporarily)
    //gRAA_sys[1]->SetPointError(ipt, extmp, pytmp*relsys);
    if (istate==1) gRAA_sys[1]->SetPointError(ipt, exsys_1s[ipt], pytmp*relsys);
    else gRAA_sys[1]->SetPointError(ipt, exsys_2s[ipt], pytmp*relsys);
  }
  
  ////////////////////////////////////////////////////////////////
  
  //// graph style 
  if(istate==1){
    SetGraphStyleOpen(gRAA[0], 4, 4, 0);
    SetGraphStyleSys(gRAA_sys[0], 4); 
    SetGraphStyle(gRAA[1], 0, 0); 
    SetGraphStyleSys(gRAA_sys[1], 0); 
  }
  else if(istate==2){
    SetGraphStyleOpen(gRAA[0], 4, 4, 1);
    SetGraphStyleSys(gRAA_sys[0], 4); 
    SetGraphStyle(gRAA[1], 1, 1); 
    SetGraphStyleSys(gRAA_sys[1], 1); 
  }

  
  //// latex for text
  TLatex* globtex = new TLatex();
  globtex->SetNDC();
  globtex->SetTextAlign(12); //left-center
  globtex->SetTextFont(42);
  globtex->SetTextSize(0.040);
  
  //// legend
  TLegend *leg= new TLegend(0.55, 0.57, 0.95, 0.74);
  SetLegendStyle(leg);
  leg -> SetHeader(Form("#varUpsilon(%dS)",istate));
  leg -> AddEntry(gRAA[0],"#sqrt{s_{NN}} = 2.76 TeV","lp");
  leg -> AddEntry(gRAA[1],"#sqrt{s_{NN}} = 5.02 TeV","lp");

  TLegendEntry *header = (TLegendEntry*)leg->GetListOfPrimitives()->First();
  header->SetTextSize(0.046);
  header->SetTextFont(62);
  //// axis et. al
  gRAA_sys[0]->GetXaxis()->SetTitle("|y|");
  gRAA_sys[0]->GetXaxis()->CenterTitle();
  gRAA_sys[0]->GetYaxis()->SetTitle("R_{AA}");
  gRAA_sys[0]->GetYaxis()->CenterTitle();
  gRAA_sys[0]->GetXaxis()->SetLimits(0.,xmax);
  gRAA_sys[0]->SetMinimum(0.0);
  gRAA_sys[0]->SetMaximum(1.3);
  /// for rap
  gRAA_sys[0]->GetXaxis()->SetNdivisions(505);
 
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  for (int is=0; is<nfile; is++){
    if ( is==0) gRAA_sys[is]->Draw("A5");
    else gRAA_sys[is]->Draw("5");
    gRAA[is]->Draw("P");
	}
  dashedLine(0.,1.,xmax,1.,1,1);
  leg->Draw();

  //// draw text
  double sz_init = 0.927; double sz_step = 0.0535;
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
  globtex->DrawLatex(0.22, sz_init-sz_step, "p_{T} < 30 GeV");
//  globtex->DrawLatex(0.22, sz_init-sz_step, "|y|^{#mu#mu} < 2.4");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 2.4");
  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 0-100%");
  
  double TAA_unc_Global_Hi = 0.028;
  double TAA_unc_Global_Lo = 0.034;

  double sys_global_val_Hi = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+TAA_unc_Global_Hi*TAA_unc_Global_Hi+nMB_unc*nMB_unc);
  double sys_global_val_Lo = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+TAA_unc_Global_Lo*TAA_unc_Global_Lo+nMB_unc*nMB_unc);
  double sys_global_y_Hi = sys_global_val_Hi;
  double sys_global_y_Lo = sys_global_val_Lo;
  double sys_global_x = 0.06;
  //double sys_global_val = TMath::Sqrt(lumi_unc_pp*lumi_unc_pp+lumi_unc_aa*lumi_unc_aa);
  double sys_global_y_15001 = 0.075;
  TBox *globalUncBox = new TBox(xmax-sys_global_x*2,1-sys_global_y_Lo,xmax-sys_global_x-0.001,1+sys_global_y_Hi);
  if(istate==1){
    globalUncBox -> SetLineColor(kRed-2);
    globalUncBox -> SetFillColorAlpha(kPink-6,0.6);
  }
  else if(istate==2){
    globalUncBox -> SetLineColor(kBlue-3);
    globalUncBox -> SetFillColorAlpha(kBlue-3,0.6);
  }
  globalUncBox -> SetLineWidth(1);
  globalUncBox -> Draw("l same");
  
  TBox *globalUncBox_15001 = new TBox(xmax-sys_global_x,1-sys_global_y_15001,xmax,1+sys_global_y_15001);
  globalUncBox_15001 -> SetLineColor(kBlack);
  globalUncBox_15001 -> SetFillColorAlpha(kGray+3,0.6);
  globalUncBox_15001 -> SetLineWidth(1);
  globalUncBox_15001 -> Draw("l same");
  
  
  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs(Form("plots/%dS_comp15001_RAA_vs_rap_asym_Constrain.pdf",istate));
  c1->SaveAs(Form("plots/%dS_comp15001_RAA_vs_rap_asym_Constrain.png",istate));

/*
	///////////////////////////////////////////////////////////////////
	//// save as a root file
	TFile *outFile = new TFile("RAA_vs_rap.root", "RECREATE");
	outFile->cd();
	for (int is=0; is<nfile; is++){
		gRAA_sys[is]->Write();	
		gRAA[is]->Write();	
	}
	outFile->Close();
*/	
	return;

} // end of main func.

