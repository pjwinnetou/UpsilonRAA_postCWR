#include "SONGKYO.h"
#include "tdrstyle.C"
#include "CMS_lumi_raaCent.C"
#include "../cutsAndBin.h"
#include "../commonUtility.h"

void strickland_RAA_pt_isArrow_asym_postCWRConstrain(bool isArrow=true)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text
  int iPeriod = 101; // 1: pp, 2: pPb, 3: PbPb, 100: RAA vs cent, 101: RAA vs pt or rap
  int iPos = 33;
  
  const int nState = 3; // Y(1S), Y(2S), and Y(3S)
  double xmax = 30.0;
//  double relsys = 0.1;

  double exsys_1s[6] =  {1., 1., 1., 1.5, 1.5, 9.};
  double exsys_2s[3] =  {2., 2.5, 10.5};
  double exsys_3s[2] =  {3., 12.};

  ////////////////////////////////////////////////////////////////
  //// read input file : value & stat.
  TFile* fIn[nState];
	TGraphErrors* gRAA[nState];
	TGraphErrors* gRAA_sys[nState];
  for (int is=0; is<nState; is++){
  	fIn[is] = new TFile(Form("Ups_%d_RAA.root",is+1),"READ");
    gRAA[is]=(TGraphErrors*)fIn[is]->Get("gRAA_pt");
    gRAA_sys[is]=(TGraphErrors*)fIn[is]->Get("gRAA_pt");
    //cout << "gRAA["<<is<<"] = " <<gRAA[is] << endl;
  }
  //// read input file : syst.
  TFile* fInSys[nState];
  TH1D* hSys[nState];
  int npoint[nState];
  for (int is=0; is<nState; is++){
  	fInSys[is] = new TFile(Form("../Systematic/postCWR_mergedSys_constrain_ups%ds_asymHi.root",is+1),"READ");
    hSys[is]=(TH1D*)fInSys[is]->Get("hptRAA_merged");
    npoint[is] = hSys[is]->GetSize()-2;
    cout << "*** Y("<<is+1<<") : # of point = " << npoint[is] << endl;
  } 
  
  //// set bin width and calculate systematic uncertainties
  double pxtmp, pytmp, extmp, eytmp;
  double relsys;

  for (int is=0; is<nState; is++){
    cout << is+1 <<"th state***************" << endl;
    if (npoint[is] != gRAA[is]->GetN()) {cout << "Error!! data file and syst. file have different binnig!" << endl; return; }
    for (int ipt=0; ipt< npoint[is] ; ipt++) { //bin by bin
      pxtmp=0; pytmp=0; extmp=0; eytmp=0; relsys=0;
      gRAA[is]->GetPoint(ipt, pxtmp, pytmp); 
      extmp=gRAA[is]->GetErrorX(ipt);
      eytmp=gRAA[is]->GetErrorY(ipt);
      relsys=hSys[is]->GetBinContent(ipt+1);
      cout << ipt <<"th bin RAA value = " << pytmp << endl;
      cout << ipt <<"th bin stat. = " << eytmp << endl;
      //cout << ipt <<"th bin rel. syst. = " << relsys << endl;
      cout << ipt <<"th bin syst. = " << pytmp*relsys << endl; 
      //// 1) remove ex from gRAA
      gRAA[is]->SetPointError(ipt, 0, eytmp);
      //// 2) set ey for gRAA_sys
      //gRAA_sys[is]->SetPointError(ipt, extmp, pytmp*relsys);
      if (is==0) gRAA_sys[is]->SetPointError(ipt, exsys_1s[ipt], pytmp*relsys);
      else if (is==1) gRAA_sys[is]->SetPointError(ipt, exsys_2s[ipt], pytmp*relsys);
      else gRAA_sys[is]->SetPointError(ipt, exsys_3s[ipt], pytmp*relsys);
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


  if (isArrow == true){
        gRAA_sys[2]->SetPoint(0,-10,-10);
        gRAA_sys[2]->SetPointError(0,0,0);
        gRAA_sys[2]->SetPoint(1,-11,-11);
        gRAA_sys[2]->SetPointError(1,0,0);
        gRAA_sys[2]->SetPoint(2,-12,-12);
        gRAA_sys[2]->SetPointError(2,0,0);
        gRAA[2]->SetPoint(0,-10,-10);
        gRAA[2]->SetPointError(0,0,0);
        gRAA[2]->SetPoint(1,-11,-11);
        gRAA[2]->SetPointError(1,0,0);
        gRAA[2]->SetPoint(2,-12,-12);
        gRAA[2]->SetPointError(2,0,0);
        gRAA_sys[2]->GetHistogram()->GetXaxis()->SetLimits(0,30);
        gRAA_sys[2]->GetHistogram()->GetXaxis()->SetRangeUser(0,30);
        gRAA_sys[2]->SetMinimum(0.0);
        gRAA_sys[2]->SetMaximum(1.3);
        gRAA[2]->GetHistogram()->GetXaxis()->SetRangeUser(0,30);
        gRAA[2]->GetHistogram()->GetXaxis()->SetLimits(0,30);
        gRAA[2]->SetMinimum(0.0);
        gRAA[2]->SetMaximum(1.3);
      }
  
  //// draw  
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.067);
  for (int is=0; is<nState; is++){
    if ( is==0) {gRAA_sys[is]->Draw("A5");}
    else if(is==ulstate && isArrow==true) {
      for(int ipt=0;ipt<n3s;ipt++){
        box68per[ipt]->Draw("l");
      }
      gRAA_sys[is]->Draw("5");
    }
    else {gRAA_sys[is]->Draw("5");}
  }
  for(int is=0;is<nState;is++){
    if(is==ulstate && isArrow==true) {
      for(int ipt=0;ipt<n3s;ipt++) {
        arr95per[ipt]->Draw();
      }
      gRAA[is]->Draw("P");
    }
    else {gRAA[is]->Draw("P");}
  }
  
  dashedLine(0.,1.,xmax,1.,1,1);
  TLegend *leg= new TLegend(0.64, 0.62, 0.855, 0.74);
  SetLegendStyle(leg);
  TLegend *leg_up= new TLegend(0.64, 0.50, 0.85, 0.62);
  SetLegendStyle(leg_up);

  TArrow *arrLeg = new TArrow(18.59,0.532,18.59,0.582,0.02,"<-|");
  arrLeg->SetLineColor(kGreen+2);
  arrLeg->SetLineWidth(2);

  if (isArrow==false) { 
    for (int is=0; is<nState; is++){
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
    ent=leg_up->AddEntry("ent"," #varUpsilon(3S) 95\% CL","f");
    ent->SetLineColor(kWhite);
//    leg_up->SetTextSize(0.03);
    leg->Draw("same");
    leg_up->Draw("same");
    arrLeg->Draw();
  }


  //// draw text
  double sz_init = 0.925; double sz_step = 0.0535;
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu} > 4 GeV/c");
//  globtex->DrawLatex(0.22, sz_init, "p_{T}^{#mu#mu} < 30 GeV/c");
  globtex->DrawLatex(0.22, sz_init-sz_step, "|y| < 2.4");
//  globtex->DrawLatex(0.22, sz_init-sz_step*2, "|#eta^{#mu}| < 2.4");
  globtex->DrawLatex(0.22, sz_init-sz_step*2, "Cent. 0-100%");
  
  TFile *fstrickland = new TFile("TheoryCurve/StrickLand_RAA_5023.root","READ");
  
  TGraphErrors *gRAA_1S_strickland[3]; 
  TGraphErrors *gRAA_2S_strickland[3]; 
  
  for(int i=0;i<3;i++)
  {
    gRAA_1S_strickland[i] = (TGraphErrors*) fstrickland-> Get(Form("RAA_strick_pt_1S_%d",i));
    gRAA_2S_strickland[i] = (TGraphErrors*) fstrickland-> Get(Form("RAA_strick_pt_2S_%d",i));
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
    gRAA_1S_strickland[i]->Draw("same");
    gRAA_2S_strickland[i]->Draw("same");
  }

  const int npt = 3;
  
  TGraphErrors* g3S = new TGraphErrors("TheoryCurve/Y3Spt5023Xi0.tsv","%lg %lg %lg %lg","\t");


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

  double x[npt] = {2.5,3.5,4.0};

  TGraphErrors* g1t = new TGraphErrors(g3S->GetN(),g3S->GetX(),y3S_1,x,0);
  TGraphErrors* g2t = new TGraphErrors(g3S->GetN(),g3S->GetX(),y3S_2,x,0);
  TGraphErrors* g3t = new TGraphErrors(g3S->GetN(),g3S->GetX(),y3S_3,x,0);



  g1t->Draw("samep");
  g1t->SetLineWidth(3);
  g1t->SetLineStyle(3);
  g1t->SetLineColor(kGreen+2);
  g1t->SetMarkerSize(0);
  g2t->Draw("SAMEp");
  g2t->SetLineWidth(3);
  g2t->SetLineStyle(1);
  g2t->SetLineColor(kGreen+2);
  g2t->SetMarkerSize(0);
  g3t->Draw("SAMEp");
  g3t->SetLineWidth(3);
  g3t->SetLineColor(kGreen+2);
  g3t->SetLineStyle(8);
  g3t->SetMarkerSize(0);


  TLegend *leg_strick= new TLegend(0.2, 0.516, 0.4, 0.646);
  SetLegendStyle(leg_strick);
  leg_strick->SetTextSize(0.040);
  leg_strick->AddEntry(gRAA_1S_strickland[2],"Y(1S)","l");
  leg_strick->AddEntry(gRAA_2S_strickland[2],"Y(2S)","l");
//  leg_strick->Draw("same");

  double line_y = 0.741;
  double line_y_diff = 0.07;
  double line_y_diff_in = 0.02;
  double line_x_end = 4.4;
  double line_x_start = 2.5;

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

/*

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
*/
  drawText2("4#pi #eta/s=1", line_x_end+1, line_y-0.0355, 22);
  drawText2("4#pi #eta/s=2", line_x_end+1, line_y-line_y_diff*1-0.0355 - line_y_diff_in, 22);
  drawText2("4#pi #eta/s=3", line_x_end+1, line_y-line_y_diff*2.04-0.0355 - line_y_diff_in*2, 22);

  drawText2("Krouppa, Strickland",line_x_start-0.1,line_y+0.06,22);





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
  
  CMS_lumi_raaCent( c1, iPeriod, iPos );

	c1->Update();
  c1->SaveAs(Form("plots/Strickland_RAA_vs_pt_isArrow%d_asym_postCWRConstrain.pdf",(int)isArrow));
  c1->SaveAs(Form("plots/Strickland_RAA_vs_pt_isArrow%d_asym_postCWRConstrain.png",(int)isArrow));

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

