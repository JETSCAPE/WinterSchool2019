#include"Custom.h"

// define a function of calculating the ratio between two root graphs

TH1D * GetRatioAsHist(string Name, TGraphErrors* Numerator, TGraphErrors* Denominator)
{
  int N1=Numerator->GetN(); double *X1=Numerator->GetX();double *Y1=Numerator->GetY();double *EY1=Numerator->GetEY(); double *EX1=Numerator->GetEX();
  int N2=Denominator->GetN(); double *X2=Denominator->GetX();double *Y2=Denominator->GetY();double *EY2=Denominator->GetEY(); double *EX2=Denominator->GetEX();

  double * XBinArray = new double [N1+1]; // boundaries of bins of the new histogram
  for(int i=0;i<N1;i++)
    { XBinArray[i] = X1[i] - EX1[i];
    }
  XBinArray[N1] = X1[N1-1] + EX1[N1-1];
  TH1D *G = new TH1D(Name.c_str(), Name.c_str(), N1, XBinArray);
  if(N1==N2)
    {   for(int i=0;i<N1;i++)
	{ 
	  double Ratio = Y1[i]/Y2[i];
	  double RatioError = Ratio*TMath::Sqrt( pow(EY1[i]/Y1[i],2.0) + pow(EY2[i]/Y2[i],2.0) );
	  G->SetBinContent(i+1,Ratio);
	  G->SetBinError(i+1,RatioError);
	  cout<<"XBinNum = "<<i+1<<"\t"<<"X[i] ="<<X1[i]<<"\t"<<Ratio<<endl;
	}
    }
  else
    {cout<<"Bins in x-axis are not the same "<<endl;
    }
  return G;
}

void RatioPPVsPbPb()
{
  CustomGlobal();
  char FileName[5000];
  
  TFile *SpectraPP  = new TFile("CombinedSpectra-pp.root");
  TFile *SpectraPbPb  = new TFile("CombinedSpectra-PbPb.root");


  TGraphErrors *JetCrossSectionPP          = (TGraphErrors *)SpectraPP->FindObjectAny("TotalDifferentialJetCrossSection")->Clone();
  TGraphErrors *SingleHadronYieldPP        = (TGraphErrors *)SpectraPP->FindObjectAny("TotalDifferentialSingleHadronYield")->Clone();
  TGraphErrors *JetShapePP                 = (TGraphErrors *)SpectraPP->FindObjectAny("TotalJetShape")->Clone();
  TGraphErrors *JetFFPP = (TGraphErrors *)SpectraPP->FindObjectAny("TotalJetFragmentationFunction")->Clone();
  TGraphErrors *JetMassPP                  = (TGraphErrors *)SpectraPP->FindObjectAny("TotalJetMass")->Clone();

  TGraphErrors *JetCrossSectionPbPb          = (TGraphErrors *)SpectraPbPb->FindObjectAny("TotalDifferentialJetCrossSection")->Clone();
  TGraphErrors *SingleHadronYieldPbPb        = (TGraphErrors *)SpectraPbPb->FindObjectAny("TotalDifferentialSingleHadronYield")->Clone();
  TGraphErrors *JetShapePbPb                 = (TGraphErrors *)SpectraPbPb->FindObjectAny("TotalJetShape")->Clone();
  TGraphErrors *JetFFPbPb = (TGraphErrors *)SpectraPbPb->FindObjectAny("TotalJetFragmentationFunction")->Clone();
  TGraphErrors *JetMassPbPb                  = (TGraphErrors *)SpectraPbPb->FindObjectAny("TotalJetMass")->Clone();

  TH1D *MCJetRAA         = GetRatioAsHist("MCJetRAA", JetCrossSectionPbPb, JetCrossSectionPP);
  TH1D *MCSingleHadronRAA= GetRatioAsHist("MCSingleHadronRAA", SingleHadronYieldPbPb, SingleHadronYieldPP );
  TH1D *MCJetShapeRatio  = GetRatioAsHist("MCJetShapeRatio", JetShapePbPb, JetShapePP);
  TH1D *MCJetFFRatio     = GetRatioAsHist("MCJetFFRatio", JetFFPbPb, JetFFPP);
  TH1D *MCJetMassRatio   = GetRatioAsHist("MCJetMassRatio", JetMassPbPb, JetMassPP);

  TCanvas *canvas = new TCanvas("canvas","canvas",800,700);
  canvas->Divide(2,3);
  canvas->cd(1);   
  MCJetRAA->Draw("AXIS"); Custom(MCJetRAA,2);
  //Sys->Draw("p2 same");Custom(Sys,kBlack);Sys->SetFillColorAlpha(kGray,0.01);Sys->SetFillStyle(1001);
  //Stat->Draw("P same");Custom(Stat,kBlack); Stat->SetLineWidth(2);
  MCJetRAA->Draw("Hist E X0 ][  same");MCJetRAA->Draw("Hist ][ same");
  MCJetRAA->GetXaxis()->SetRangeUser(50,300);MCJetRAA->GetYaxis()->SetRangeUser(0.0, 1.0);
  MCJetRAA->GetXaxis()->SetTitle("Jet #it{p}_{T}  (GeV)");MCJetRAA->GetYaxis()->SetTitle("Jet R_{AA}");
  gPad->SetMargin(0.12,0.03,0.15,0.2);

  
  TLegend *leg1 = new TLegend(0.16,0.16,0.52 ,0.37);
  leg1->SetTextFont(22);leg1->SetTextSize(0.06);
  //leg1->AddEntry(Sys,"CMS 2.76 TeV","lepf"); 
  leg1->AddEntry(MCJetRAA, "JETSCAPE","l");
  leg1->SetFillStyle(0);  leg1->Draw();
  TPaveText *Pt1 = new TPaveText(50,1.05,300,1.2,"NB");
  Pt1->AddText("Jet R_{AA} at 2.76 TeV, R=0.3, |#eta_{jet}|<2.0");
  Pt1->SetTextFont(22);Pt1->SetTextSize(0.06);Pt1->Draw();

  canvas->cd(2);
  MCSingleHadronRAA->Draw("AXIS"); Custom(MCSingleHadronRAA,2);
  //Sys->Draw("p2 same");Custom(Sys,kBlack);Sys->SetFillColorAlpha(kGray,0.01);Sys->SetFillStyle(1001);
  //Stat->Draw("P same");Custom(Stat,kBlack); Stat->SetLineWidth(2);
  MCSingleHadronRAA->Draw("Hist E X0 ][  same");MCSingleHadronRAA->Draw("Hist ][ same");
  MCSingleHadronRAA->GetXaxis()->SetRangeUser(0,100);MCSingleHadronRAA->GetYaxis()->SetRangeUser(0.0, 1.0);
  MCSingleHadronRAA->GetXaxis()->SetTitle("#it{p}_{T}  (GeV)");MCSingleHadronRAA->GetYaxis()->SetTitle(" Hadron R_{AA}");
  gPad->SetMargin(0.12,0.03,0.15,0.2);


  TLegend *leg2 = new TLegend(0.16,0.86,0.52 ,0.37);
  leg2->SetTextFont(22);leg2->SetTextSize(0.06);
  //leg2->AddEntry(G5,"CMS 2.76 TeV","lepf"); 
  leg2->AddEntry(MCSingleHadronRAA, "JETSCAPE","l");
  leg2->SetFillStyle(0);  leg2->Draw();
  TPaveText *Pt2 = new TPaveText(5,1.05,100,1.2,"NB");
  Pt2->AddText("Single hadron R_{AA} at 2.76 TeV,  |#eta|<1.0");
  Pt2->SetTextFont(22);Pt2->SetTextSize(0.06);Pt2->Draw();
  
  canvas->cd(3);
  MCJetShapeRatio->Draw("AXIS"); Custom(MCJetShapeRatio,2);
  //Sys->Draw("p2 same");Custom(Sys,kBlack);Sys->SetFillColorAlpha(kGray,0.01);Sys->SetFillStyle(1001);
  //Stat->Draw("P same");Custom(Stat,kBlack); Stat->SetLineWidth(2);
  MCJetShapeRatio->Draw("Hist E X0 ][  same");MCJetShapeRatio->Draw("Hist ][ same");
  MCJetShapeRatio->GetXaxis()->SetRangeUser(0,0.35);MCJetShapeRatio->GetYaxis()->SetRangeUser(0.5, 1.3);
  MCJetShapeRatio->GetXaxis()->SetTitle("r");MCJetShapeRatio->GetYaxis()->SetTitle("#rho_{PbPb}/#rho_{pp}");
  gPad->SetMargin(0.12,0.03,0.15,0.2);


  TLegend *leg3 = new TLegend(0.16,0.16,0.52 ,0.37);
  leg3->SetTextFont(22);leg3->SetTextSize(0.06);
  //leg3->AddEntry(G5,"CMS 2.76 TeV","lepf"); 
  leg3->AddEntry(MCJetShapeRatio, "JETSCAPE","l");
  leg3->SetFillStyle(0);  leg3->Draw();
  TPaveText *Pt3 = new TPaveText(0, 1.32, 0.35, 1.42,"NB");
  Pt3->AddText("JetShape ratio PbPb/pp at 2.76 TeV, p^{jet}_{T}>100 GeV, |#eta|<2.0");
  Pt3->SetTextFont(22);Pt3->SetTextSize(0.06);Pt3->Draw();

  canvas->cd(4);
  MCJetFFRatio->Draw("AXIS"); Custom(MCJetFFRatio,2);
  //Sys->Draw("p2 same");Custom(Sys,kBlack);Sys->SetFillColorAlpha(kGray,0.01);Sys->SetFillStyle(1001);
  //Stat->Draw("P same");Custom(Stat,kBlack); Stat->SetLineWidth(2);
  MCJetFFRatio->Draw("Hist E X0 ][  same");MCJetFFRatio->Draw("Hist ][ same");
  MCJetFFRatio->GetXaxis()->SetRangeUser(0.0,1.0);MCJetFFRatio->GetYaxis()->SetRangeUser(0.0, 1.5);
  MCJetFFRatio->GetXaxis()->SetTitle("Z");MCJetFFRatio->GetYaxis()->SetTitle(" Jet FF ratio PbPb/pp ");
  gPad->SetMargin(0.12,0.03,0.15,0.2);


  TLegend *leg4 = new TLegend(0.16,0.16,0.52 ,0.37);
  leg4->SetTextFont(22);leg4->SetTextSize(0.06);
  //leg4->AddEntry(G5,"CMS 2.76 TeV","lepf"); 
  leg4->AddEntry(MCJetFFRatio, "JETSCAPE","l");
  leg4->SetFillStyle(0);  leg4->Draw();
  TPaveText *Pt4 = new TPaveText(0,1.52,1,1.72,"NB");
  Pt4->AddText("Jet Fragmentation function ratio PbPb/pp at 2.76TeV, |#eta|<2.0");
  Pt4->SetTextFont(22);Pt4->SetTextSize(0.06);Pt4->Draw();

  canvas->cd(5);
  JetMassPbPb->Draw("AP"); Custom(JetMassPbPb,2);
  JetMassPP->Draw("P same"); Custom(JetMassPP,3);
  JetMassPbPb->GetXaxis()->SetRangeUser(0,30);JetMassPbPb->GetYaxis()->SetRangeUser(-0.05, 0.15);
  JetMassPbPb->GetXaxis()->SetTitle("Jet Mass (GeV)");JetMassPbPb->GetYaxis()->SetTitle("1/N dN/dM_{jet} (GeV^{-1})");
 //MCJetMassRatio->Draw("AXIS"); Custom(MCJetMassRatio,2);
  //Sys->Draw("p2 same");Custom(Sys,kBlack);Sys->SetFillColorAlpha(kGray,0.01);Sys->SetFillStyle(1001);
  //Stat->Draw("P same");Custom(Stat,kBlack); Stat->SetLineWidth(2);
  //MCJetMassRatio->Draw("Hist E X0 ][  same");MCJetMassRatio->Draw("Hist ][ same");
  //MCJetMassRatio->GetXaxis()->SetRangeUser(50,300);MCJetMassRatio->GetYaxis()->SetRangeUser(0.0, 1.5);
  //MCJetMassRatio->GetXaxis()->SetTitle("Jet Mass");MCJetMassRatio->GetYaxis()->SetTitle("Jet R_{AA}");
  gPad->SetMargin(0.12,0.03,0.15,0.2);


  TLegend *leg5 = new TLegend(0.36,0.16,0.72 ,0.37);
  leg5->SetTextFont(22);leg5->SetTextSize(0.06);
  //leg5->AddEntry(G5,"CMS 2.76 TeV","lepf"); 
  leg5->AddEntry(JetMassPP, "JETSCAPE pp","l");
  leg5->AddEntry(JetMassPbPb, "JETSCAPE PbPb","l");
  leg5->SetFillStyle(0);  leg5->Draw();
  TPaveText *Pt5 = new TPaveText(0,0.155,30,0.19,"NB");
  Pt5->AddText("JetMass at 2.76 TeV, |#eta|<2.0");
  Pt5->SetTextFont(22);Pt5->SetTextSize(0.06);Pt5->Draw();
 
  sprintf(FileName,"PlotRatioJETSCAPE.png");
  //canvas->SaveAs(FileName);

    
  

}
