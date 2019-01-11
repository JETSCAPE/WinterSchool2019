void DrawJetscapeLogo(double x1, double y1, double x2, double y2)
{

  // TPad * pad1 = new TPad("pad1","pad1",x1,y1,x2,y2);pad1->Draw();pad1->cd();
  //TImage *img = TImage::Open("/home/amit/Dropbox/JetscapeLogo.jpg");
  //img->Draw();
}

void Custom(TH1D *h1)
{
gPad->SetTickx();
gPad->SetTicky();
// h1->SetLineColor(color);
//h1->GetYaxis()->SetTitle("dN/dE per 1GeV");h1->GetXaxis()->SetTitle("E (GeV)");
h1->GetYaxis()->CenterTitle();
h1->GetXaxis()->CenterTitle();
//h1->GetYaxis()->SetRangeUser(0.001,10.0);
//h1->GetXaxis()->SetRangeUser(0.0,60);
h1->GetXaxis()->SetNdivisions(307);
h1->GetYaxis()->SetNdivisions(307);
h1->GetYaxis()->SetTitleSize(0.05);h1->GetYaxis()->SetTitleFont(22);
h1->GetXaxis()->SetTitleSize(0.05);h1->GetXaxis()->SetTitleFont(22);
h1->GetYaxis()->SetLabelSize(0.05);h1->GetYaxis()->SetLabelFont(22);
h1->GetXaxis()->SetLabelSize(0.05);h1->GetXaxis()->SetLabelFont(22);
h1->GetYaxis()->SetTitleOffset(1.2);
h1->GetXaxis()->SetTitleOffset(1.2);
h1->GetYaxis()->SetTickSize(0.04);
h1->GetXaxis()->SetTickSize(0.04);
h1->SetMarkerSize(3);h1->SetLineWidth(3);

}

void Custom(TH1D *h1, int color)
{
gPad->SetTickx();
gPad->SetTicky();
h1->SetLineColor(color);
h1->SetMarkerColor(color);
//h1->GetYaxis()->SetTitle("dN/dE per 1GeV");h1->GetXaxis()->SetTitle("E (GeV)");
h1->GetYaxis()->CenterTitle();
h1->GetXaxis()->CenterTitle();
//h1->GetYaxis()->SetRangeUser(0.001,10.0);
//h1->GetXaxis()->SetRangeUser(0.0,60);
h1->GetXaxis()->SetNdivisions(307);
h1->GetYaxis()->SetNdivisions(307);
h1->GetYaxis()->SetTitleSize(0.07);h1->GetYaxis()->SetTitleFont(22);
h1->GetXaxis()->SetTitleSize(0.07);h1->GetXaxis()->SetTitleFont(22);
h1->GetYaxis()->SetLabelSize(0.07);h1->GetYaxis()->SetLabelFont(22);
h1->GetXaxis()->SetLabelSize(0.07);h1->GetXaxis()->SetLabelFont(22);
h1->GetYaxis()->SetTitleOffset(0.8);
h1->GetXaxis()->SetTitleOffset(1.);
h1->GetYaxis()->SetTickSize(0.04);
h1->GetXaxis()->SetTickSize(0.04);
h1->SetMarkerSize(3);h1->SetLineWidth(3);

}


void CustomGlobal()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetLineWidth(3);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetEndErrorSize(2);
}

void Custom(TH1D & h, int color)
{
  gPad->SetTickx();
  gPad->SetTicky();
  h.SetMarkerColor(color);
  h.SetLineColor(color);
  //h.GetYaxis()->SetTitle("dN/dE per 1GeV");h.GetXaxis()->SetTitle("E (GeV)");
  h.GetYaxis()->CenterTitle();
  h.GetXaxis()->CenterTitle();
  h.GetYaxis()->SetRangeUser(0.001,10.0);
  h.GetXaxis()->SetRangeUser(0.0,60);
  h.GetXaxis()->SetNdivisions(307);
  h.GetYaxis()->SetNdivisions(307);
  h.GetYaxis()->SetTitleSize(0.05);h.GetYaxis()->SetTitleFont(22);
  h.GetXaxis()->SetTitleSize(0.05);h.GetXaxis()->SetTitleFont(22);
  h.GetYaxis()->SetLabelSize(0.05);h.GetYaxis()->SetLabelFont(22);
  h.GetXaxis()->SetLabelSize(0.05);h.GetXaxis()->SetLabelFont(22);
  h.GetYaxis()->SetTitleOffset(1.5);
  h.GetXaxis()->SetTitleOffset(1.2);
  h.GetYaxis()->SetTickSize(0.04);
  h.GetXaxis()->SetTickSize(0.04);
  h.SetMarkerSize(3);h.SetLineWidth(3);

}

void Custom(TGraph *h1)
{
  gPad->SetTickx();
  gPad->SetTicky();
  h1->SetLineColor(kBlack);  
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetTitleSize(0.05);h1->GetYaxis()->SetTitleFont(22);
  h1->GetXaxis()->SetTitleSize(0.05);h1->GetXaxis()->SetTitleFont(22);
  h1->GetYaxis()->SetLabelSize(0.05);h1->GetYaxis()->SetLabelFont(22);
  h1->GetXaxis()->SetLabelSize(0.05);h1->GetXaxis()->SetLabelFont(22);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTickSize(0.04);
  h1->GetXaxis()->SetTickSize(0.04);
  h1->SetMarkerStyle(8);
  h1->SetMarkerSize(1);h1->SetLineWidth(3);
  h1->SetFillColor(0);
}
void Custom(TGraph *h1, int colornum)
{
  gPad->SetTickx();
  gPad->SetTicky();
  h1->SetLineColor(colornum);
  h1->SetMarkerColor(colornum);
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetTitleSize(0.05);h1->GetYaxis()->SetTitleFont(22);
  h1->GetXaxis()->SetTitleSize(0.05);h1->GetXaxis()->SetTitleFont(22);
  h1->GetYaxis()->SetLabelSize(0.05);h1->GetYaxis()->SetLabelFont(22);
  h1->GetXaxis()->SetLabelSize(0.05);h1->GetXaxis()->SetLabelFont(22);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTickSize(0.04);
  h1->GetXaxis()->SetTickSize(0.04);
  h1->SetMarkerStyle(8);
  h1->SetMarkerSize(1);h1->SetLineWidth(3);
  h1->SetFillColor(0);
}


void Custom(TGraphErrors *h1)
{
  gPad->SetTickx();
  gPad->SetTicky();
  h1->SetLineColor(kRed);
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetTitleSize(0.05);h1->GetYaxis()->SetTitleFont(22);
  h1->GetXaxis()->SetTitleSize(0.05);h1->GetXaxis()->SetTitleFont(22);
  h1->GetYaxis()->SetLabelSize(0.05);h1->GetYaxis()->SetLabelFont(22);
  h1->GetXaxis()->SetLabelSize(0.05);h1->GetXaxis()->SetLabelFont(22);
  h1->GetYaxis()->SetTitleOffset(1.2);
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTickSize(0.04);
  h1->GetXaxis()->SetTickSize(0.04);
  h1->SetMarkerStyle(8);
  h1->SetMarkerSize(1);h1->SetLineWidth(3);

}

void Custom(TGraphErrors *h1, int colornum)
{
  gPad->SetTickx();
  gPad->SetTicky();
  h1->SetLineColor(colornum);
  h1->SetMarkerColor(colornum);
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->CenterTitle();
  h1->GetXaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetNdivisions(307);
  h1->GetYaxis()->SetTitleSize(0.07);h1->GetYaxis()->SetTitleFont(22);
  h1->GetXaxis()->SetTitleSize(0.07);h1->GetXaxis()->SetTitleFont(22);
  h1->GetYaxis()->SetLabelSize(0.07);h1->GetYaxis()->SetLabelFont(22);
  h1->GetXaxis()->SetLabelSize(0.07);h1->GetXaxis()->SetLabelFont(22);
  h1->GetYaxis()->SetTitleOffset(0.8);
  h1->GetXaxis()->SetTitleOffset(1.2);
  h1->GetYaxis()->SetTickSize(0.04);
  h1->GetXaxis()->SetTickSize(0.04);
  h1->SetMarkerStyle(8);
  h1->SetMarkerSize(1);h1->SetLineWidth(3);
  h1->SetFillColor(0);
}
