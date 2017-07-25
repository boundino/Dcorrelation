#ifndef _XJJROOTUTI_H_
#define _XJJROOTUTI_H_

#include "rootheaders.h"

namespace xjjrootuti
{
  void setgstyle();
  void sethempty(TH2F* hempty, float xoffset=0, float yoffset=0);
  void setth1(TH1F* hempty, float xoffset=0, float yoffset=0);
  void setth1(TH1D* hempty, float xoffset=0, float yoffset=0);
  void settgr(TGraph* gr, float xoffset=0, float yoffset=0);
  void settgr(TGraphErrors* gr, float xoffset=0, float yoffset=0);
  void settgr(TGraphAsymmErrors* gr, float xoffset=0, float yoffset=0);
  void drawCMS(TString collision);
  void settex(TLatex* tex, Float_t size=0.04, Int_t align=12);
  void setleg(TLegend* leg, Float_t size=0.04);
  void setbranchaddress(TTree* nt, const char* bname, void* addr);
}

void xjjrootuti::setgstyle()
{
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetTextSize(0.05);
  gStyle->SetTextFont(42);
  gStyle->SetPadRightMargin(0.043);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.145);
  gStyle->SetTitleX(.0f);
}

void xjjrootuti::sethempty(TH2F* hempty, float xoffset/*=0*/, float yoffset/*=0*/)
{
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetXaxis()->SetTitleOffset(1.10+xoffset);
  hempty->GetYaxis()->SetTitleOffset(1.30+yoffset);
  hempty->GetXaxis()->SetTitleSize(0.05);
  hempty->GetYaxis()->SetTitleSize(0.05);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.045);
  hempty->GetYaxis()->SetLabelSize(0.045);
  hempty->SetStats(0);
}

void xjjrootuti::setth1(TH1F* hempty, float xoffset/*=0*/, float yoffset/*=0*/)
{
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetXaxis()->SetTitleOffset(1.10+xoffset);
  hempty->GetYaxis()->SetTitleOffset(1.30+yoffset);
  hempty->GetXaxis()->SetTitleSize(0.05);
  hempty->GetYaxis()->SetTitleSize(0.05);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.045);
  hempty->GetYaxis()->SetLabelSize(0.045);
  hempty->SetStats(0);
  hempty->SetMarkerStyle(20);
  hempty->SetMarkerSize(1.1);
  hempty->SetMarkerColor(kBlack);
  hempty->SetLineColor(kBlack);
}

void xjjrootuti::setth1(TH1D* hempty, float xoffset/*=0*/, float yoffset/*=0*/)
{
  hempty->GetXaxis()->CenterTitle();
  hempty->GetYaxis()->CenterTitle();
  hempty->GetXaxis()->SetTitleOffset(1.10+xoffset);
  hempty->GetYaxis()->SetTitleOffset(1.30+yoffset);
  hempty->GetXaxis()->SetTitleSize(0.05);
  hempty->GetYaxis()->SetTitleSize(0.05);
  hempty->GetXaxis()->SetTitleFont(42);
  hempty->GetYaxis()->SetTitleFont(42);
  hempty->GetXaxis()->SetLabelFont(42);
  hempty->GetYaxis()->SetLabelFont(42);
  hempty->GetXaxis()->SetLabelSize(0.045);
  hempty->GetYaxis()->SetLabelSize(0.045);
  hempty->SetStats(0);
  hempty->SetMarkerStyle(20);
  hempty->SetMarkerSize(1.1);
  hempty->SetMarkerColor(kBlack);
  hempty->SetLineColor(kBlack);
}

void xjjrootuti::settgr(TGraph* gr, float xoffset/*=0*/, float yoffset/*=0*/)
{
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleOffset(1.10+xoffset);
  gr->GetYaxis()->SetTitleOffset(1.30+yoffset);
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelSize(0.045);
  gr->GetYaxis()->SetLabelSize(0.045);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.1);
  gr->SetMarkerColor(kBlack);
  gr->SetLineColor(kBlack);
}

void xjjrootuti::settgr(TGraphErrors* gr, float xoffset/*=0*/, float yoffset/*=0*/)
{
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleOffset(1.10+xoffset);
  gr->GetYaxis()->SetTitleOffset(1.30+yoffset);
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelSize(0.045);
  gr->GetYaxis()->SetLabelSize(0.045);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.1);
  gr->SetMarkerColor(kBlack);
  gr->SetLineColor(kBlack);
}

void xjjrootuti::settgr(TGraphAsymmErrors* gr, float xoffset/*=0*/, float yoffset/*=0*/)
{
  gr->GetXaxis()->CenterTitle();
  gr->GetYaxis()->CenterTitle();
  gr->GetXaxis()->SetTitleOffset(1.10+xoffset);
  gr->GetYaxis()->SetTitleOffset(1.30+yoffset);
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTitleFont(42);
  gr->GetYaxis()->SetTitleFont(42);
  gr->GetXaxis()->SetLabelFont(42);
  gr->GetYaxis()->SetLabelFont(42);
  gr->GetXaxis()->SetLabelSize(0.045);
  gr->GetYaxis()->SetLabelSize(0.045);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.1);
  gr->SetMarkerColor(kBlack);
  gr->SetLineColor(kBlack);
}

void xjjrootuti::drawCMS(TString collision)
{
  TLatex* texCms = new TLatex(0.18,0.93, "#scale[1.25]{CMS} Preliminary");
  texCms->SetNDC();
  texCms->SetTextAlign(12);
  texCms->SetTextSize(0.04);
  texCms->SetTextFont(42);
  texCms->Draw();
  TLatex* texCol = new TLatex(0.96,0.93, Form("%s #sqrt{s_{NN}} = 5.02 TeV",collision.Data()));
  texCol->SetNDC();
  texCol->SetTextAlign(32);
  texCol->SetTextSize(0.04);
  texCol->SetTextFont(42);
  texCol->Draw();
}

void xjjrootuti::settex(TLatex* tex, Float_t size/*=0.04*/, Int_t align/*=12*/)
{
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetTextAlign(align);
  tex->SetTextSize(size);
}

void xjjrootuti::setleg(TLegend* leg, Float_t size/*=0.04*/)
{
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(size);
}

void xjjrootuti::setbranchaddress(TTree* nt, const char* bname, void* addr)
{
  nt->SetBranchStatus(bname, 1);
  nt->SetBranchAddress(bname, addr);
}

/*
void removeError(TH1F* h)
{
  for(int i=1;i<=h->GetNbinsX();i++)
    {
      h->SetBinError(i,0);
    }	
}

// divide by bin width
void divideBinWidth(TH1* h)
{
  //h->Sumw2();
  for(int i=1;i<=h->GetNbinsX();i++)
    {
      Float_t val = h->GetBinContent(i);
      Float_t valErr = h->GetBinError(i);
      val/=h->GetBinWidth(i);
      valErr/=h->GetBinWidth(i);
      h->SetBinContent(i,val);
      h->SetBinError(i,valErr);
    }
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

// make a histogram from TF1 function
TH1F* functionHist(TF1* f, TH1F* h, TString fHistname)
{
  TH1F* hF = (TH1F*)h->Clone(fHistname);
  for (int i=1;i<=h->GetNbinsX();i++)
    {
      Double_t var = f->Integral(h->GetBinLowEdge(i),h->GetBinLowEdge(i+1))/h->GetBinWidth(i);
      hF->SetBinContent(i,var);
      hF->SetBinError(i,0);
    }
  return hF;
}

TLegend* myLegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
  TLegend* leg = new TLegend(x1,y1,x2,y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  return leg; 
}
*/


#endif
