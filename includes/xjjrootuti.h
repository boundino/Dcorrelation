#ifndef _XJJROOTUTI_H_
#define _XJJROOTUTI_H_

#include "rootheaders.h"

namespace xjjrootuti
{
  const float x_tex_left_top = 0.23;
  const float y_tex_left_top = 0.91;
  const float dy_tex_left_top = 0.06;

  void setgstyle();
  void sethempty(TH2F* hempty, float xoffset=0, float yoffset=0);
  template <class T>
  void setthgr(T* hempty, float xoffset=0, float yoffset=0);
  template <class T>
  void setthgrstyle(T* h, Color_t mlcolor, Style_t mstyle=-1, Size_t msize=-1, Style_t lstyle=-1, Width_t lwidth=-1);
  void drawCMS(TString collision);
  void settex(TLatex* tex, Float_t size=0.04, Int_t align=12);
  void settexndraw(TLatex* tex, Float_t size=0.04, Int_t align=12);
  void setleg(TLegend* leg, Float_t size=0.04);
  void setlegndraw(TLegend* leg, Float_t size=0.04);
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

template <class T>
void xjjrootuti::setthgr(T* hempty, float xoffset/*=0*/, float yoffset/*=0*/)
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

template <class T>
void xjjrootuti::setthgrstyle(T* h, Color_t mlcolor, Style_t mstyle/*=-1*/, Size_t msize/*=-1*/, Style_t lstyle/*=-1*/, Width_t lwidth/*=-1*/)
{
  h->SetMarkerColor(mlcolor);
  h->SetLineColor(mlcolor);
  if(mstyle>=0) h->SetMarkerStyle(mstyle);
  if(lstyle>=0) h->SetLineStyle(lstyle);
  if(msize>=0)  h->SetMarkerSize(msize);
  if(lwidth>=0) h->SetLineWidth(lwidth);
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

void xjjrootuti::settexndraw(TLatex* tex, Float_t size/*=0.04*/, Int_t align/*=12*/)
{
  xjjrootuti::settex(tex, size, align);
  tex->Draw();
}

void xjjrootuti::setleg(TLegend* leg, Float_t size/*=0.04*/)
{
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(size);
}

void xjjrootuti::setlegndraw(TLegend* leg, Float_t size/*=0.04*/)
{
  xjjrootuti::setleg(leg, size);
  leg->Draw();
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

*/


#endif
