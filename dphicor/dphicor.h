#ifndef _DPHICOR_H_
#define _DPHICOR_H_

#include "../includes/prefilters.h"
#include "../includes/xjjcuti.h"
#include "../includes/xjjrootuti.h"
#include "../includes/readD.h"
#include "../includes/dfitter.h"
#include "../includes/thgrstyle.h"
#include <TFile.h>
#include <TH2F.h>
#include <TColor.h>

TString histname[] = {"all_all", "all_signal", "signal_all", "signal_signal", "sideband_all", "sideband_signal"};
const int nhist = sizeof(histname)/sizeof(histname[0]);
Bool_t  histsave[nhist] = {true,      false,        false,         false,           true,           false};

//
std::map<TString, xjjroot::thgrstyle> histstyle = 
  {
    {"hdphi_all_all",             xjjroot::thgrstyle(kBlack,                       20,  1.1,  kBlack,                       1,  1,  -1,        -1,  -1,    "p")},     //  black
    {"hdphi_all_signal",          xjjroot::thgrstyle(kAzure-6,                     20,  1.1,  kAzure-6,                     1,  1,  -1,        -1,  -1,    "p")},     //  blue
    {"hdphi_signal_signal",       xjjroot::thgrstyle(kGreen+3,                     20,  1.1,  kGreen+3,                     1,  1,  -1,        -1,  -1,    "p")},     //  green
    {"hdphi_all_all_hist",        xjjroot::thgrstyle(-1,                           -1,  -1,   kBlack,                       1,  2,  kGray+1,   -1,  1001,  "hist")},  //  black
    {"hdphi_all_signal_hist",     xjjroot::thgrstyle(-1,                           -1,  -1,   kAzure-6,                     1,  2,  kAzure-5,  -1,  1001,  "hist")},  //  blue
    {"hdphi_signal_signal_hist",  xjjroot::thgrstyle(-1,                           -1,  -1,   kGreen+3,                     1,  2,  kGreen-5,  -1,  1001,  "hist")},  //  green
    {"hdphi_all_all_fit",         xjjroot::thgrstyle(kOrange,                      20,  1.1,  kOrange,                      1,  1,  -1,        -1,  -1,    "pe")},    //  yellow
    {"hdphi_subtract_signal",     xjjroot::thgrstyle(TColor::GetColor("#ff8faf"),  20,  1.1,  TColor::GetColor("#ff8faf"),  1,  1,  -1,        -1,  -1,    "pe")},    //  pink
    {"hdphi_subtract_all_fit",    xjjroot::thgrstyle(TColor::GetColor("#ed5e5e"),  20,  1.1,  TColor::GetColor("#ed5e5e"),  1,  1,  -1,        -1,  -1,    "pe")}     //  red
  };

std::map<TString, TString> histleg = 
  {
    {"hdphi_all_all",             "all D_{lead}, all D"}, 
    {"hdphi_all_signal",          "all D_{lead}, g-mat D"}, 
    {"hdphi_signal_signal",       "g-mat D_{lead}, g-mat D"}, 
    {"hdphi_all_all_hist",        "all D_{lead}, all D"}, 
    {"hdphi_all_signal_hist",     "all D_{lead}, g-mat D"}, 
    {"hdphi_signal_signal_hist",  "g-mat D_{lead}, g-mat D"}, 
    {"hdphi_all_all_fit",         "all D_{lead}, fit D"}, 
    {"hdphi_subtract_signal",     "bkgsub D_{lead}, g-mat D"}, 
    {"hdphi_subtract_all_fit",    "bkgsub D_{lead}, fit D"}
  };

TH1D* rebindiffhist(const TH1D* h, Int_t nnewbin, Double_t* newbin, TString newname)
{
  TH1D* h_mulbin = (TH1D*)h->Clone("h_mulbin");
  xjjroot::multiplebinwid(h_mulbin);
  TH1D* h_rebin = (TH1D*)h_mulbin->Rebin(nnewbin, newname, newbin);
  xjjroot::dividebinwid(h_rebin);
  delete h_mulbin;
  return h_rebin;
}

#endif
