#ifndef _STYBKG_H_
#define _STYBKG_H_

#include "../includes/xjjcuti.h"
#include "../includes/xjjrootuti.h"
#include "../includes/readD.h"
#include "../includes/dfitter.h"
#include "../includes/thgrstyle.h"
#include <TFile.h>
#include <TH2F.h>
#include <TColor.h>

const int MAX_XB = 20000;
const int MAX_GEN = 6000;
const double MASS_DZERO = 1.8649;

const Double_t n_hist_dzero = 60;
const Double_t min_hist_dzero = 1.7;
const Double_t max_hist_dzero = 2.0;
const Double_t binwid_hist_dzero = (max_hist_dzero-min_hist_dzero)/n_hist_dzero;

//
Double_t dmass_sideband_l = 0.07;
Double_t dmass_sideband_h = 0.12;

const int nPtBins = 10;
Double_t ptBins[nPtBins+1] = {2, 3, 4, 5, 6, 8, 10, 12.5, 15, 20, 999};
const int nDphiBins_fine = 20; // 50
Double_t minDphi = 0;
Double_t maxDphi = M_PI;
const int nDphiBins = 11;
Double_t fphiBins[nDphiBins+1] = {0., 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
Double_t dphiBins[nDphiBins+1]; 
const int nCoBins = 2;
std::map<TString, int> collsyst_list = {{"pp", 0}, {"PbPb", 1}};
TString histname[] = {"incl_all_signal", "incl_signal_signal", "incl_swap_signal", "incl_comb_signal", "sideband_all_signal", "sideband_signal_signal", "sideband_swap_signal", "sideband_comb_signal"};
const int nhist = sizeof(histname)/sizeof(histname[0]);
// Bool_t  histsave[nhist] = {true,      false,        true,         false,           true,           false};

//
std::map<TString, xjjroot::thgrstyle> histstyle = 
  {
    {"hdphi_incl_signal_signal_hist",            xjjroot::thgrstyle(-1,         -1,  -1,   kOrange+2,  1,  2,  kOrange-2,  -1,  1001,  "hist")},  
    {"hdphi_incl_signalNswap_signal_hist",       xjjroot::thgrstyle(-1,         -1,  -1,   kGreen+3,   1,  2,  kGreen-2,   -1,  1001,  "hist")},  
    {"hdphi_incl_signalNswapNcomb_signal_hist",  xjjroot::thgrstyle(-1,         -1,  -1,   kAzure+3,   1,  2,  kAzure+5,   -1,  1001,  "hist")},  
    {"hdphi_incl_all_signal",                    xjjroot::thgrstyle(kBlack,     20,  1.1,  kBlack,     1,  1,  -1,         -1,  -1,    "pe")},
    {"hdphi_incl_signal_signal_norm",            xjjroot::thgrstyle(kOrange-3,  20,  1.1,  kOrange-3,  1,  1,  -1,         -1,  -1,    "pe")},
    {"hdphi_incl_swap_signal_norm",              xjjroot::thgrstyle(kGreen+4,   20,  1.1,  kGreen+4,   1,  1,  -1,         -1,  -1,    "pe")},
    {"hdphi_incl_comb_signal_norm",              xjjroot::thgrstyle(kBlue,      20,  1.1,  kBlue,      1,  1,  -1,         -1,  -1,    "pe")},
    {"hdphi_sideband_signalNswap_signal_norm",   xjjroot::thgrstyle(kViolet+8,  20,  1.1,  kViolet+8,  1,  1,  -1,         -1,  -1,    "pe")},
    {"hdphi_sideband_comb_signal_norm",          xjjroot::thgrstyle(kAzure+3,   20,  1.1,  kAzure+3,   1,  1,  -1,         -1,  -1,    "pe")},
    {"hdphi_incl_signalNswap_signal_norm_hist",  xjjroot::thgrstyle(-1,         -1,  -1,   kViolet+8,  1,  2,  kViolet-9,  -1,  1001,  "hist")},
    {"hdphi_incl_comb_signal_norm_hist",         xjjroot::thgrstyle(-1,         20,  1.1,  kAzure+3,   1,  2,  kAzure+5,   -1,  1001,  "hist")},
    {"hdphi_incl_all_signal_subtract",           xjjroot::thgrstyle(kRed+1,     20,  1.1,  kRed+1,     1,  1,  -1,         -1,  -1,    "pe")}
  };

std::map<TString, TString> histleg = 
  {
    {"hdphi_incl_signal_signal_hist",            "signal D_{lead}"}, 
    {"hdphi_incl_signalNswap_signal_hist",       "swap D_{lead}"},
    {"hdphi_incl_signalNswapNcomb_signal_hist",  "comb D_{lead}"},
    {"hdphi_incl_all_signal",                    "all D_{lead}"},
    {"hdphi_incl_signal_signal_norm",            "signal D_{lead}"},
    {"hdphi_incl_swap_signal_norm",              "swap D_{lead}"},
    {"hdphi_incl_comb_signal_norm",              "comb D_{lead}"},
    {"hdphi_sideband_signalNswap_signal_norm",   "sideband signal D_{lead}"},
    {"hdphi_sideband_comb_signal_norm",          "sideband bkg D_{lead}"},
    {"hdphi_incl_signalNswap_signal_norm_hist",  "all signal D_{lead}"},
    {"hdphi_incl_comb_signal_norm_hist",         "all bkg D_{lead}"},
    {"hdphi_incl_all_signal_subtract",           "bkgsub D_{lead}"}
  };


//
std::vector<TString> cutval_list_skim_pp = {"pBeamScrapingFilter", "pPAprimaryVertexFilter"};
std::vector<TString> cutval_list_skim_PbPb = {"pclusterCompatibilityFilter", "pprimaryVertexFilter", "phfCoincFilter3"};

Double_t cutval_list_trkPt[nCoBins][nPtBins]    = {{1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0},   {1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0}};
Double_t cutval_list_trkEta[nCoBins][nPtBins]   = {{1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5},   {1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5}};
Double_t cutval_list_trkPtErr[nCoBins][nPtBins] = {{0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3},   {0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3}};
Double_t cutval_list_Dy[nCoBins][nPtBins]       = {{1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0},   {1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0}};
Double_t cutval_list_Dalpha[nCoBins][nPtBins]   = {{0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12},  {0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12}};
Double_t cutval_list_Dsvpv[nCoBins][nPtBins]    = {{4.62,  4.80,  4.63,  4.53,  4.09,  4.02,  3.66,  3.70,  3.53,  3.00},  {4.62,  4.80,  4.63,  4.53,  4.09,  4.02,  3.66,  3.70,  3.53,  3.00}};
Double_t cutval_list_Dchi2cl[nCoBins][nPtBins]  = {{0.161, 0.197, 0.141, 0.172, 0.120, 0.098, 0.099, 0.084, 0.047, 0.050}, {0.161, 0.197, 0.141, 0.172, 0.120, 0.098, 0.099, 0.084, 0.047, 0.050}};

TString cutval_list_hlt[nCoBins] = {"HLT_DmesonPPTrackingGlobal_Dpt15_v1", "HLT_HIDmesonHITrackingGlobal_Dpt20_v1"};
std::vector<TString> cutval_list_skim[nCoBins] = {cutval_list_skim_pp, cutval_list_skim_PbPb};

/* ----- */
Double_t cutval_trkPt;
Double_t cutval_trkEta;
Double_t cutval_trkPtErr;
Double_t cutval_Dy;
Double_t cutval_Dsvpv;
Double_t cutval_Dalpha;
Double_t cutval_Dchi2cl;
TString  cutval_hlt;
std::vector<TString> cutval_skim;

//

int initcutval(TString collisionsyst)
{
  if(collsyst_list.find(collisionsyst)==collsyst_list.end())
    {
      std::cout<<"\033[1;31merror:\033[0m invalid \"collisionsyst\" - initcutval()"<<std::endl;
      return 1;
    }
  int icollsyst = collsyst_list[collisionsyst];
  cutval_hlt = cutval_list_hlt[icollsyst];
  cutval_skim = cutval_list_skim[icollsyst];
}

int initcutval_ptdep(TString collisionsyst, int ipt)
{
  if(collsyst_list.find(collisionsyst)==collsyst_list.end())
    {
      std::cout<<"\033[1;31merror:\033[0m invalid \"collisionsyst\" - initcutval()"<<std::endl;
      return 1;
    }
  if(ipt<0 || ipt>=nPtBins) return 2;

  int icollsyst = collsyst_list[collisionsyst];
  cutval_trkPt = cutval_list_trkPt[icollsyst][ipt];
  cutval_trkEta = cutval_list_trkEta[icollsyst][ipt];
  cutval_trkPtErr = cutval_list_trkPtErr[icollsyst][ipt];
  cutval_Dy = cutval_list_Dy[icollsyst][ipt];
  cutval_Dsvpv = cutval_list_Dsvpv[icollsyst][ipt];
  cutval_Dalpha = cutval_list_Dalpha[icollsyst][ipt];
  cutval_Dchi2cl = cutval_list_Dchi2cl[icollsyst][ipt];
  return 0;
}

void initbinning()
{
  // for(int i=0;i<=nDphiBins;i++) dphiBins[i] = minDphi+i*(maxDphi-minDphi)/nDphiBins;
  dphiBins[0] = minDphi;
  for(int i=1;i<=nDphiBins;i++) dphiBins[i] = dphiBins[i-1]+fphiBins[i]*(maxDphi-minDphi);
}

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
