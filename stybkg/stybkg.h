#ifndef _STYBKG_H_
#define _STYBKG_H_

#include "../includes/prefilters.h"
#include "../includes/xjjcuti.h"
#include "../includes/xjjrootuti.h"
#include "../includes/readD.h"
#include "../includes/dfitter.h"
#include "../includes/thgrstyle.h"
#include <TFile.h>
#include <TH2F.h>
#include <TColor.h>


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


#endif
