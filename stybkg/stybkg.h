#ifndef _STYBKG_H_
#define _STYBKG_H_

#include "../includes/prefilters.h"
#include "../includes/xjjcuti.h"
#include "../includes/xjjrootuti.h"
#include "../includes/readD.h"
#include "../includes/dfitter.h"
#include "../includes/thgrstyle.h"
#include <TFile.h>
#include <TKey.h>
#include <TH2F.h>
#include <TColor.h>
#include <vector>
#include <map>

std::vector<TString> histname = {"incl_all_signal", "incl_signal_signal", "incl_swap_signal", "incl_comb_signal", "sideband_all_signal", "sideband_signal_signal", "sideband_swap_signal", "sideband_comb_signal"};
const int nhist = histname.size();

//
std::map<TString, TH1D*> histptr;
std::map<TString, xjjroot::thgrstyle> histstyle = 
  {
    {"hdphi_incl_signal_signal_hist",            xjjroot::thgrstyle(-1,         -1,  -1,   kOrange+2,  1,  2,  kOrange-2,  -1,  1001,  std::vector<TString>{"hist", "signal D_{lead}",          "f"})},  
    {"hdphi_incl_signalNswap_signal_hist",       xjjroot::thgrstyle(-1,         -1,  -1,   kGreen+3,   1,  2,  kGreen-2,   -1,  1001,  std::vector<TString>{"hist", "swap D_{lead}",            "f"})},  
    {"hdphi_incl_signalNswapNcomb_signal_hist",  xjjroot::thgrstyle(-1,         -1,  -1,   kAzure+3,   1,  2,  kAzure+5,   -1,  1001,  std::vector<TString>{"hist", "comb D_{lead}",            "f"})},  
    {"hdphi_incl_all_signal",                    xjjroot::thgrstyle(kBlack,     20,  1.1,  kBlack,     1,  1,  -1,         -1,  -1,    std::vector<TString>{"pe",   "all D_{lead}",             "p"})},
    {"hdphi_incl_signal_signal_norm",            xjjroot::thgrstyle(kOrange-3,  20,  1.1,  kOrange-3,  1,  1,  -1,         -1,  -1,    std::vector<TString>{"pe",   "signal D_{lead}",          "p"})},
    {"hdphi_incl_swap_signal_norm",              xjjroot::thgrstyle(kGreen+4,   20,  1.1,  kGreen+4,   1,  1,  -1,         -1,  -1,    std::vector<TString>{"pe",   "swap D_{lead}",            "p"})},
    {"hdphi_incl_comb_signal_norm",              xjjroot::thgrstyle(kBlue,      20,  1.1,  kBlue,      1,  1,  -1,         -1,  -1,    std::vector<TString>{"pe",   "comb D_{lead}",            "p"})},
    {"hdphi_sideband_signalNswap_signal_norm",   xjjroot::thgrstyle(kViolet+8,  20,  1.1,  kViolet+8,  1,  1,  -1,         -1,  -1,    std::vector<TString>{"pe",   "sideband signal D_{lead}", "p"})},
    {"hdphi_sideband_comb_signal_norm",          xjjroot::thgrstyle(kAzure+3,   20,  1.1,  kAzure+3,   1,  1,  -1,         -1,  -1,    std::vector<TString>{"pe",   "sideband bkg D_{lead}",    "p"})},
    {"hdphi_incl_signalNswap_signal_norm_hist",  xjjroot::thgrstyle(-1,         -1,  -1,   kViolet+8,  1,  2,  kViolet-9,  -1,  1001,  std::vector<TString>{"hist", "all signal D_{lead}",      "f"})},
    {"hdphi_incl_comb_signal_norm_hist",         xjjroot::thgrstyle(-1,         20,  1.1,  kAzure+3,   1,  2,  kAzure+5,   -1,  1001,  std::vector<TString>{"hist", "all bkg D_{lead}",         "f"})},
    {"hdphi_incl_all_signal_subtract",           xjjroot::thgrstyle(kRed+1,     20,  1.1,  kRed+1,     1,  1,  -1,         -1,  -1,    std::vector<TString>{"pe",   "bkgsub D_{lead}",          "p"})}
  };

//
TH1D* hmassSignalLD;                                                             //
TH1D* hmassSwappedLD;                                                            //
TH1D* hmassLD;                                                                   //
std::vector<TH1D*> ahdphi(nhist, 0);                                             // [nhist];
std::vector<TH1D*> ahdphi_plot;

//
int createhists(Option_t* option)
{
  TString opt = option;
  opt.ToLower();

  if(opt.Contains("savetpl"))
    {
      hmassSignalLD = new TH1D("hmassSignalLD", ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
      hmassSwappedLD = new TH1D("hmassSwappedLD", ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
      return 0;
    }
  if(opt.Contains("savehist"))
    {
      hmassLD = new TH1D("hmassLD", ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
      for(int l=0;l<nhist;l++)
        ahdphi[l] = new TH1D(Form("hdphi_%s",histname[l].Data()), ";#Delta#phi (rad);Entries (rad^{-1})", nDphiBins, dphiBins.data());
      return 0;
    } 

  std::cout<<"error: invalid option for createhists()"<<std::endl;
  return 1;
}

int writehists(Option_t* option)
{
  TString opt = option;
  opt.ToLower();

  if(opt.Contains("savetpl"))
    {
      hmassSignalLD->Write();
      hmassSwappedLD->Write();
      return 0;
    }
  if(opt.Contains("savehist"))
    {
      hmassLD->Write();
      for(int l=0;l<nhist;l++)
        ahdphi[l]->Write();
      return 0;
    }
  if(opt.Contains("usehist"))
    {
      for(std::vector<TH1D*>::iterator it=ahdphi_plot.begin(); it!=ahdphi_plot.end(); it++)
        (*it)->Write();
      return 0;
    }

  std::cout<<"error: invalid option for writehists()"<<std::endl;
  return 1;
}

int gethists(std::vector<TFile*> inf, Option_t* option)
{
  TString opt = option;
  opt.ToLower();

  if(opt.Contains("usehist"))
    {
      if(inf.size()!=2) return 3;
      hmassLD = (TH1D*)inf[1]->Get("hmassLD");
      hmassSignalLD = (TH1D*)inf[0]->Get("hmassSignalLD");
      hmassSwappedLD = (TH1D*)inf[0]->Get("hmassSwappedLD");
      for(int l=0;l<nhist;l++)
        ahdphi[l] = (TH1D*)inf[1]->Get(Form("hdphi_%s",histname[l].Data()));
      return 0;
    }

  if(opt.Contains("plothist"))
    {
      if(inf.size()!=1) return 3;
      TIter next(inf[0]->GetListOfKeys());
      TKey* key;
      while((key = (TKey*)next()))
        {
          if(TString(key->GetClassName(), 4)!="TH1D") continue;
          histptr.insert(std::pair<TString, TH1D*>(key->GetName(), (TH1D*)inf[0]->Get(key->GetName())));
        }
      return 0;
    }

  std::cout<<"error: invalid option for gethists()"<<std::endl;
  return 1;
}


#endif
