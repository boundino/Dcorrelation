#ifndef _DPHICOR_H_
#define _DPHICOR_H_

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

std::vector<TString> histname = {"all_all", "all_signal", "signal_all", "signal_signal", "sideband_all", "sideband_signal"};
const int nhist = histname.size();
std::vector<Bool_t> histsave = {true,      false,        false,         false,           true,           false};

//
std::map<TString, TH1D*> histptr;
std::map<TString, xjjroot::thgrstyle> histstyle = 
  {
    {"hdphi_all_all",             xjjroot::thgrstyle(kBlack,                       20,  1.1,  kBlack,                       1,  1,  -1,        -1,  -1,    std::vector<TString>{"p",    "all D_{lead}, all D",      "p"})},   //  black
    {"hdphi_all_signal",          xjjroot::thgrstyle(kAzure-6,                     20,  1.1,  kAzure-6,                     1,  1,  -1,        -1,  -1,    std::vector<TString>{"p",    "all D_{lead}, g-mat D",    "p"})},   //  blue
    {"hdphi_signal_signal",       xjjroot::thgrstyle(kGreen+3,                     20,  1.1,  kGreen+3,                     1,  1,  -1,        -1,  -1,    std::vector<TString>{"p",    "g-mat D_{lead}, g-mat D",  "p"})},   //  green
    {"hdphi_all_all_hist",        xjjroot::thgrstyle(-1,                           -1,  -1,   kBlack,                       1,  2,  kGray+1,   -1,  1001,  std::vector<TString>{"hist", "all D_{lead}, all D",      "f"})},   //  black
    {"hdphi_all_signal_hist",     xjjroot::thgrstyle(-1,                           -1,  -1,   kAzure-6,                     1,  2,  kAzure-5,  -1,  1001,  std::vector<TString>{"hist", "all D_{lead}, g-mat D",    "f"})},   //  blue
    {"hdphi_signal_signal_hist",  xjjroot::thgrstyle(-1,                           -1,  -1,   kGreen+3,                     1,  2,  kGreen-5,  -1,  1001,  std::vector<TString>{"hist", "g-mat D_{lead}, g-mat D",  "f"})},   //  green
    {"hdphi_all_all_fit",         xjjroot::thgrstyle(kOrange,                      20,  1.1,  kOrange,                      1,  1,  -1,        -1,  -1,    std::vector<TString>{"pe",   "all D_{lead}, fit D",      "p"})},   //  yellow
    {"hdphi_subtract_signal",     xjjroot::thgrstyle(TColor::GetColor("#ff8faf"),  20,  1.1,  TColor::GetColor("#ff8faf"),  1,  1,  -1,        -1,  -1,    std::vector<TString>{"pe",   "bkgsub D_{lead}, g-mat D", "p"})},   //  pink
    {"hdphi_subtract_all_fit",    xjjroot::thgrstyle(TColor::GetColor("#ed5e5e"),  20,  1.1,  TColor::GetColor("#ed5e5e"),  1,  1,  -1,        -1,  -1,    std::vector<TString>{"pe",   "bkgsub D_{lead}, fit D",   "p"})}    //  red
  };

//
TH1D* hmassSignalLD;                                                             // 
TH1D* hmassSwappedLD;                                                            // 
std::vector<TH1D*> ahmassSignal(nDphiBins, 0);                                   // [nDphiBins];
std::vector<TH1D*> ahmassSwapped(nDphiBins, 0);                                  // [nDphiBins];

TH1D* hmassLD;                                                                   // 
std::vector<TH1D*> ahdphi(nhist, 0);                                             // [nhist];
std::vector<TH1D*> ahmassLD(nhist, 0);                                           // [nhist];
std::vector<std::vector<TH1D*>> ahmass(nhist, std::vector<TH1D*>(nDphiBins, 0)); // [nhist][nDphiBins];

std::vector<TH1D*> ahdphi_fit(nhist, 0);                                         // [nhist];
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
      for(int i=0;i<nDphiBins;i++)
        {
          ahmassSignal[i] = new TH1D(Form("hmassSignal_%d",i), ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
          ahmassSwapped[i] = new TH1D(Form("hmassSwapped_%d",i), ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
        }
      return 0;
    }
  if(opt.Contains("savehist"))
    {
      hmassLD = new TH1D("hmassLD", ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
      for(int l=0;l<nhist;l++)
        {
          ahdphi[l] = new TH1D(Form("hdphi_%s",histname[l].Data()), ";#Delta#phi (rad);Entries (rad^{-1})", nDphiBins, dphiBins.data());
          ahmassLD[l] = new TH1D(Form("hmassLD_%s",histname[l].Data()), ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
          for(int i=0;i<nDphiBins;i++) ahmass[l].at(i) = new TH1D(Form("hmass_%s_%d",histname[l].Data(),i), ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
        }
      return 0;
    }
  if(opt.Contains("usehist"))
    {
      for(int l=0;l<nhist;l++) 
        ahdphi_fit[l] = new TH1D(Form("hdphi_%s_fit",histname[l].Data()), ";#Delta#phi (rad);Entries (rad^{-1})", nDphiBins, dphiBins.data());
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
      for(int i=0;i<nDphiBins;i++)
        {
          ahmassSignal[i]->Write();
          ahmassSwapped[i]->Write();
        }
      return 0;
    }
  if(opt.Contains("savehist"))
    {
      hmassLD->Write();
      for(int l=0;l<nhist;l++)
        {
          ahdphi[l]->Write();
          ahmassLD[l]->Write();
          if(!histsave[l]) continue;
          for(int i=0;i<nDphiBins;i++) ahmass.at(l).at(i)->Write();
        }
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
        {
          ahdphi[l] = (TH1D*)inf[1]->Get(Form("hdphi_%s",histname[l].Data()));
          if(!histsave[l]) continue;
          for(int i=0;i<nDphiBins;i++)
            ahmass[l].at(i) = (TH1D*)inf[1]->Get(Form("hmass_%s_%d",histname[l].Data(),i));
        }
      for(int i=0;i<nDphiBins;i++)
        {
          ahmassSignal[i] = (TH1D*)inf[0]->Get(Form("hmassSignal_%d",i));
          ahmassSwapped[i] = (TH1D*)inf[0]->Get(Form("hmassSwapped_%d",i));
        }
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
