#ifndef _STYBKG_H_
#define _STYBKG_H_

#include "../includes/prefilters.h"
#include "../includes/xjjcuti.h"
#include "../includes/xjjrootuti.h"
#include "../includes/readD.h"
#include "../includes/thgrstyle.h"
#include <TFile.h>
#include <TKey.h>
#include <TH2F.h>
#include <TColor.h>
#include <TCanvas.h>
#include <vector>
#include <map>

std::vector<TString> histnameD = {"DD_Reco", "DDbar_Reco"};
const int nhistD = histnameD.size();
std::vector<TString> histnameG = {"DD_Gen", "DDbar_Gen"};
const int nhistG = histnameG.size();

//
std::map<TString, TH1D*> histptr;
std::map<TString, xjjroot::thgrstyle> histstyle = 
  {
    {"hdphi_DD_Reco", xjjroot::thgrstyle(kAzure+3, 24, 1.1, kAzure+3, 1, 2, -1, -1, -1, std::vector<TString>{"pe", "Reco DD (#bar{D}#bar{D})", "p"})}, 
    {"hdphi_DDbar_Reco", xjjroot::thgrstyle(kAzure+3, 20, 1.1, kAzure+3, 1, 2, -1, -1, -1, std::vector<TString>{"pe", "Reco D#bar{D} (#bar{D}D)", "p"})}, 
    {"hdphi_DD_Gen", xjjroot::thgrstyle(kGreen+4, 25, 1.1, kGreen+4, 1, 2, -1, -1, -1, std::vector<TString>{"pe", "Gen DD (#bar{D}#bar{D})", "p"})}, 
    {"hdphi_DDbar_Gen", xjjroot::thgrstyle(kGreen+4, 21, 1.1, kGreen+4, 1, 2, -1, -1, -1, std::vector<TString>{"pe", "Gen D#bar{D} (#bar{D}D)", "p"})}, 
  };

//
std::vector<TH1D*> ahdphiD(nhistD, 0); // [nhistD];
std::vector<TH1D*> ahdphiG(nhistG, 0); // [nhistG];
std::vector<TH1D*> ahdphi_plot;

//
int createhists(Option_t* option)
{
  TString opt = option;
  opt.ToLower();

  if(opt.Contains("savehist"))
    {
      for(int l=0;l<nhistD;l++)
        ahdphiD[l] = new TH1D(Form("hdphi_%s",histnameD[l].Data()), ";#Delta#phi (rad);Entries (rad^{-1})", nDphiBins, dphiBins.data());
      for(int l=0;l<nhistG;l++)
        ahdphiG[l] = new TH1D(Form("hdphi_%s",histnameG[l].Data()), ";#Delta#phi (rad);Entries (rad^{-1})", nDphiBins, dphiBins.data());
      return 0;
    } 

  std::cout<<"error: invalid option for createhists()"<<std::endl;
  return 1;
}

int writehists(Option_t* option)
{
  TString opt = option;
  opt.ToLower();

  if(opt.Contains("savehist"))
    {
      for(int l=0;l<nhistD;l++)
        ahdphiD[l]->Write();
      for(int l=0;l<nhistG;l++)
        ahdphiG[l]->Write();
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
      if(inf.size()!=1) return 3;
      for(int l=0;l<nhistD;l++)
        ahdphiD[l] = (TH1D*)inf[0]->Get(Form("hdphi_%s",histnameD[l].Data()));
      for(int l=0;l<nhistG;l++)
        ahdphiG[l] = (TH1D*)inf[0]->Get(Form("hdphi_%s",histnameG[l].Data()));
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
