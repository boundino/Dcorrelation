#ifndef _PREFILTERS_H_
#define _PREFILTERS_H_

#include <vector>
#include <map>
#include <iostream>
#include <TString.h>

//
Double_t dmass_sideband_l = 0.07;
Double_t dmass_sideband_h = 0.12;

std::vector<Double_t> ptBins = {0, 3, 4, 5, 6, 8, 10, 12.5, 15, 20, 999};
const int nPtBins = ptBins.size()-1;

Double_t minDphi = 0;
Double_t maxDphi = M_PI;
const int nDphiBins = 11;
std::vector<Double_t> fphiBins = {0., 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
std::vector<Double_t> dphiBins; 

const int nCoBins = 2;
std::map<TString, int> collsyst_list = {{"pp", 0}, {"PbPb", 1}};

std::vector<Double_t> ptasymBins = {0, 0.2, 0.4, 0.6, 0.8, 1};
const int nPtasymBins = ptasymBins.size()-1;

TString tMC[] = {"data", "MC"};

//
const int MAX_XB = 20000;
const int MAX_GEN = 6000;
const double MASS_DZERO = 1.8649;

const Double_t n_hist_dzero = 60;
const Double_t min_hist_dzero = 1.7;
const Double_t max_hist_dzero = 2.0;
const Double_t binwid_hist_dzero = (max_hist_dzero-min_hist_dzero)/n_hist_dzero;

//
std::vector<std::vector<Double_t>> cutval_list_trkPt    = {std::vector<Double_t>{0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5},   std::vector<Double_t>{0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5}};
std::vector<std::vector<Double_t>> cutval_list_trkEta   = {std::vector<Double_t>{1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5},   std::vector<Double_t>{1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5}};
std::vector<std::vector<Double_t>> cutval_list_trkPtErr = {std::vector<Double_t>{0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3},   std::vector<Double_t>{0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3}};
std::vector<std::vector<Double_t>> cutval_list_Dy       = {std::vector<Double_t>{1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0},   std::vector<Double_t>{1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0}};
std::vector<std::vector<Double_t>> cutval_list_Dalpha   = {std::vector<Double_t>{0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12},  std::vector<Double_t>{0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12}};
std::vector<std::vector<Double_t>> cutval_list_Dsvpv    = {std::vector<Double_t>{4.62,  4.80,  4.63,  4.53,  4.09,  4.02,  3.66,  3.70,  3.53,  3.00},  std::vector<Double_t>{4.62,  4.80,  4.63,  4.53,  4.09,  4.02,  3.66,  3.70,  3.53,  3.00}};
std::vector<std::vector<Double_t>> cutval_list_Dchi2cl  = {std::vector<Double_t>{0.161, 0.197, 0.141, 0.172, 0.120, 0.098, 0.099, 0.084, 0.047, 0.050}, std::vector<Double_t>{0.161, 0.197, 0.141, 0.172, 0.120, 0.098, 0.099, 0.084, 0.047, 0.050}};

std::vector<std::vector<TString>> cutval_list_hlt       = {std::vector<TString>{"HLT_DmesonPPTrackingGlobal_Dpt15_v1"},                                 std::vector<TString>{"HLT_HIDmesonHITrackingGlobal_Dpt20_v1"}};
std::vector<std::vector<TString>> cutval_list_skim      = {std::vector<TString>{"pBeamScrapingFilter", "pPAprimaryVertexFilter"},                       std::vector<TString>{"pclusterCompatibilityFilter", "pprimaryVertexFilter", "phfCoincFilter3"}};

/* ----- */
Double_t cutval_trkPt;
Double_t cutval_trkEta;
Double_t cutval_trkPtErr;
Double_t cutval_Dy;
Double_t cutval_Dsvpv;
Double_t cutval_Dalpha;
Double_t cutval_Dchi2cl;
std::vector<TString> cutval_hlt;
std::vector<TString> cutval_skim;
Int_t npass;

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
  return 0;
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
  cutval_trkPt = cutval_list_trkPt[icollsyst].at(ipt);
  cutval_trkEta = cutval_list_trkEta[icollsyst].at(ipt);
  cutval_trkPtErr = cutval_list_trkPtErr[icollsyst].at(ipt);
  cutval_Dy = cutval_list_Dy[icollsyst].at(ipt);
  cutval_Dsvpv = cutval_list_Dsvpv[icollsyst].at(ipt);
  cutval_Dalpha = cutval_list_Dalpha[icollsyst].at(ipt);
  cutval_Dchi2cl = cutval_list_Dchi2cl[icollsyst].at(ipt);
  return 0;
}

int initbinning()
{
  if(nDphiBins!=fphiBins.size()-1)
    {
      std::cout<<"\033[1;31merror:\033[0m \"fphiBins.size()\" is not consistent with \"nDphiBins\"."<<std::endl;
      return 1;
    }

  Double_t frac = 0;
  for(int i=0;i<=nDphiBins;i++) dphiBins.push_back((frac=frac+fphiBins[i])*(maxDphi-minDphi));
  if(nDphiBins!=dphiBins.size()-1)
    {
      std::cout<<"\033[1;31merror:\033[0m \"dphiBins.size()\" is not consistent with \"nDphiBins\"."<<std::endl;
      return 1;
    }
  return 0;
}

#endif
