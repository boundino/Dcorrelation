#include "../includes/xjjuti.h"
#include "../includes/xjjrootuti.h"
#include "../includes/readD.h"
#include "../includes/fit.h"

const int MAX_XB = 20000;
const int MAX_GEN = 6000;

//
const int nPtBins = 10;
Float_t ptBins[nPtBins+1] = {2, 3, 4, 5, 6, 8, 10, 12.5, 15, 20, 999};
const int nDphiBins_fine = 50;
const int nDphiBins = 10;
Float_t minDphi = 0;
Float_t maxDphi = M_PI;
Float_t dphiBins[nDphiBins+1];
const int nCoBins = 2;
std::map<TString, int> collsyst_list = {{"pp", 0}, {"PbPb", 1}};
const int nhist = 4;
TString histname[nhist] = {"hdphi_all_all", "hdphi_signal_signal", "hdphi_all_signal", "hdphi_signal_all"};
TString histleg[nhist]  = {"all D_{lead}, all D", "g-mat D_{lead}, g-mat D", "all D_{lead}, g-mat D", "g-mat D_{lead}, all D"};
Color_t hcolor[nhist]   = {kBlack, kGreen+3, kAzure-6, kMagenta+3};

std::vector<TString> cutval_list_skim_pp = {"pBeamScrapingFilter", "pPAprimaryVertexFilter"};
std::vector<TString> cutval_list_skim_PbPb = {"pclusterCompatibilityFilter", "pprimaryVertexFilter", "phfCoincFilter3"};

Float_t cutval_list_trkPt[nCoBins][nPtBins]    = {{1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0},   {1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0}};
Float_t cutval_list_trkEta[nCoBins][nPtBins]   = {{1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5},   {1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5,   1.5}};
Float_t cutval_list_trkPtErr[nCoBins][nPtBins] = {{0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3},   {0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3,   0.3}};
Float_t cutval_list_Dy[nCoBins][nPtBins]       = {{1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0},   {1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0}};
Float_t cutval_list_Dalpha[nCoBins][nPtBins]   = {{0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12},  {0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12}};
Float_t cutval_list_Dsvpv[nCoBins][nPtBins]    = {{4.62,  4.80,  4.63,  4.53,  4.09,  4.02,  3.66,  3.70,  3.53,  3.00},  {4.62,  4.80,  4.63,  4.53,  4.09,  4.02,  3.66,  3.70,  3.53,  3.00}};
Float_t cutval_list_Dchi2cl[nCoBins][nPtBins]  = {{0.161, 0.197, 0.141, 0.172, 0.120, 0.098, 0.099, 0.084, 0.047, 0.050}, {0.161, 0.197, 0.141, 0.172, 0.120, 0.098, 0.099, 0.084, 0.047, 0.050}};

TString cutval_list_hlt[nCoBins] = {"HLT_DmesonPPTrackingGlobal_Dpt15_v1", "HLT_HIDmesonHITrackingGlobal_Dpt20_v1"};
std::vector<TString> cutval_list_skim[nCoBins] = {cutval_list_skim_pp, cutval_list_skim_PbPb};

/* ----- */
Float_t cutval_trkPt;
Float_t cutval_trkEta;
Float_t cutval_trkPtErr;
Float_t cutval_Dy;
Float_t cutval_Dsvpv;
Float_t cutval_Dalpha;
Float_t cutval_Dchi2cl;
TString cutval_hlt;
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
