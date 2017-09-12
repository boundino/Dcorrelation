// #include "TRandom3.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"

#include "D_track_tree.h"
#include "D_tree.h"
#include "G_tree.h"
#include "track_tree.h"
#include "genpart_tree.h"
#include "jet_tree.h"

#include "Corrections/getTrkCorr.h"

#include <stdlib.h>
#include <stdint.h>
#include <functional>
#include <algorithm>

// static const double pi = 3.141592653589793238462643383279502884;
double getAngleToEP(double angle);
float getTrkWeight(TrkCorr* trkCorr, int itrk, int hiBin, jetTree* jt_trkcorr, trackTree* tt);

int D_track_skim(std::string input, std::string output, 
                 bool isPP, bool isMC, 
                 // float trkptmin = 10, 
                 int start = 0, int end = -1) 
{
  bool isHI = !isPP;

  /**********************************************************
   * CREATE OUTPUT TREE
   **********************************************************/
  TFile* foutput = new TFile(output.c_str(), "recreate");

  TTree* outtree = new TTree("dtrk", "");
  DTrackTree dtrk(outtree);
  /**********************************************************
   * OPEN INPUT FILES
   **********************************************************/
  TFile* finput = TFile::Open(input.c_str(), "read");

#define _SET_BRANCH_ADDRESS(tree, branch, var) {        \
    tree->SetBranchStatus(#branch, 1);                  \
    tree->SetBranchAddress(#branch, &var);              \
  }

  TTree* event_tree = (TTree*)finput->Get("hiEvtAnalyzer/HiTree");
  if (!event_tree) { printf("Could not access event tree!\n"); return 1; }
  event_tree->SetBranchStatus("*", 0);
  int hiBin;
  float vz;
  float hiEvtPlanes[29];
  _SET_BRANCH_ADDRESS(event_tree, run, dtrk.run);
  _SET_BRANCH_ADDRESS(event_tree, evt, dtrk.evt);
  _SET_BRANCH_ADDRESS(event_tree, lumi, dtrk.lumi);
  _SET_BRANCH_ADDRESS(event_tree, hiBin, hiBin);
  _SET_BRANCH_ADDRESS(event_tree, vz, vz);
  _SET_BRANCH_ADDRESS(event_tree, weight, dtrk.weight);
  _SET_BRANCH_ADDRESS(event_tree, pthat, dtrk.pthat);
  _SET_BRANCH_ADDRESS(event_tree, hiEvtPlanes, hiEvtPlanes);

  TTree* skim_tree = (TTree*)finput->Get("skimanalysis/HltTree");
  if(!skim_tree) { printf("Could not access skim tree!\n"); return 1; }
  skim_tree->SetBranchStatus("*", 0);
  // int pcollisionEventSelection;
  // int HBHENoiseFilterResultRun2Loose;
  int pPAprimaryVertexFilter;
  int pBeamScrapingFilter;
  int pclusterCompatibilityFilter;
  int pprimaryVertexFilter;
  int phfCoincFilter3;
  // _SET_BRANCH_ADDRESS(skim_tree, pcollisionEventSelection, pcollisionEventSelection);
  // _SET_BRANCH_ADDRESS(skim_tree, HBHENoiseFilterResultRun2Loose, HBHENoiseFilterResultRun2Loose);
  if(isPP)
    {
      _SET_BRANCH_ADDRESS(skim_tree, pPAprimaryVertexFilter, pPAprimaryVertexFilter);
      _SET_BRANCH_ADDRESS(skim_tree, pBeamScrapingFilter, pBeamScrapingFilter);
    }
  else
    {
      _SET_BRANCH_ADDRESS(skim_tree, pclusterCompatibilityFilter, pclusterCompatibilityFilter);
      _SET_BRANCH_ADDRESS(skim_tree, pprimaryVertexFilter, pprimaryVertexFilter);
      _SET_BRANCH_ADDRESS(skim_tree, phfCoincFilter3, phfCoincFilter3);
    }
  TTree* hlt_tree = (TTree*)finput->Get("hltanalysis/HltTree");
  if(!hlt_tree) { printf("Could not access hlt tree!\n"); return 1; }
  // hlt_tree->SetBranchStatus("*", 0);

  foutput->cd();  
  TTree* hlt_new = hlt_tree->CloneTree(0);
  hlt_new->SetName("hlt");

  TTree* D_tree = (TTree*)finput->Get("Dfinder/ntDkpi");
  if(!D_tree) {printf("Could not access D tree!\n"); return 1; }
  D_tree->SetBranchStatus("*", 0);
  DTree dt(D_tree);

  TTree* G_tree = (TTree*)finput->Get("Dfinder/ntGen");
  if(!G_tree) {printf("Could not access G tree!\n"); return 1; }
  G_tree->SetBranchStatus("*", 0);
  GTree gt(G_tree);

  TTree* jet_tree_for_trk_corr = isPP ? (TTree*)finput->Get("ak4CaloJetAnalyzer/t") : (TTree*)finput->Get("akPu4CaloJetAnalyzer/t");
  if (!jet_tree_for_trk_corr) { printf("Could not access jet tree for track corrections!\n"); return 1; }
  jet_tree_for_trk_corr->SetBranchStatus("*", 0);
  jetTree jt_trkcorr(jet_tree_for_trk_corr);

  TTree* track_tree = isPP ? (TTree*)finput->Get("ppTrack/trackTree") : (TTree*)finput->Get("anaTrack/trackTree");
  if(!track_tree) { printf("Could not access track tree!\n"); return 1; }
  track_tree->SetBranchStatus("*", 0);
  trackTree tt(track_tree);

  TTree* genpart_tree = (TTree*)finput->Get("HiGenParticleAna/hi");
  if(!genpart_tree) { printf("Could not access gen tree!\n"); isMC = false; }
  genpartTree gpt;
  if(isMC) 
    {
      genpart_tree->SetBranchStatus("*", 0);
      gpt.read_tree(genpart_tree);
    }

  /**********************************************************
   * OPEN CORRECTION FILES
   **********************************************************/
  TrkCorr* trkCorr;
  if(isHI)
    trkCorr = new TrkCorr("Corrections/TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/");
  else
    trkCorr = new TrkCorr("Corrections/TrkCorr_July22_Iterative_pp_eta2p4/");

  /**********************************************************
   * BEGIN EVENT LOOP
   **********************************************************/
  int nevents = event_tree->GetEntries();
  for (int j = start; j < nevents; j++) 
    {
      dtrk.clear_vectors();

      skim_tree->GetEntry(j);
      event_tree->GetEntry(j);

      hlt_tree->GetEntry(j);
      if(j%1000==0) { printf("processing event: %i / %i\n", j, nevents); }
      if(j==end) { printf("done: %i\n", end); break; }

      if(fabs(vz) > 15) continue;
      if(!isPP) 
        { // HI event selection
          if(pclusterCompatibilityFilter < 1 || pprimaryVertexFilter < 1 || phfCoincFilter3 < 1) continue;
        } 
      else 
        { // pp event selection
          if(pPAprimaryVertexFilter < 1 || pBeamScrapingFilter < 1) continue;
        }

      //! D cuts and selection
      D_tree->GetEntry(j);
      for (int id = 0; id < dt.Dsize; ++id) 
        {
          dtrk.copy_index(dt, id);
        }
      dtrk.copy_variables(dt);

      //! Gen cuts and selection
      G_tree->GetEntry(j);
      for (int id = 0; id < gt.Gsize; ++id) 
        {
          dtrk.copy_index_gen(gt, id);
        }
      dtrk.copy_variables_gen(gt);

      int nTrk = 0;

      jet_tree_for_trk_corr->GetEntry(j);
      float maxJetPt = -999;
      for (int k = 0; k < jt_trkcorr.nref; k++) {
        if (TMath::Abs(jt_trkcorr.jteta[k]) > 2) continue;
        if (jt_trkcorr.jtpt[k] > maxJetPt) maxJetPt = jt_trkcorr.jtpt[k];
      }

      float maxTrkPt = -999;
      //! (2.4) Begin track cuts and selection
      track_tree->GetEntry(j);
      for(int itrk = 0; itrk < tt.nTrk; ++itrk) 
        {
          if(tt.trkPt[itrk] < 1 || tt.trkPt[itrk] > 300 || fabs(tt.trkEta[itrk]) > 2.4) continue;
          if(tt.highPurity[itrk] != 1) continue;
          if(tt.trkPtError[itrk] / tt.trkPt[itrk] > 0.1 || TMath::Abs(tt.trkDz1[itrk] / tt.trkDzError1[itrk]) > 3 || TMath::Abs(tt.trkDxy1[itrk] / tt.trkDxyError1[itrk]) > 3) continue;
          if(!isPP && tt.trkChi2[itrk] / (float)tt.trkNdof[itrk] / (float)tt.trkNlayer[itrk] > 0.15) continue;
          if(!isPP && tt.trkNHit[itrk] < 11) continue;

          float Et = (tt.pfHcal[itrk] + tt.pfEcal[itrk]) / TMath::CosH(tt.trkEta[itrk]);
          if(!(tt.trkPt[itrk] < 20 || (Et > 0.5 * tt.trkPt[itrk]))) continue;
          if(tt.trkPt[itrk] > maxTrkPt) maxTrkPt = tt.trkPt[itrk];
          float trkWeight = 0;
          if(isPP) trkWeight = getTrkWeight(trkCorr, itrk, 0, &jt_trkcorr, &tt);
          else trkWeight = getTrkWeight(trkCorr, itrk, hiBin, &jt_trkcorr, &tt);

          dtrk.trkPt.push_back(tt.trkPt[itrk]);
          dtrk.trkPtError.push_back(tt.trkPtError[itrk]);
          dtrk.trkNHit.push_back(tt.trkNHit[itrk]);
          dtrk.trkNlayer.push_back(tt.trkNlayer[itrk]);
          dtrk.trkEta.push_back(tt.trkEta[itrk]);
          dtrk.trkPhi.push_back(tt.trkPhi[itrk]);
          dtrk.trkCharge.push_back(tt.trkCharge[itrk]);
          dtrk.highPurity.push_back(tt.highPurity[itrk]);
          dtrk.trkChi2.push_back(tt.trkChi2[itrk]);
          dtrk.trkNdof.push_back(tt.trkNdof[itrk]);
          dtrk.trkDxy1.push_back(tt.trkDxy1[itrk]);
          dtrk.trkDxyError1.push_back(tt.trkDxyError1[itrk]);
          dtrk.trkDz1.push_back(tt.trkDz1[itrk]);
          dtrk.trkDzError1.push_back(tt.trkDzError1[itrk]);
          dtrk.pfEcal.push_back(tt.pfEcal[itrk]);
          dtrk.pfHcal.push_back(tt.pfHcal[itrk]);
          dtrk.trkWeight.push_back(trkWeight);
          nTrk++;
        }
      dtrk.nTrk = nTrk;
      //! End track selection

      if(isMC)
        {
          genpart_tree->GetEntry(j);
          dtrk.mult = gpt.mult;
          for(int igenp = 0; igenp < gpt.mult; ++igenp) 
            {
              dtrk.pt.push_back((*gpt.pt)[igenp]);
              dtrk.eta.push_back((*gpt.eta)[igenp]);
              dtrk.phi.push_back((*gpt.phi)[igenp]);
              dtrk.pdg.push_back((*gpt.pdg)[igenp]);
              dtrk.chg.push_back((*gpt.chg)[igenp]);
              dtrk.matchingID.push_back((*gpt.matchingID)[igenp]);
              dtrk.sube.push_back((*gpt.sube)[igenp]);
            }
        }

      //

      dtrk.isPP = isPP;
      dtrk.hiBin = hiBin;
      dtrk.vz = vz;
      memcpy(dtrk.hiEvtPlanes, hiEvtPlanes, 29 * sizeof(float));

      outtree->Fill();
      hlt_new->Fill();
    }

  foutput->cd();
  outtree->Write("", TObject::kOverwrite);
  hlt_new->Write("", TObject::kOverwrite);
  foutput->Write("", TObject::kOverwrite);
  foutput->Close();

  return 0;
}

double getAngleToEP(double angle) 
{
  angle = (angle > TMath::Pi()) ? 2 * TMath::Pi() - angle : angle;
  return (angle > TMath::Pi() / 2) ? TMath::Pi() - angle : angle;
}

float getTrkWeight(TrkCorr* trkCorr, int itrk, int hiBin, jetTree* jt_trkcorr, trackTree* tt) 
{
  float rmin = 999;
  for (int k = 0; k < jt_trkcorr->nref; k++) 
    {
      if (jt_trkcorr->jtpt[k] < 50) break;
      if ((TMath::Abs(jt_trkcorr->chargedSum[k] / jt_trkcorr->rawpt[k]) < 0.01) || (TMath::Abs(jt_trkcorr->jteta[k] > 2))) continue;
      float R = TMath::Power(jt_trkcorr->jteta[k] - tt->trkEta[itrk], 2) + TMath::Power(TMath::ACos(TMath::Cos(jt_trkcorr->jtphi[k] - tt->trkPhi[itrk])), 2);
      if (rmin * rmin > R) rmin = TMath::Power(R, 0.5);
    }
  return trkCorr->getTrkCorr(tt->trkPt[itrk], tt->trkEta[itrk], tt->trkPhi[itrk], hiBin, rmin);
}

int main(int argc, char* argv[]) 
{
  if(argc == 5)
    return D_track_skim(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]));
  else if(argc == 6)
    return D_track_skim(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
  else if(argc == 7)
    return D_track_skim(argv[1], argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
  else
    printf("Usage: ./D_track_skim.exe [input] [output] [isPP] [isMC] [start] [end]\n");

  return 1;
}
