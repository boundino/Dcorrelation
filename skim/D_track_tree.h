#ifndef _D_TRACK_TREE_H
#define _D_TRACK_TREE_H

#include "TTree.h"

#include "D_tree.h"
#include "G_tree.h"

#include <vector>

class DTrackTree 
{
public:
  DTrackTree() 
  {
    isPP = 0;
    run = 0;
    evt = 0;
    lumi = 0;
    hiBin = -1;
    vz = -99;
    weight = -1;
    pthat = -1;

    nTrk = 0;
    mult = 0;

    Dsize = -99;
    PVz = -99;
  }
  ~DTrackTree() {};

  DTrackTree(TTree* t) : DTrackTree() 
  {
    this->create_tree(t);
  }

  void create_tree(TTree* t);

  void set_hlt_tree(TTree* ht, Bool_t isPP);

  void clear_vectors();
  void copy_variables(DTree dt);
  void copy_variables_gen(GTree gt);
  void copy_index(DTree dt, int i);
  void copy_index_gen(GTree gt, int i);

  int isPP;
  uint32_t run;
  unsigned long long evt;
  uint32_t lumi;
  int hiBin;
  float vz;
  float weight;
  float pthat;

  float hiEvtPlanes[29];

  //

  float PVz;
  int Dsize;
  std::vector<int>      Dindex;
  std::vector<int>      Dtype;
  std::vector<float>    Dmass;
  std::vector<float>    Dpt;
  std::vector<float>    Deta;
  std::vector<float>    Dphi;
  std::vector<float>    Dy;
  std::vector<float>    Dchi2ndf;
  std::vector<float>    Dchi2cl;
  std::vector<float>    Ddtheta;
  std::vector<float>    Dalpha;
  std::vector<float>    DsvpvDistance;
  std::vector<float>    DsvpvDisErr;
  std::vector<float>    DlxyBS;
  std::vector<float>    DlxyBSErr;
  std::vector<float>    Dtrk1Pt;
  std::vector<float>    Dtrk2Pt;
  std::vector<float>    Dtrk1Eta;
  std::vector<float>    Dtrk2Eta;
  std::vector<float>    Dtrk1Phi;
  std::vector<float>    Dtrk2Phi;
  std::vector<float>    Dtrk1PtErr;
  std::vector<float>    Dtrk2PtErr;
  std::vector<float>    Dtrk1MassHypo;
  std::vector<float>    Dtrk2MassHypo;
  std::vector<bool>     Dtrk1highPurity;
  std::vector<bool>     Dtrk2highPurity;
  std::vector<int>      Dtrk1Idx;
  std::vector<int>      Dtrk2Idx;
  std::vector<float>    Dtrk1Y;
  std::vector<float>    Dtrk2Y;
  std::vector<float>    Dgen;
  std::vector<int>      DgenIndex;
  std::vector<int>      DgennDa;
  std::vector<float>    Dgenpt;
  std::vector<float>    Dgeneta;
  std::vector<float>    Dgenphi;
  std::vector<float>    Dgeny;
  std::vector<int>      DgencollisionId;
  std::vector<float>    DgenBAncestorpt;
  std::vector<int>      DgenBAncestorpdgId;

  int Gsize;
  std::vector<float>    Gy;
  std::vector<float>    Geta;
  std::vector<float>    Gphi;
  std::vector<float>    Gpt;
  std::vector<int>      GpdgId;
  std::vector<int>      GcollisionId;
  std::vector<int>      GisSignal;
  std::vector<float>    GBAncestorpt;
  std::vector<int>      GBAncestorpdgId;
  std::vector<float>    Gtk1pt;
  std::vector<float>    Gtk1eta;
  std::vector<float>    Gtk1y;
  std::vector<float>    Gtk1phi;
  std::vector<float>    Gtk2pt;
  std::vector<float>    Gtk2eta;
  std::vector<float>    Gtk2y;
  std::vector<float>    Gtk2phi;

  int nTrk;
  std::vector<float>    trkPt;
  std::vector<float>    trkPtError;
  std::vector<uint8_t>  trkNHit;
  std::vector<uint8_t>  trkNlayer;
  std::vector<float>    trkEta;
  std::vector<float>    trkPhi;
  std::vector<int>      trkCharge;
  std::vector<bool>     highPurity;
  std::vector<float>    trkChi2;
  std::vector<uint8_t>  trkNdof;
  std::vector<float>    trkDxy1;
  std::vector<float>    trkDxyError1;
  std::vector<float>    trkDz1;
  std::vector<float>    trkDzError1;
  std::vector<float>    pfEcal;
  std::vector<float>    pfHcal;
  std::vector<float>    trkWeight;

  int mult;
  std::vector<float>    pt;
  std::vector<float>    eta;
  std::vector<float>    phi;
  std::vector<int>      pdg;
  std::vector<int>      chg;
  std::vector<int>      matchingID;
  std::vector<int>      sube;

};

void DTrackTree::create_tree(TTree* t) 
{
  t->Branch("isPP", &isPP, "isPP/I");
  t->Branch("run", &run, "run/i");
  t->Branch("evt", &evt, "evt/l");
  t->Branch("lumi", &lumi, "lumi/i");
  t->Branch("hiBin", &hiBin, "hiBin/I");
  t->Branch("vz", &vz, "vz/F");
  t->Branch("weight", &weight, "weight/F");
  t->Branch("pthat", &pthat, "pthat/F");

  t->Branch("hiEvtPlanes", hiEvtPlanes, "hiEvtPlanes[29]/F");

  t->Branch("PVz", &PVz, "PVz/F");
  t->Branch("Dsize", &Dsize, "Dsize/I");
  t->Branch("Dindex", &Dindex);
  t->Branch("Dtype", &Dtype);
  t->Branch("Dmass", &Dmass);
  t->Branch("Dpt", &Dpt);
  t->Branch("Deta", &Deta);
  t->Branch("Dphi", &Dphi);
  t->Branch("Dy", &Dy);
  t->Branch("Dchi2ndf", &Dchi2ndf);
  t->Branch("Dchi2cl", &Dchi2cl);
  t->Branch("Ddtheta", &Ddtheta);
  t->Branch("Dalpha", &Dalpha);
  t->Branch("DsvpvDistance", &DsvpvDistance);
  t->Branch("DsvpvDisErr", &DsvpvDisErr);
  t->Branch("DlxyBS", &DlxyBS);
  t->Branch("DlxyBSErr", &DlxyBSErr);
  t->Branch("Dtrk1Pt", &Dtrk1Pt);
  t->Branch("Dtrk2Pt", &Dtrk2Pt);
  t->Branch("Dtrk1Eta", &Dtrk1Eta);
  t->Branch("Dtrk2Eta", &Dtrk2Eta);
  t->Branch("Dtrk1Phi", &Dtrk1Phi);
  t->Branch("Dtrk2Phi", &Dtrk2Phi);
  t->Branch("Dtrk1PtErr", &Dtrk1PtErr);
  t->Branch("Dtrk2PtErr", &Dtrk2PtErr);
  t->Branch("Dtrk1MassHypo", &Dtrk1MassHypo);
  t->Branch("Dtrk2MassHypo", &Dtrk2MassHypo);
  t->Branch("Dtrk1highPurity", &Dtrk1highPurity);
  t->Branch("Dtrk2highPurity", &Dtrk2highPurity);
  t->Branch("Dtrk1Idx", &Dtrk1Idx);
  t->Branch("Dtrk2Idx", &Dtrk2Idx);
  t->Branch("Dtrk1Y", &Dtrk1Y);
  t->Branch("Dtrk2Y", &Dtrk2Y);
  t->Branch("Dgen", &Dgen);
  t->Branch("DgenIndex", &DgenIndex);
  t->Branch("DgennDa", &DgennDa);
  t->Branch("Dgenpt", &Dgenpt);
  t->Branch("Dgeneta", &Dgeneta);
  t->Branch("Dgenphi", &Dgenphi);
  t->Branch("Dgeny", &Dgeny);
  t->Branch("DgencollisionId", &DgencollisionId);
  t->Branch("DgenBAncestorpt", &DgenBAncestorpt);
  t->Branch("DgenBAncestorpdgId", &DgenBAncestorpdgId);

  t->Branch("Gsize",&Gsize,"Gsize/I");
  t->Branch("Gy",&Gy);
  t->Branch("Geta",&Geta);
  t->Branch("Gphi",&Gphi);
  t->Branch("Gpt",&Gpt);
  t->Branch("GpdgId",&GpdgId);
  t->Branch("GcollisionId",&GcollisionId);
  t->Branch("GisSignal",&GisSignal);
  t->Branch("GBAncestorpt",&GBAncestorpt);
  t->Branch("GBAncestorpdgId",&GBAncestorpdgId);
  t->Branch("Gtk1pt",&Gtk1pt);
  t->Branch("Gtk1eta",&Gtk1eta);
  t->Branch("Gtk1y",&Gtk1y);
  t->Branch("Gtk1phi",&Gtk1phi);
  t->Branch("Gtk2pt",&Gtk2pt);
  t->Branch("Gtk2eta",&Gtk2eta);
  t->Branch("Gtk2y",&Gtk2y);
  t->Branch("Gtk2phi",&Gtk2phi);

  t->Branch("nTrk", &nTrk, "nTrk/I");
  t->Branch("trkPt", &trkPt);
  t->Branch("trkEta", &trkEta);
  t->Branch("trkPhi", &trkPhi);
  t->Branch("trkCharge", &trkCharge);
  t->Branch("trkPtError", &trkPtError);
  t->Branch("trkNHit", &trkNHit);
  t->Branch("trkNlayer", &trkNlayer);
  t->Branch("highPurity", &highPurity);
  t->Branch("trkChi2", &trkChi2);
  t->Branch("trkNdof", &trkNdof);
  t->Branch("trkDxy1", &trkDxy1);
  t->Branch("trkDxyError1", &trkDxyError1);
  t->Branch("trkDz1", &trkDz1);
  t->Branch("trkDzError1", &trkDzError1);
  t->Branch("pfEcal", &pfEcal);
  t->Branch("pfHcal", &pfHcal);
  t->Branch("trkWeight", &trkWeight);

  t->Branch("mult", &mult, "mult/I");
  t->Branch("pt", &pt);
  t->Branch("eta", &eta);
  t->Branch("phi", &phi);
  t->Branch("pdg", &pdg);
  t->Branch("chg", &chg);
  t->Branch("matchingID", &matchingID);
  t->Branch("sube", &sube);

}

void DTrackTree::copy_variables_gen(GTree gt)
{
  Gsize = gt.Gsize;
}

void DTrackTree::copy_variables(DTree dt) 
{
  PVz = dt.PVz;
  Dsize = dt.Dsize;
}

void DTrackTree::copy_index_gen(GTree gt, int i)
{
  Gy.push_back(gt.Gy[i]);
  Geta.push_back(gt.Geta[i]);
  Gphi.push_back(gt.Gphi[i]);
  Gpt.push_back(gt.Gpt[i]);
  GpdgId.push_back(gt.GpdgId[i]);
  GcollisionId.push_back(gt.GcollisionId[i]);
  GisSignal.push_back(gt.GisSignal[i]);
  GBAncestorpt.push_back(gt.GBAncestorpt[i]);
  GBAncestorpdgId.push_back(gt.GBAncestorpdgId[i]);
  Gtk1pt.push_back(gt.Gtk1pt[i]);
  Gtk1eta.push_back(gt.Gtk1eta[i]);
  Gtk1y.push_back(gt.Gtk1y[i]);
  Gtk1phi.push_back(gt.Gtk1phi[i]);
  Gtk2pt.push_back(gt.Gtk2pt[i]);
  Gtk2eta.push_back(gt.Gtk2eta[i]);
  Gtk2y.push_back(gt.Gtk2y[i]);
  Gtk2phi.push_back(gt.Gtk2phi[i]);
}

void DTrackTree::copy_index(DTree dt, int i) 
{
  Dindex.push_back(dt.Dindex[i]);
  Dtype.push_back(dt.Dtype[i]);
  Dmass.push_back(dt.Dmass[i]);
  Dpt.push_back(dt.Dpt[i]);
  Deta.push_back(dt.Deta[i]);
  Dphi.push_back(dt.Dphi[i]);
  Dy.push_back(dt.Dy[i]);
  Dchi2ndf.push_back(dt.Dchi2ndf[i]);
  Dchi2cl.push_back(dt.Dchi2cl[i]);
  Ddtheta.push_back(dt.Ddtheta[i]);
  Dalpha.push_back(dt.Dalpha[i]);
  DsvpvDistance.push_back(dt.DsvpvDistance[i]);
  DsvpvDisErr.push_back(dt.DsvpvDisErr[i]);
  DlxyBS.push_back(dt.DlxyBS[i]);
  DlxyBSErr.push_back(dt.DlxyBSErr[i]);
  Dtrk1Pt.push_back(dt.Dtrk1Pt[i]);
  Dtrk2Pt.push_back(dt.Dtrk2Pt[i]);
  Dtrk1Eta.push_back(dt.Dtrk1Eta[i]);
  Dtrk2Eta.push_back(dt.Dtrk2Eta[i]);
  Dtrk1Phi.push_back(dt.Dtrk1Phi[i]);
  Dtrk2Phi.push_back(dt.Dtrk2Phi[i]);
  Dtrk1PtErr.push_back(dt.Dtrk1PtErr[i]);
  Dtrk2PtErr.push_back(dt.Dtrk2PtErr[i]);
  Dtrk1MassHypo.push_back(dt.Dtrk1MassHypo[i]);
  Dtrk2MassHypo.push_back(dt.Dtrk2MassHypo[i]);
  Dtrk1highPurity.push_back(dt.Dtrk1highPurity[i]);
  Dtrk2highPurity.push_back(dt.Dtrk2highPurity[i]);
  Dtrk1Idx.push_back(dt.Dtrk1Idx[i]);
  Dtrk2Idx.push_back(dt.Dtrk2Idx[i]);
  Dtrk1Y.push_back(dt.Dtrk1Y[i]);
  Dtrk2Y.push_back(dt.Dtrk2Y[i]);
  Dgen.push_back(dt.Dgen[i]);
  DgenIndex.push_back(dt.DgenIndex[i]);
  DgennDa.push_back(dt.DgennDa[i]);
  Dgenpt.push_back(dt.Dgenpt[i]);
  Dgeneta.push_back(dt.Dgeneta[i]);
  Dgenphi.push_back(dt.Dgenphi[i]);
  Dgeny.push_back(dt.Dgeny[i]);
  DgencollisionId.push_back(dt.DgencollisionId[i]);
  DgenBAncestorpt.push_back(dt.DgenBAncestorpt[i]);
  DgenBAncestorpdgId.push_back(dt.DgenBAncestorpdgId[i]);
}

void DTrackTree::clear_vectors() 
{
  Dindex.clear();
  Dtype.clear();
  Dmass.clear();
  Dpt.clear();
  Deta.clear();
  Dphi.clear();
  Dy.clear();
  Dchi2ndf.clear();
  Dchi2cl.clear();
  Ddtheta.clear();
  Dalpha.clear();
  DsvpvDistance.clear();
  DsvpvDisErr.clear();
  DlxyBS.clear();
  DlxyBSErr.clear();
  Dtrk1Pt.clear();
  Dtrk2Pt.clear();
  Dtrk1Eta.clear();
  Dtrk2Eta.clear();
  Dtrk1Phi.clear();
  Dtrk2Phi.clear();
  Dtrk1PtErr.clear();
  Dtrk2PtErr.clear();
  Dtrk1MassHypo.clear();
  Dtrk2MassHypo.clear();
  Dtrk1highPurity.clear();
  Dtrk2highPurity.clear();
  Dtrk1Idx.clear();
  Dtrk2Idx.clear();
  Dtrk1Y.clear();
  Dtrk2Y.clear();
  Dgen.clear();
  DgenIndex.clear();
  DgennDa.clear();
  Dgenpt.clear();
  Dgeneta.clear();
  Dgenphi.clear();
  Dgeny.clear();
  DgencollisionId.clear();
  DgenBAncestorpt.clear();
  DgenBAncestorpdgId.clear();

  Gy.clear();
  Geta.clear();
  Gphi.clear();
  Gpt.clear();
  GpdgId.clear();
  GcollisionId.clear();
  GisSignal.clear();
  GBAncestorpt.clear();
  GBAncestorpdgId.clear();
  Gtk1pt.clear();
  Gtk1eta.clear();
  Gtk1y.clear();
  Gtk1phi.clear();
  Gtk2pt.clear();
  Gtk2eta.clear();
  Gtk2y.clear();
  Gtk2phi.clear();

  trkPt.clear();
  trkPtError.clear();
  trkNHit.clear();
  trkNlayer.clear();
  trkEta.clear();
  trkPhi.clear();
  trkCharge.clear();
  highPurity.clear();
  trkChi2.clear();
  trkNdof.clear();
  trkDxy1.clear();
  trkDxyError1.clear();
  trkDz1.clear();
  trkDzError1.clear();
  pfEcal.clear();
  pfHcal.clear();
  trkWeight.clear();

  pt.clear();
  eta.clear();
  phi.clear();
  pdg.clear();
  chg.clear();
  matchingID.clear();
  sube.clear();

}

void DTrackTree::set_hlt_tree(TTree* ht, Bool_t isPP)
{
  /*
    ht->SetBranchStatus("*", 0);
    if(isPP) //pp
    {
    ht->SetBranchStatus("HLT_L1MinimumBias*",1);
    ht->SetBranchStatus("HLT_Dmeson*",1);
    ht->SetBranchStatus("HLT_AK4CaloJet40_Eta5p1_v1*",1);
    ht->SetBranchStatus("HLT_AK4CaloJet60_Eta5p1_v1*",1);
    ht->SetBranchStatus("HLT_AK4CaloJet80_Eta5p1_v1*",1);
    ht->SetBranchStatus("L1_SingleJet16_BptxAND*",1);
    ht->SetBranchStatus("L1_SingleJet24_BptxAND*",1);
    ht->SetBranchStatus("L1_SingleJet28_BptxAND*",1);
    ht->SetBranchStatus("L1_SingleJet40_BptxAND*",1);
    ht->SetBranchStatus("L1_SingleJet48_BptxAND*",1);      
    }
    else //PbPb
    {
    ht->SetBranchStatus("HLT_HIL1MinimumBias*",1);
    ht->SetBranchStatus("HLT_HIDmeson*",1);
    ht->SetBranchStatus("HLT_HIPuAK4CaloJet*",1);
    ht->SetBranchStatus("L1_MinimumBiasHF1_AND",1);
    ht->SetBranchStatus("L1_MinimumBiasHF2_AND",1);
    ht->SetBranchStatus("L1_SingleS1Jet16_BptxAND*", 1);
    ht->SetBranchStatus("L1_SingleS1Jet28_BptxAND*", 1);
    ht->SetBranchStatus("L1_SingleS1Jet32_BptxAND*", 1);
    ht->SetBranchStatus("L1_SingleJet44_BptxAND*",   1);
    ht->SetBranchStatus("L1_SingleS1Jet52_BptxAND*", 1);
    ht->SetBranchStatus("L1_SingleS1Jet56_BptxAND*", 1);
    ht->SetBranchStatus("L1_SingleS1Jet64_BptxAND*", 1);
    }
  */
}

#endif
