#ifndef _TRACK_TREE_H
#define _TRACK_TREE_H

#include "TTree.h"

class trackTree 
{
public:
  trackTree() 
  {
    nTrk = -99;
    trkPt = 0;
    trkEta = 0;
    trkPhi = 0;
    trkCharge = 0;
    trkPtError = 0;
    trkNHit = 0;
    trkNlayer = 0;
    highPurity = 0;
    trkChi2 = 0;
    trkNdof = 0;
    trkDxy1 = 0;
    trkDxyError1 = 0;
    trkDz1 = 0;
    trkDzError1 = 0;
    pfEcal = 0;
    pfHcal = 0;
  }
  ~trackTree() 
  {
    clear_memory();
  };

  trackTree(TTree* t) : trackTree() 
  {
    this->read_tree(t);
  }

  void init_memory();
  void clear_memory();
  void read_tree(TTree* t);

  int nTrk;
  float* trkPt;
  float* trkEta;
  float* trkPhi;
  int* trkCharge;
  float* trkPtError;
  uint8_t* trkNHit;
  uint8_t* trkNlayer;
  bool* highPurity;
  float* trkChi2;
  uint8_t* trkNdof;
  float* trkDxy1;
  float* trkDxyError1;
  float* trkDz1;
  float* trkDzError1;
  float* pfEcal;
  float* pfHcal;
};

void trackTree::init_memory() 
{
  trkPt = (float*) malloc(sizeof(float) * 100000);
  trkEta = (float*) malloc(sizeof(float) * 100000);
  trkPhi = (float*) malloc(sizeof(float) * 100000);
  trkCharge = (int*) malloc(sizeof(int) * 100000);
  trkPtError = (float*) malloc(sizeof(float) * 100000);
  trkNHit = (uint8_t*) malloc(sizeof(uint8_t) * 100000);
  trkNlayer = (uint8_t*) malloc(sizeof(uint8_t) * 100000);
  highPurity = (bool*) malloc(sizeof(bool) * 100000);
  trkChi2 = (float*) malloc(sizeof(float) * 100000);
  trkNdof = (uint8_t*) malloc(sizeof(uint8_t) * 100000);
  trkDxy1 = (float*) malloc(sizeof(float) * 100000);
  trkDxyError1 = (float*) malloc(sizeof(float) * 100000);
  trkDz1 = (float*) malloc(sizeof(float) * 100000);
  trkDzError1 = (float*) malloc(sizeof(float) * 100000);
  pfEcal = (float*) malloc(sizeof(float) * 100000);
  pfHcal = (float*) malloc(sizeof(float) * 100000);
}

void trackTree::clear_memory() 
{
  free(trkPt);
  free(trkEta);
  free(trkPhi);
  free(trkCharge);
  free(trkPtError);
  free(trkNHit);
  free(trkNlayer);
  free(highPurity);
  free(trkChi2);
  free(trkNdof);
  free(trkDxy1);
  free(trkDxyError1);
  free(trkDz1);
  free(trkDzError1);
  free(pfEcal);
  free(pfHcal);
}

void trackTree::read_tree(TTree* t) 
{
  init_memory();

  t->SetBranchStatus("nTrk", 1);
  t->SetBranchStatus("trkPt", 1);
  t->SetBranchStatus("trkEta", 1);
  t->SetBranchStatus("trkPhi", 1);
  t->SetBranchStatus("trkCharge", 1);
  t->SetBranchStatus("trkPtError", 1);
  t->SetBranchStatus("trkNHit", 1);
  t->SetBranchStatus("trkNlayer", 1);
  t->SetBranchStatus("highPurity", 1);
  t->SetBranchStatus("trkChi2", 1);
  t->SetBranchStatus("trkNdof", 1);
  t->SetBranchStatus("trkDxy1", 1);
  t->SetBranchStatus("trkDxyError1", 1);
  t->SetBranchStatus("trkDz1", 1);
  t->SetBranchStatus("trkDzError1", 1);
  t->SetBranchStatus("pfEcal", 1);
  t->SetBranchStatus("pfHcal", 1);

  t->SetBranchAddress("nTrk", &nTrk);
  t->SetBranchAddress("trkPt", trkPt);
  t->SetBranchAddress("trkEta", trkEta);
  t->SetBranchAddress("trkPhi", trkPhi);
  t->SetBranchAddress("trkCharge", trkCharge);
  t->SetBranchAddress("trkPtError", trkPtError);
  t->SetBranchAddress("trkNHit", trkNHit);
  t->SetBranchAddress("trkNlayer", trkNlayer);
  t->SetBranchAddress("highPurity", highPurity);
  t->SetBranchAddress("trkChi2", trkChi2);
  t->SetBranchAddress("trkNdof", trkNdof);
  t->SetBranchAddress("trkDxy1", trkDxy1);
  t->SetBranchAddress("trkDxyError1", trkDxyError1);
  t->SetBranchAddress("trkDz1", trkDz1);
  t->SetBranchAddress("trkDzError1", trkDzError1);
  t->SetBranchAddress("pfEcal", pfEcal);
  t->SetBranchAddress("pfHcal", pfHcal);
}

#endif
