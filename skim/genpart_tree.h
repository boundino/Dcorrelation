#ifndef _GENPART_TREE
#define _GENPART_TREE

#include "TTree.h"

class genpartTree 
{
public:
  genpartTree() 
  {
    mult = -99;
    pt = 0;
    eta = 0;
    phi = 0;
    pdg = 0;
    chg = 0;
    matchingID = 0;
    sube = 0;
  }
  ~genpartTree() {};

  genpartTree(TTree* t) : genpartTree() 
  {
    this->read_tree(t);
  }

  void read_tree(TTree* t);

  int mult;
  std::vector<float>* pt;
  std::vector<float>* eta;
  std::vector<float>* phi;
  std::vector<int>* pdg;
  std::vector<int>* chg;
  std::vector<int>* matchingID;
  std::vector<int>* sube;
};

void genpartTree::read_tree(TTree* t) 
{
  t->SetBranchStatus("mult", 1);
  t->SetBranchStatus("pt", 1);
  t->SetBranchStatus("eta", 1);
  t->SetBranchStatus("phi", 1);
  t->SetBranchStatus("pdg", 1);
  t->SetBranchStatus("chg", 1);
  t->SetBranchStatus("matchingID", 1);
  t->SetBranchStatus("sube", 1);

  t->SetBranchAddress("mult", &mult);
  t->SetBranchAddress("pt", &pt);
  t->SetBranchAddress("eta", &eta);
  t->SetBranchAddress("phi", &phi);
  t->SetBranchAddress("pdg", &pdg);
  t->SetBranchAddress("chg", &chg);
  t->SetBranchAddress("matchingID", &matchingID);
  t->SetBranchAddress("sube", &sube);
}

#endif
