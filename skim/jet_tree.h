#ifndef _JET_TREE_H
#define _JET_TREE_H

#include "TTree.h"

class jetTree 
{
public:
  jetTree() 
  {
    nref = -99;
    rawpt = 0;
    jtpt = 0;
    jteta = 0;
    jtphi = 0;
    refpt = 0;
    refeta = 0;
    refphi = 0;
    refparton_flavor = 0;
    subid = 0;
    chargedSum = 0;

    ngen = -99;
    genpt = 0;
    geneta = 0;
    genphi = 0;
    gensubid = 0;
  }
  ~jetTree() 
  {
    clear_memory();
  };

  jetTree(TTree* t) : jetTree() 
  {
    this->read_tree(t);
  }

  void init_memory();
  void clear_memory();
  void read_tree(TTree* t);

  int nref;
  float* rawpt;
  float* jtpt;
  float* jteta;
  float* jtphi;
  float* refpt;
  float* refeta;
  float* refphi;
  int* refparton_flavor;
  int* subid;
  float* chargedSum;

  int ngen;
  float* genpt;
  float* geneta;
  float* genphi;
  int* gensubid;
};

void jetTree::init_memory() 
{
  rawpt = (float*) malloc(sizeof(float) * 300);
  jtpt = (float*) malloc(sizeof(float) * 300);
  jteta = (float*) malloc(sizeof(float) * 300);
  jtphi = (float*) malloc(sizeof(float) * 300);
  refpt = (float*) malloc(sizeof(float) * 300);
  refeta = (float*) malloc(sizeof(float) * 300);
  refphi = (float*) malloc(sizeof(float) * 300);
  refparton_flavor = (int*) malloc(sizeof(int) * 300);
  subid = (int*) malloc(sizeof(int) * 300);
  chargedSum = (float*) malloc(sizeof(float) * 300);

  genpt = (float*) malloc(sizeof(float) * 5000);
  geneta = (float*) malloc(sizeof(float) * 5000);
  genphi = (float*) malloc(sizeof(float) * 5000);
  gensubid = (int*) malloc(sizeof(int) * 5000);
}

void jetTree::clear_memory() 
{
  free(rawpt);
  free(jtpt);
  free(jteta);
  free(jtphi);
  free(refpt);
  free(refeta);
  free(refphi);
  free(refparton_flavor);
  free(subid);
  free(chargedSum);

  free(genpt);
  free(geneta);
  free(genphi);
  free(gensubid);
}

void jetTree::read_tree(TTree* t) 
{
  init_memory();

  t->SetBranchStatus("nref", 1);
  t->SetBranchStatus("rawpt", 1);
  t->SetBranchStatus("jtpt", 1);
  t->SetBranchStatus("jteta", 1);
  t->SetBranchStatus("jtphi", 1);
  t->SetBranchStatus("refpt", 1);
  t->SetBranchStatus("refeta", 1);
  t->SetBranchStatus("refphi", 1);
  t->SetBranchStatus("refparton_flavor", 1);
  t->SetBranchStatus("chargedSum", 1);
  t->SetBranchStatus("subid", 1);

  t->SetBranchStatus("ngen", 1);
  t->SetBranchStatus("genpt", 1);
  t->SetBranchStatus("geneta", 1);
  t->SetBranchStatus("genphi", 1);
  t->SetBranchStatus("gensubid", 1);

  t->SetBranchAddress("nref", &nref);
  t->SetBranchAddress("rawpt", rawpt);
  t->SetBranchAddress("jtpt", jtpt);
  t->SetBranchAddress("jteta", jteta);
  t->SetBranchAddress("jtphi", jtphi);
  t->SetBranchAddress("refpt", refpt);
  t->SetBranchAddress("refeta", refeta);
  t->SetBranchAddress("refphi", refphi);
  t->SetBranchAddress("refparton_flavor", refparton_flavor);
  t->SetBranchAddress("chargedSum", chargedSum);
  t->SetBranchAddress("subid", subid);

  t->SetBranchAddress("ngen", &ngen);
  t->SetBranchAddress("genpt", genpt);
  t->SetBranchAddress("geneta", geneta);
  t->SetBranchAddress("genphi", genphi);
  t->SetBranchAddress("gensubid", gensubid);
}

#endif
