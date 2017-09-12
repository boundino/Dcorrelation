#ifndef _G_TREE_H
#define _G_TREE_H

#include "TTree.h"
#define MAX_GEN      6000

class GTree
{
 public:
  GTree()
    {
      Gsize = -99;
    }
  ~GTree() {};
  
 GTree(TTree* t) : GTree()
    {
      this->read_tree(t);
    }
  void read_tree(TTree* t);

  Int_t     Gsize;
  Float_t   Gy[MAX_GEN];
  Float_t   Geta[MAX_GEN];
  Float_t   Gphi[MAX_GEN];
  Float_t   Gpt[MAX_GEN];
  Int_t     GpdgId[MAX_GEN];
  Int_t     GcollisionId[MAX_GEN];
  Int_t     GisSignal[MAX_GEN];
  Float_t   GBAncestorpt[MAX_GEN];
  Int_t     GBAncestorpdgId[MAX_GEN];
  Float_t   Gtk1pt[MAX_GEN];
  Float_t   Gtk1eta[MAX_GEN];
  Float_t   Gtk1y[MAX_GEN];
  Float_t   Gtk1phi[MAX_GEN];
  Float_t   Gtk2pt[MAX_GEN];
  Float_t   Gtk2eta[MAX_GEN];
  Float_t   Gtk2y[MAX_GEN];
  Float_t   Gtk2phi[MAX_GEN];
};

void GTree::read_tree(TTree* t)
{
  t->SetBranchStatus("Gsize", 1);
  t->SetBranchStatus("Gy", 1);
  t->SetBranchStatus("Geta", 1);
  t->SetBranchStatus("Gphi", 1);
  t->SetBranchStatus("Gpt", 1);
  t->SetBranchStatus("GpdgId", 1);
  t->SetBranchStatus("GcollisionId", 1);
  t->SetBranchStatus("GisSignal", 1);
  t->SetBranchStatus("GBAncestorpt", 1);
  t->SetBranchStatus("GBAncestorpdgId", 1);
  t->SetBranchStatus("Gtk1pt", 1);
  t->SetBranchStatus("Gtk1eta", 1);
  t->SetBranchStatus("Gtk1y", 1);
  t->SetBranchStatus("Gtk1phi", 1);
  t->SetBranchStatus("Gtk2pt", 1);
  t->SetBranchStatus("Gtk2eta", 1);
  t->SetBranchStatus("Gtk2y", 1);
  t->SetBranchStatus("Gtk2phi", 1);

  t->SetBranchAddress("Gsize",&Gsize);
  t->SetBranchAddress("Gy",Gy);
  t->SetBranchAddress("Geta",Geta);
  t->SetBranchAddress("Gphi",Gphi);
  t->SetBranchAddress("Gpt",Gpt);
  t->SetBranchAddress("GpdgId",GpdgId);
  t->SetBranchAddress("GcollisionId",GcollisionId);
  t->SetBranchAddress("GisSignal",GisSignal);
  t->SetBranchAddress("GBAncestorpt",GBAncestorpt);
  t->SetBranchAddress("GBAncestorpdgId",GBAncestorpdgId);
  t->SetBranchAddress("Gtk1pt",Gtk1pt);
  t->SetBranchAddress("Gtk1eta",Gtk1eta);
  t->SetBranchAddress("Gtk1y",Gtk1y);
  t->SetBranchAddress("Gtk1phi",Gtk1phi);
  t->SetBranchAddress("Gtk2pt",Gtk2pt);
  t->SetBranchAddress("Gtk2eta",Gtk2eta);
  t->SetBranchAddress("Gtk2y",Gtk2y);
  t->SetBranchAddress("Gtk2phi",Gtk2phi);
}

#endif
