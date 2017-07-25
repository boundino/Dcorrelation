#ifndef _READD_H
#define _READD_H

#include <TTree.h>

class readD
{
 private:

  int              MAX_ND;
  int              MAX_NG;
  float            cut_trkPt;
  float            cut_trkEta;
  float            cut_trkPtErr;
  float            cut_Dy;
  float            cut_Dsvpv;
  float            cut_Dalpha;
  float            cut_Dchi2cl;
  float            cut_Dpt_min;
  float            cut_Dpt_max;

  void init();
  void setbranchaddress(TTree* nt, const char* bname, void* addr);

public:

  int              hiBin;
  float            vz;
  float            weight;
  float            pthatweight;
  float            maxDgenpt;
  float            pthat;
  int              RunNo;
  int              EvtNo;
  int              LumiNo;
  int              Dsize;
  float            PVx;
  float            PVy;
  float            PVz;
  float            PVnchi2;

  int*             Dindex;
  float*           Dmass;
  float*           Dpt;
  float*           Deta;
  float*           Dphi;
  float*           Dy;
  float*           Dchi2ndf;
  float*           Dchi2cl;
  float*           Ddtheta;
  float*           Dlxy;
  float*           Dalpha;
  float*           DsvpvDistance;
  float*           DsvpvDisErr;
  float*           DsvpvDistance_2D;
  float*           DsvpvDisErr_2D;
  float*           DlxyBS;
  float*           DlxyBSErr;
  float*           DMaxDoca;
  float*           Dtrk1Pt;
  float*           Dtrk2Pt;
  float*           Dtrk1Eta;
  float*           Dtrk2Eta;
  float*           Dtrk1Phi;
  float*           Dtrk2Phi;
  float*           Dtrk1PtErr;
  float*           Dtrk2PtErr;
  float*           Dtrk1Dxy;
  float*           Dtrk2Dxy;
  float*           Dtrk1PixelHit;
  float*           Dtrk2PixelHit;
  float*           Dtrk1StripHit;
  float*           Dtrk2StripHit;
  float*           Dtrk1nStripLayer;
  float*           Dtrk2nStripLayer;
  float*           Dtrk1nPixelLayer;
  float*           Dtrk2nPixelLayer;
  float*           Dtrk1Chi2ndf;
  float*           Dtrk2Chi2ndf;
  float*           Dtrk1MassHypo;
  float*           Dtrk2MassHypo;
  int*             Dtrk1Algo;
  int*             Dtrk2Algo;
  bool*            Dtrk1highPurity;
  bool*            Dtrk2highPurity;
  float*           Dtrk1EtaErr;
  float*           Dtrk2EtaErr;
  float*           Dtrk1PhiErr;
  float*           Dtrk2PhiErr;
  float*           Dtrk1Y;
  float*           Dtrk2Y;

  float*           Dgen;

  int              Gsize;
  int*             GpdgId;
  int*             GisSignal;
  float*           Gy;
  float*           Geta;
  float*           Gphi;
  float*           Gpt;

  //
  readD(int MAX_ND_, int MAX_NG_)
    {
      MAX_ND = MAX_ND_;
      MAX_NG = MAX_NG_;
      cut_trkPt = -99;
      cut_trkEta = -99;
      cut_trkPtErr = -99;
      cut_Dy = -99;
      cut_Dsvpv = -99;
      cut_Dalpha = -99;
      cut_Dchi2cl = -99;
      cut_Dpt_min = -99;
      cut_Dpt_max = -99;
    }
  ~readD() {};

  void setbranchesaddress(TTree* nt, TTree* ntGen, TTree* ntHi);

  void settrkcut(float _cut_trkPt, float _cut_trkEta, float _cut_trkPtErr);
  void setDcut(float _cut_Dy, float _cut_Dsvpv, float _cut_Dalpha, float _cut_Dchi2cl, float _cut_Dpt_min=-1, float _cut_Dpt_max=1.e+10);
  void settmvacut(float _cut_Dchi2cl, float _cut_Dsvpv); // update
  bool isselected(int j);
};


void readD::init()
{
  hiBin = -99;
  vz = -99;
  weight = -99;
  pthat = -99;
  pthatweight = -99;
  maxDgenpt = -99;
  RunNo = -99;
  EvtNo = -99;
  LumiNo = -99;
  Dsize = -99;
  PVx = -99;
  PVy = -99;
  PVz = -99;
  PVnchi2 = -99;
  Gsize = -99;

  //
  Dindex = new int[MAX_ND];
  Dmass = new float[MAX_ND];
  Dpt = new float[MAX_ND];
  Deta = new float[MAX_ND];
  Dphi = new float[MAX_ND];
  Dy = new float[MAX_ND];
  Dchi2ndf = new float[MAX_ND];
  Dchi2cl = new float[MAX_ND];
  Ddtheta = new float[MAX_ND];
  Dlxy = new float[MAX_ND];
  Dalpha = new float[MAX_ND];
  DsvpvDistance = new float[MAX_ND];
  DsvpvDisErr = new float[MAX_ND];
  DsvpvDistance_2D = new float[MAX_ND];
  DsvpvDisErr_2D = new float[MAX_ND];
  DlxyBS = new float[MAX_ND];
  DlxyBSErr = new float[MAX_ND];
  DMaxDoca = new float[MAX_ND];
  Dtrk1Pt = new float[MAX_ND];
  Dtrk2Pt = new float[MAX_ND];
  Dtrk1Eta = new float[MAX_ND];
  Dtrk2Eta = new float[MAX_ND];
  Dtrk1Phi = new float[MAX_ND];
  Dtrk2Phi = new float[MAX_ND];
  Dtrk1PtErr = new float[MAX_ND];
  Dtrk2PtErr = new float[MAX_ND];
  Dtrk1Dxy = new float[MAX_ND];
  Dtrk2Dxy = new float[MAX_ND];
  Dtrk1PixelHit = new float[MAX_ND];
  Dtrk2PixelHit = new float[MAX_ND];
  Dtrk1StripHit = new float[MAX_ND];
  Dtrk2StripHit = new float[MAX_ND];
  Dtrk1nStripLayer = new float[MAX_ND];
  Dtrk2nStripLayer = new float[MAX_ND];
  Dtrk1nPixelLayer = new float[MAX_ND];
  Dtrk2nPixelLayer = new float[MAX_ND];
  Dtrk1Chi2ndf = new float[MAX_ND];
  Dtrk2Chi2ndf = new float[MAX_ND];
  Dtrk1MassHypo = new float[MAX_ND];
  Dtrk2MassHypo = new float[MAX_ND];
  Dtrk1Algo = new int[MAX_ND];
  Dtrk2Algo = new int[MAX_ND];
  Dtrk1highPurity = new bool[MAX_ND];
  Dtrk2highPurity = new bool[MAX_ND];
  Dtrk1EtaErr = new float[MAX_ND];
  Dtrk2EtaErr = new float[MAX_ND];
  Dtrk1PhiErr = new float[MAX_ND];
  Dtrk2PhiErr = new float[MAX_ND];
  Dtrk1Y = new float[MAX_ND];
  Dtrk2Y = new float[MAX_ND];

  Dgen = new float[MAX_ND];

  GpdgId = new int[MAX_NG];
  GisSignal = new int[MAX_NG];
  Gy = new float[MAX_NG];
  Geta = new float[MAX_NG];
  Gphi = new float[MAX_NG];
  Gpt = new float[MAX_NG];
}

void readD::setbranchaddress(TTree* nt, const char* bname, void* addr)
{
  nt->SetBranchStatus(bname, 1);
  nt->SetBranchAddress(bname, addr);
}

void readD::setbranchesaddress(TTree* nt, TTree* ntGen, TTree* ntHi)
{
  init();
  nt->SetBranchStatus("*", 0);
  ntGen->SetBranchStatus("*", 0);
  ntHi->SetBranchStatus("*", 0);

  setbranchaddress(ntHi,  "hiBin",            &hiBin);
  setbranchaddress(ntHi,  "vz",               &vz);
  setbranchaddress(ntHi,  "weight",           &weight);
  setbranchaddress(ntHi,  "pthat",            &pthat);

  setbranchaddress(ntHi,  "pthatweight",      &pthatweight);
  setbranchaddress(ntHi,  "maxDgenpt",        &maxDgenpt);

  setbranchaddress(nt,    "RunNo",            &RunNo);
  setbranchaddress(nt,    "EvtNo",            &EvtNo);
  setbranchaddress(nt,    "LumiNo",           &LumiNo);
  setbranchaddress(nt,    "Dsize",            &Dsize);
  setbranchaddress(nt,    "PVx",              &PVx);
  setbranchaddress(nt,    "PVy",              &PVy);
  setbranchaddress(nt,    "PVz",              &PVz);
  setbranchaddress(nt,    "PVnchi2",          &PVnchi2);

  setbranchaddress(nt,    "Dindex",           Dindex);
  setbranchaddress(nt,    "Dmass",            Dmass);
  setbranchaddress(nt,    "Dpt",              Dpt);
  setbranchaddress(nt,    "Deta",             Deta);
  setbranchaddress(nt,    "Dphi",             Dphi);
  setbranchaddress(nt,    "Dy",               Dy);
  setbranchaddress(nt,    "Dchi2ndf",         Dchi2ndf);
  setbranchaddress(nt,    "Dchi2cl",          Dchi2cl);
  setbranchaddress(nt,    "Ddtheta",          Ddtheta);
  setbranchaddress(nt,    "Dlxy",             Dlxy);
  setbranchaddress(nt,    "Dalpha",           Dalpha);
  setbranchaddress(nt,    "DsvpvDistance",    DsvpvDistance);
  setbranchaddress(nt,    "DsvpvDisErr",      DsvpvDisErr);
  setbranchaddress(nt,    "DsvpvDistance_2D", DsvpvDistance_2D);
  setbranchaddress(nt,    "DsvpvDisErr_2D",   DsvpvDisErr_2D);
  setbranchaddress(nt,    "DlxyBS",           DlxyBS);
  setbranchaddress(nt,    "DlxyBSErr",        DlxyBSErr);
  setbranchaddress(nt,    "DMaxDoca",         DMaxDoca);
  setbranchaddress(nt,    "Dtrk1Pt",          Dtrk1Pt);
  setbranchaddress(nt,    "Dtrk2Pt",          Dtrk2Pt);
  setbranchaddress(nt,    "Dtrk1Eta",         Dtrk1Eta);
  setbranchaddress(nt,    "Dtrk2Eta",         Dtrk2Eta);
  setbranchaddress(nt,    "Dtrk1Phi",         Dtrk1Phi);
  setbranchaddress(nt,    "Dtrk2Phi",         Dtrk2Phi);
  setbranchaddress(nt,    "Dtrk1PtErr",       Dtrk1PtErr);
  setbranchaddress(nt,    "Dtrk2PtErr",       Dtrk2PtErr);
  setbranchaddress(nt,    "Dtrk1Dxy",         Dtrk1Dxy);
  setbranchaddress(nt,    "Dtrk2Dxy",         Dtrk2Dxy);
  setbranchaddress(nt,    "Dtrk1PixelHit",    Dtrk1PixelHit);
  setbranchaddress(nt,    "Dtrk2PixelHit",    Dtrk2PixelHit);
  setbranchaddress(nt,    "Dtrk1StripHit",    Dtrk1StripHit);
  setbranchaddress(nt,    "Dtrk2StripHit",    Dtrk2StripHit);
  setbranchaddress(nt,    "Dtrk1nStripLayer", Dtrk1nStripLayer);
  setbranchaddress(nt,    "Dtrk2nStripLayer", Dtrk2nStripLayer);
  setbranchaddress(nt,    "Dtrk1nPixelLayer", Dtrk1nPixelLayer);
  setbranchaddress(nt,    "Dtrk2nPixelLayer", Dtrk2nPixelLayer);
  setbranchaddress(nt,    "Dtrk1Chi2ndf",     Dtrk1Chi2ndf);
  setbranchaddress(nt,    "Dtrk2Chi2ndf",     Dtrk2Chi2ndf);
  setbranchaddress(nt,    "Dtrk1MassHypo",    Dtrk1MassHypo);
  setbranchaddress(nt,    "Dtrk2MassHypo",    Dtrk2MassHypo);
  setbranchaddress(nt,    "Dtrk1Algo",        Dtrk1Algo);
  setbranchaddress(nt,    "Dtrk2Algo",        Dtrk2Algo);
  setbranchaddress(nt,    "Dtrk1highPurity",  Dtrk1highPurity);
  setbranchaddress(nt,    "Dtrk2highPurity",  Dtrk2highPurity);
  setbranchaddress(nt,    "Dtrk1EtaErr",      Dtrk1EtaErr);
  setbranchaddress(nt,    "Dtrk2EtaErr",      Dtrk2EtaErr);
  setbranchaddress(nt,    "Dtrk1PhiErr",      Dtrk1PhiErr);
  setbranchaddress(nt,    "Dtrk2PhiErr",      Dtrk2PhiErr);
  setbranchaddress(nt,    "Dtrk1Y",           Dtrk1Y);
  setbranchaddress(nt,    "Dtrk2Y",           Dtrk2Y);
  setbranchaddress(nt,    "Dgen",             Dgen);

  setbranchaddress(ntGen, "Gsize",            &Gsize);
  setbranchaddress(ntGen, "GpdgId",           GpdgId);
  setbranchaddress(ntGen, "GisSignal",        GisSignal);
  setbranchaddress(ntGen, "Gy",               Gy);
  setbranchaddress(ntGen, "Geta",             Geta);
  setbranchaddress(ntGen, "Gphi",             Gphi);
  setbranchaddress(ntGen, "Gpt",              Gpt);

}

void readD::settrkcut(float _cut_trkPt, float _cut_trkEta, float _cut_trkPtErr)
{
  cut_trkPt = _cut_trkPt;
  cut_trkEta = _cut_trkEta;
  cut_trkPtErr = _cut_trkPtErr;
}

void readD::setDcut(float _cut_Dy, float _cut_Dsvpv, float _cut_Dalpha, float _cut_Dchi2cl, float _cut_Dpt_min/*=-1*/, float _cut_Dpt_max/*=1.e+10*/)
{
  cut_Dy = _cut_Dy;
  cut_Dsvpv = _cut_Dsvpv;
  cut_Dalpha = _cut_Dalpha;
  cut_Dchi2cl = _cut_Dchi2cl;
  cut_Dpt_min = _cut_Dpt_min;
  cut_Dpt_max = _cut_Dpt_max;
}

void readD::settmvacut(float _cut_Dchi2cl, float _cut_Dsvpv)
{
  cut_Dchi2cl = _cut_Dchi2cl;
  cut_Dsvpv = _cut_Dsvpv;
}

bool readD::isselected(int j)
{
  if(Dtrk1Pt[j] > cut_trkPt && Dtrk2Pt[j] > cut_trkPt &&
     TMath::Abs(Dtrk1Eta[j]) < cut_trkEta && TMath::Abs(Dtrk2Eta[j]) < cut_trkEta &&
     (Dtrk1PtErr[j]/Dtrk1Pt[j]) < cut_trkPtErr && (Dtrk2PtErr[j]/Dtrk2Pt[j]) < cut_trkPtErr &&
     Dtrk1highPurity[j] && Dtrk2highPurity[j] &&
     TMath::Abs(Dy[j]) < cut_Dy &&
     (DsvpvDistance[j]/DsvpvDisErr[j]) > cut_Dsvpv &&
     Dalpha[j] < cut_Dalpha &&
     Dchi2cl[j] > cut_Dchi2cl &&
     Dpt[j] > cut_Dpt_min && Dpt[j] < cut_Dpt_max) return true;
  else return false;
}


#endif
