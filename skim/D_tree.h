#ifndef _D_TREE_H
#define _D_TREE_H

#include "TTree.h"
#define MAX_XB       20000

class DTree 
{
public:
  DTree() 
  {
    RunNo = -99;
    EvtNo = -99;
    LumiNo = -99;
    Dsize = -99;
    PVx = -99;
    PVy = -99;
    PVz = -99;
    PVnchi2 = -99;
    BSx = -99;
    BSy = -99;
    BSz = -99;
    PVxE = -99;
    PVyE = -99;
    PVzE = -99;
    BSxErr = -99;
    BSyErr = -99;
    BSzErr = -99;
    BSdxdz = -99;
    BSdydz = -99;
    BSdxdzErr = -99;
    BSdydzErr = -99;
    BSWidthX = -99;
    BSWidthXErr = -99;
    BSWidthY = -99;
    BSWidthYErr = -99;
  }
  ~DTree() {};

  DTree(TTree* t) : DTree() 
  {
    this->read_tree(t);
  }

  void read_tree(TTree* t);

  Int_t           RunNo;
  Int_t           EvtNo;
  Int_t           LumiNo;
  Int_t           Dsize;
  Float_t         PVx;
  Float_t         PVy;
  Float_t         PVz;
  Float_t         PVnchi2;
  Float_t         BSx;
  Float_t         BSy;
  Float_t         BSz;
  Float_t         PVxE;
  Float_t         PVyE;
  Float_t         PVzE;
  Float_t         BSxErr;
  Float_t         BSyErr;
  Float_t         BSzErr;
  Float_t         BSdxdz;
  Float_t         BSdydz;
  Float_t         BSdxdzErr;
  Float_t         BSdydzErr;
  Float_t         BSWidthX;
  Float_t         BSWidthXErr;
  Float_t         BSWidthY;
  Float_t         BSWidthYErr;
  Int_t           Dindex[MAX_XB];                    // [Dsize]
  Int_t           Dtype[MAX_XB];                     // [Dsize]
  Float_t         Dmass[MAX_XB];                     // [Dsize]
  Float_t         Dpt[MAX_XB];                       // [Dsize]
  Float_t         Deta[MAX_XB];                      // [Dsize]
  Float_t         Dphi[MAX_XB];                      // [Dsize]
  Float_t         Dy[MAX_XB];                        // [Dsize]
  Float_t         DvtxX[MAX_XB];                     // [Dsize]
  Float_t         DvtxY[MAX_XB];                     // [Dsize]
  Float_t         Dd0[MAX_XB];                       // [Dsize]
  Float_t         Dd0Err[MAX_XB];                    // [Dsize]
  Float_t         Ddxyz[MAX_XB];                     // [Dsize]
  Float_t         DdxyzErr[MAX_XB];                  // [Dsize]
  Float_t         Dchi2ndf[MAX_XB];                  // [Dsize]
  Float_t         Dchi2cl[MAX_XB];                   // [Dsize]
  Float_t         Ddtheta[MAX_XB];                   // [Dsize]
  Float_t         Dlxy[MAX_XB];                      // [Dsize]
  Float_t         Dalpha[MAX_XB];                    // [Dsize]
  Float_t         DsvpvDistance[MAX_XB];             // [Dsize]
  Float_t         DsvpvDisErr[MAX_XB];               // [Dsize]
  Float_t         DsvpvDistance_2D[MAX_XB];          // [Dsize]
  Float_t         DsvpvDisErr_2D[MAX_XB];            // [Dsize]
  Float_t         DtktkRes_chi2ndf[MAX_XB];          // [Dsize]
  Float_t         DtktkRes_chi2cl[MAX_XB];           // [Dsize]
  Float_t         DtktkRes_alpha[MAX_XB];            // [Dsize]
  Float_t         DtktkRes_svpvDistance[MAX_XB];     // [Dsize]
  Float_t         DtktkRes_svpvDisErr[MAX_XB];       // [Dsize]
  Float_t         DlxyBS[MAX_XB];                    // [Dsize]
  Float_t         DlxyBSErr[MAX_XB];                 // [Dsize]
  Float_t         DMaxDoca[MAX_XB];                  // [Dsize]
  Float_t         Dtrk1Pt[MAX_XB];                   // [Dsize]
  Float_t         Dtrk2Pt[MAX_XB];                   // [Dsize]
  Float_t         Dtrk1Eta[MAX_XB];                  // [Dsize]
  Float_t         Dtrk2Eta[MAX_XB];                  // [Dsize]
  Float_t         Dtrk1Phi[MAX_XB];                  // [Dsize]
  Float_t         Dtrk2Phi[MAX_XB];                  // [Dsize]
  Float_t         Dtrk1PtErr[MAX_XB];                // [Dsize]
  Float_t         Dtrk2PtErr[MAX_XB];                // [Dsize]
  Float_t         Dtrk1Dxy[MAX_XB];                  // [Dsize]
  Float_t         Dtrk2Dxy[MAX_XB];                  // [Dsize]
  Float_t         Dtrk1PixelHit[MAX_XB];             // [Dsize]
  Float_t         Dtrk2PixelHit[MAX_XB];             // [Dsize]
  Float_t         Dtrk1StripHit[MAX_XB];             // [Dsize]
  Float_t         Dtrk2StripHit[MAX_XB];             // [Dsize]
  Float_t         Dtrk1nStripLayer[MAX_XB];          // [Dsize]
  Float_t         Dtrk2nStripLayer[MAX_XB];          // [Dsize]
  Float_t         Dtrk1nPixelLayer[MAX_XB];          // [Dsize]
  Float_t         Dtrk2nPixelLayer[MAX_XB];          // [Dsize]
  Float_t         Dtrk1Chi2ndf[MAX_XB];              // [Dsize]
  Float_t         Dtrk2Chi2ndf[MAX_XB];              // [Dsize]
  Float_t         Dtrk1MassHypo[MAX_XB];             // [Dsize]
  Float_t         Dtrk2MassHypo[MAX_XB];             // [Dsize]
  Int_t           Dtrk1Algo[MAX_XB];                 // [Dsize]
  Int_t           Dtrk2Algo[MAX_XB];                 // [Dsize]
  Int_t           Dtrk1originalAlgo[MAX_XB];         // [Dsize]
  Int_t           Dtrk2originalAlgo[MAX_XB];         // [Dsize]
  Bool_t          Dtrk1highPurity[MAX_XB];           // [Dsize]
  Bool_t          Dtrk2highPurity[MAX_XB];           // [Dsize]
  Float_t         Dtrk3Pt[MAX_XB];                   // [Dsize]
  Float_t         Dtrk4Pt[MAX_XB];                   // [Dsize]
  Float_t         Dtrk3Eta[MAX_XB];                  // [Dsize]
  Float_t         Dtrk4Eta[MAX_XB];                  // [Dsize]
  Float_t         Dtrk3Phi[MAX_XB];                  // [Dsize]
  Float_t         Dtrk4Phi[MAX_XB];                  // [Dsize]
  Float_t         Dtrk3PtErr[MAX_XB];                // [Dsize]
  Float_t         Dtrk4PtErr[MAX_XB];                // [Dsize]
  Float_t         Dtrk3Dxy[MAX_XB];                  // [Dsize]
  Float_t         Dtrk4Dxy[MAX_XB];                  // [Dsize]
  Float_t         Dtrk3PixelHit[MAX_XB];             // [Dsize]
  Float_t         Dtrk4PixelHit[MAX_XB];             // [Dsize]
  Float_t         Dtrk3StripHit[MAX_XB];             // [Dsize]
  Float_t         Dtrk4StripHit[MAX_XB];             // [Dsize]
  Float_t         Dtrk3nStripLayer[MAX_XB];          // [Dsize]
  Float_t         Dtrk4nStripLayer[MAX_XB];          // [Dsize]
  Float_t         Dtrk3nPixelLayer[MAX_XB];          // [Dsize]
  Float_t         Dtrk4nPixelLayer[MAX_XB];          // [Dsize]
  Float_t         Dtrk3Chi2ndf[MAX_XB];              // [Dsize]
  Float_t         Dtrk4Chi2ndf[MAX_XB];              // [Dsize]
  Float_t         Dtrk3MassHypo[MAX_XB];             // [Dsize]
  Float_t         Dtrk4MassHypo[MAX_XB];             // [Dsize]
  Int_t           Dtrk3Algo[MAX_XB];                 // [Dsize]
  Int_t           Dtrk4Algo[MAX_XB];                 // [Dsize]
  Int_t           Dtrk3originalAlgo[MAX_XB];         // [Dsize]
  Int_t           Dtrk4originalAlgo[MAX_XB];         // [Dsize]
  Bool_t          Dtrk3highPurity[MAX_XB];           // [Dsize]
  Bool_t          Dtrk4highPurity[MAX_XB];           // [Dsize]
  Int_t           Dtrk1Idx[MAX_XB];                  // [Dsize]
  Int_t           Dtrk2Idx[MAX_XB];                  // [Dsize]
  Float_t         Dtrk1EtaErr[MAX_XB];               // [Dsize]
  Float_t         Dtrk2EtaErr[MAX_XB];               // [Dsize]
  Float_t         Dtrk1PhiErr[MAX_XB];               // [Dsize]
  Float_t         Dtrk2PhiErr[MAX_XB];               // [Dsize]
  Float_t         Dtrk1Y[MAX_XB];                    // [Dsize]
  Float_t         Dtrk2Y[MAX_XB];                    // [Dsize]
  Float_t         Dtrk1D0Err[MAX_XB];                // [Dsize]
  Float_t         Dtrk2D0Err[MAX_XB];                // [Dsize]
  Float_t         Dtrk1MVAVal[MAX_XB];               // [Dsize]
  Float_t         Dtrk2MVAVal[MAX_XB];               // [Dsize]
  Int_t           Dtrk1Quality[MAX_XB];              // [Dsize]
  Int_t           Dtrk2Quality[MAX_XB];              // [Dsize]
  Int_t           Dtrk3Idx[MAX_XB];                  // [Dsize]
  Int_t           Dtrk4Idx[MAX_XB];                  // [Dsize]
  Float_t         Dtrk3EtaErr[MAX_XB];               // [Dsize]
  Float_t         Dtrk4EtaErr[MAX_XB];               // [Dsize]
  Float_t         Dtrk3PhiErr[MAX_XB];               // [Dsize]
  Float_t         Dtrk4PhiErr[MAX_XB];               // [Dsize]
  Float_t         Dtrk3Y[MAX_XB];                    // [Dsize]
  Float_t         Dtrk4Y[MAX_XB];                    // [Dsize]
  Float_t         Dtrk3D0Err[MAX_XB];                // [Dsize]
  Float_t         Dtrk4D0Err[MAX_XB];                // [Dsize]
  Float_t         Dtrk3MVAVal[MAX_XB];               // [Dsize]
  Float_t         Dtrk4MVAVal[MAX_XB];               // [Dsize]
  Int_t           Dtrk3Quality[MAX_XB];              // [Dsize]
  Int_t           Dtrk4Quality[MAX_XB];              // [Dsize]
  Float_t         DtktkResmass[MAX_XB];              // [Dsize]
  Float_t         DtktkRespt[MAX_XB];                // [Dsize]
  Float_t         DtktkReseta[MAX_XB];               // [Dsize]
  Float_t         DtktkResphi[MAX_XB];               // [Dsize]
  Float_t         DRestrk1Pt[MAX_XB];                // [Dsize]
  Float_t         DRestrk1Eta[MAX_XB];               // [Dsize]
  Float_t         DRestrk1Phi[MAX_XB];               // [Dsize]
  Float_t         DRestrk1Y[MAX_XB];                 // [Dsize]
  Float_t         DRestrk1Dxy[MAX_XB];               // [Dsize]
  Float_t         DRestrk1D0Err[MAX_XB];             // [Dsize]
  Int_t           DRestrk1originalAlgo[MAX_XB];      // [Dsize]
  Float_t         DRestrk2Pt[MAX_XB];                // [Dsize]
  Float_t         DRestrk2Eta[MAX_XB];               // [Dsize]
  Float_t         DRestrk2Phi[MAX_XB];               // [Dsize]
  Float_t         DRestrk2Y[MAX_XB];                 // [Dsize]
  Float_t         DRestrk2Dxy[MAX_XB];               // [Dsize]
  Float_t         DRestrk2D0Err[MAX_XB];             // [Dsize]
  Int_t           DRestrk2originalAlgo[MAX_XB];      // [Dsize]
  Float_t         DRestrk3Pt[MAX_XB];                // [Dsize]
  Float_t         DRestrk3Eta[MAX_XB];               // [Dsize]
  Float_t         DRestrk3Phi[MAX_XB];               // [Dsize]
  Float_t         DRestrk3Y[MAX_XB];                 // [Dsize]
  Float_t         DRestrk3Dxy[MAX_XB];               // [Dsize]
  Float_t         DRestrk3D0Err[MAX_XB];             // [Dsize]
  Int_t           DRestrk3originalAlgo[MAX_XB];      // [Dsize]
  Float_t         DRestrk4Pt[MAX_XB];                // [Dsize]
  Float_t         DRestrk4Eta[MAX_XB];               // [Dsize]
  Float_t         DRestrk4Phi[MAX_XB];               // [Dsize]
  Float_t         DRestrk4Y[MAX_XB];                 // [Dsize]
  Float_t         DRestrk4Dxy[MAX_XB];               // [Dsize]
  Float_t         DRestrk4D0Err[MAX_XB];             // [Dsize]
  Int_t           DRestrk4originalAlgo[MAX_XB];      // [Dsize]
  Float_t         Dgen[MAX_XB];                      // [Dsize]
  Int_t           DgenIndex[MAX_XB];                 // [Dsize]
  Int_t           DgennDa[MAX_XB];                   // [Dsize]
  Float_t         Dgenpt[MAX_XB];                    // [Dsize]
  Float_t         Dgeneta[MAX_XB];                   // [Dsize]
  Float_t         Dgenphi[MAX_XB];                   // [Dsize]
  Float_t         Dgeny[MAX_XB];                     // [Dsize]
  Int_t           DgencollisionId[MAX_XB];           // [Dsize]
  Float_t         DgenBAncestorpt[MAX_XB];           // [Dsize]
  Int_t           DgenBAncestorpdgId[MAX_XB];        // [Dsize]
};

void DTree::read_tree(TTree* t) 
{
  t->SetBranchStatus("RunNo", 1);
  t->SetBranchStatus("EvtNo", 1);
  t->SetBranchStatus("LumiNo", 1);
  t->SetBranchStatus("Dsize", 1);
  t->SetBranchStatus("PVx", 1);
  t->SetBranchStatus("PVy", 1);
  t->SetBranchStatus("PVz", 1);
  t->SetBranchStatus("PVnchi2", 1);
  t->SetBranchStatus("BSx", 1);
  t->SetBranchStatus("BSy", 1);
  t->SetBranchStatus("BSz", 1);
  t->SetBranchStatus("PVxE", 1);
  t->SetBranchStatus("PVyE", 1);
  t->SetBranchStatus("PVzE", 1);
  t->SetBranchStatus("BSxErr", 1);
  t->SetBranchStatus("BSyErr", 1);
  t->SetBranchStatus("BSzErr", 1);
  t->SetBranchStatus("BSdxdz", 1);
  t->SetBranchStatus("BSdydz", 1);
  t->SetBranchStatus("BSdxdzErr", 1);
  t->SetBranchStatus("BSdydzErr", 1);
  t->SetBranchStatus("BSWidthX", 1);
  t->SetBranchStatus("BSWidthXErr", 1);
  t->SetBranchStatus("BSWidthY", 1);
  t->SetBranchStatus("BSWidthYErr", 1);

  t->SetBranchStatus("Dindex", 1);
  t->SetBranchStatus("Dtype", 1);
  t->SetBranchStatus("Dmass", 1);
  t->SetBranchStatus("Dpt", 1);
  t->SetBranchStatus("Deta", 1);
  t->SetBranchStatus("Dphi", 1);
  t->SetBranchStatus("Dy", 1);
  t->SetBranchStatus("DvtxX", 1);
  t->SetBranchStatus("DvtxY", 1);
  t->SetBranchStatus("Dd0", 1);
  t->SetBranchStatus("Dd0Err", 1);
  t->SetBranchStatus("Ddxyz", 1);
  t->SetBranchStatus("DdxyzErr", 1);
  t->SetBranchStatus("Dchi2ndf", 1);
  t->SetBranchStatus("Dchi2cl", 1);
  t->SetBranchStatus("Ddtheta", 1);
  t->SetBranchStatus("Dlxy", 1);
  t->SetBranchStatus("Dalpha", 1);
  t->SetBranchStatus("DsvpvDistance", 1);
  t->SetBranchStatus("DsvpvDisErr", 1);
  t->SetBranchStatus("DsvpvDistance_2D", 1);
  t->SetBranchStatus("DsvpvDisErr_2D", 1);
  t->SetBranchStatus("DtktkRes_chi2ndf", 1);
  t->SetBranchStatus("DtktkRes_chi2cl", 1);
  t->SetBranchStatus("DtktkRes_alpha", 1);
  t->SetBranchStatus("DtktkRes_svpvDistance", 1);
  t->SetBranchStatus("DtktkRes_svpvDisErr", 1);
  t->SetBranchStatus("DlxyBS", 1);
  t->SetBranchStatus("DlxyBSErr", 1);
  t->SetBranchStatus("DMaxDoca", 1);
  t->SetBranchStatus("Dtrk1Pt", 1);
  t->SetBranchStatus("Dtrk2Pt", 1);
  t->SetBranchStatus("Dtrk1Eta", 1);
  t->SetBranchStatus("Dtrk2Eta", 1);
  t->SetBranchStatus("Dtrk1Phi", 1);
  t->SetBranchStatus("Dtrk2Phi", 1);
  t->SetBranchStatus("Dtrk1PtErr", 1);
  t->SetBranchStatus("Dtrk2PtErr", 1);
  t->SetBranchStatus("Dtrk1Dxy", 1);
  t->SetBranchStatus("Dtrk2Dxy", 1);
  t->SetBranchStatus("Dtrk1PixelHit", 1);
  t->SetBranchStatus("Dtrk2PixelHit", 1);
  t->SetBranchStatus("Dtrk1StripHit", 1);
  t->SetBranchStatus("Dtrk2StripHit", 1);
  t->SetBranchStatus("Dtrk1nStripLayer", 1);
  t->SetBranchStatus("Dtrk2nStripLayer", 1);
  t->SetBranchStatus("Dtrk1nPixelLayer", 1);
  t->SetBranchStatus("Dtrk2nPixelLayer", 1);
  t->SetBranchStatus("Dtrk1Chi2ndf", 1);
  t->SetBranchStatus("Dtrk2Chi2ndf", 1);
  t->SetBranchStatus("Dtrk1MassHypo", 1);
  t->SetBranchStatus("Dtrk2MassHypo", 1);
  t->SetBranchStatus("Dtrk1Algo", 1);
  t->SetBranchStatus("Dtrk2Algo", 1);
  t->SetBranchStatus("Dtrk1originalAlgo", 1);
  t->SetBranchStatus("Dtrk2originalAlgo", 1);
  t->SetBranchStatus("Dtrk1highPurity", 1);
  t->SetBranchStatus("Dtrk2highPurity", 1);
  t->SetBranchStatus("Dtrk3Pt", 1);
  t->SetBranchStatus("Dtrk4Pt", 1);
  t->SetBranchStatus("Dtrk3Eta", 1);
  t->SetBranchStatus("Dtrk4Eta", 1);
  t->SetBranchStatus("Dtrk3Phi", 1);
  t->SetBranchStatus("Dtrk4Phi", 1);
  t->SetBranchStatus("Dtrk3PtErr", 1);
  t->SetBranchStatus("Dtrk4PtErr", 1);
  t->SetBranchStatus("Dtrk3Dxy", 1);
  t->SetBranchStatus("Dtrk4Dxy", 1);
  t->SetBranchStatus("Dtrk3PixelHit", 1);
  t->SetBranchStatus("Dtrk4PixelHit", 1);
  t->SetBranchStatus("Dtrk3StripHit", 1);
  t->SetBranchStatus("Dtrk4StripHit", 1);
  t->SetBranchStatus("Dtrk3nStripLayer", 1);
  t->SetBranchStatus("Dtrk4nStripLayer", 1);
  t->SetBranchStatus("Dtrk3nPixelLayer", 1);
  t->SetBranchStatus("Dtrk4nPixelLayer", 1);
  t->SetBranchStatus("Dtrk3Chi2ndf", 1);
  t->SetBranchStatus("Dtrk4Chi2ndf", 1);
  t->SetBranchStatus("Dtrk3MassHypo", 1);
  t->SetBranchStatus("Dtrk4MassHypo", 1);
  t->SetBranchStatus("Dtrk3Algo", 1);
  t->SetBranchStatus("Dtrk4Algo", 1);
  t->SetBranchStatus("Dtrk3originalAlgo", 1);
  t->SetBranchStatus("Dtrk4originalAlgo", 1);
  t->SetBranchStatus("Dtrk3highPurity", 1);
  t->SetBranchStatus("Dtrk4highPurity", 1);
  t->SetBranchStatus("Dtrk1Idx", 1);
  t->SetBranchStatus("Dtrk2Idx", 1);
  t->SetBranchStatus("Dtrk1EtaErr", 1);
  t->SetBranchStatus("Dtrk2EtaErr", 1);
  t->SetBranchStatus("Dtrk1PhiErr", 1);
  t->SetBranchStatus("Dtrk2PhiErr", 1);
  t->SetBranchStatus("Dtrk1Y", 1);
  t->SetBranchStatus("Dtrk2Y", 1);
  t->SetBranchStatus("Dtrk1D0Err", 1);
  t->SetBranchStatus("Dtrk2D0Err", 1);
  t->SetBranchStatus("Dtrk1MVAVal", 1);
  t->SetBranchStatus("Dtrk2MVAVal", 1);
  t->SetBranchStatus("Dtrk1Quality", 1);
  t->SetBranchStatus("Dtrk2Quality", 1);
  t->SetBranchStatus("Dtrk3Idx", 1);
  t->SetBranchStatus("Dtrk4Idx", 1);
  t->SetBranchStatus("Dtrk3EtaErr", 1);
  t->SetBranchStatus("Dtrk4EtaErr", 1);
  t->SetBranchStatus("Dtrk3PhiErr", 1);
  t->SetBranchStatus("Dtrk4PhiErr", 1);
  t->SetBranchStatus("Dtrk3Y", 1);
  t->SetBranchStatus("Dtrk4Y", 1);
  t->SetBranchStatus("Dtrk3D0Err", 1);
  t->SetBranchStatus("Dtrk4D0Err", 1);
  t->SetBranchStatus("Dtrk3MVAVal", 1);
  t->SetBranchStatus("Dtrk4MVAVal", 1);
  t->SetBranchStatus("Dtrk3Quality", 1);
  t->SetBranchStatus("Dtrk4Quality", 1);
  t->SetBranchStatus("DtktkResmass", 1);
  t->SetBranchStatus("DtktkRespt", 1);
  t->SetBranchStatus("DtktkReseta", 1);
  t->SetBranchStatus("DtktkResphi", 1);
  t->SetBranchStatus("DRestrk1Pt", 1);
  t->SetBranchStatus("DRestrk1Eta", 1);
  t->SetBranchStatus("DRestrk1Phi", 1);
  t->SetBranchStatus("DRestrk1Y", 1);
  t->SetBranchStatus("DRestrk1Dxy", 1);
  t->SetBranchStatus("DRestrk1D0Err", 1);
  t->SetBranchStatus("DRestrk1originalAlgo", 1);
  t->SetBranchStatus("DRestrk2Pt", 1);
  t->SetBranchStatus("DRestrk2Eta", 1);
  t->SetBranchStatus("DRestrk2Phi", 1);
  t->SetBranchStatus("DRestrk2Y", 1);
  t->SetBranchStatus("DRestrk2Dxy", 1);
  t->SetBranchStatus("DRestrk2D0Err", 1);
  t->SetBranchStatus("DRestrk2originalAlgo", 1);
  t->SetBranchStatus("DRestrk3Pt", 1);
  t->SetBranchStatus("DRestrk3Eta", 1);
  t->SetBranchStatus("DRestrk3Phi", 1);
  t->SetBranchStatus("DRestrk3Y", 1);
  t->SetBranchStatus("DRestrk3Dxy", 1);
  t->SetBranchStatus("DRestrk3D0Err", 1);
  t->SetBranchStatus("DRestrk3originalAlgo", 1);
  t->SetBranchStatus("DRestrk4Pt", 1);
  t->SetBranchStatus("DRestrk4Eta", 1);
  t->SetBranchStatus("DRestrk4Phi", 1);
  t->SetBranchStatus("DRestrk4Y", 1);
  t->SetBranchStatus("DRestrk4Dxy", 1);
  t->SetBranchStatus("DRestrk4D0Err", 1);
  t->SetBranchStatus("DRestrk4originalAlgo", 1);
  t->SetBranchStatus("Dgen", 1);
  t->SetBranchStatus("DgenIndex", 1);
  t->SetBranchStatus("DgennDa", 1);
  t->SetBranchStatus("Dgenpt", 1);
  t->SetBranchStatus("Dgeneta", 1);
  t->SetBranchStatus("Dgenphi", 1);
  t->SetBranchStatus("Dgeny", 1);
  t->SetBranchStatus("DgencollisionId", 1);
  t->SetBranchStatus("DgenBAncestorpt", 1);
  t->SetBranchStatus("DgenBAncestorpdgId", 1);

  t->SetBranchAddress("RunNo", &RunNo);
  t->SetBranchAddress("EvtNo", &EvtNo);
  t->SetBranchAddress("LumiNo", &LumiNo);
  t->SetBranchAddress("Dsize", &Dsize);
  t->SetBranchAddress("PVx", &PVx);
  t->SetBranchAddress("PVy", &PVy);
  t->SetBranchAddress("PVz", &PVz);
  t->SetBranchAddress("PVnchi2", &PVnchi2);
  t->SetBranchAddress("BSx", &BSx);
  t->SetBranchAddress("BSy", &BSy);
  t->SetBranchAddress("BSz", &BSz);
  t->SetBranchAddress("PVxE", &PVxE);
  t->SetBranchAddress("PVyE", &PVyE);
  t->SetBranchAddress("PVzE", &PVzE);
  t->SetBranchAddress("BSxErr", &BSxErr);
  t->SetBranchAddress("BSyErr", &BSyErr);
  t->SetBranchAddress("BSzErr", &BSzErr);
  t->SetBranchAddress("BSdxdz", &BSdxdz);
  t->SetBranchAddress("BSdydz", &BSdydz);
  t->SetBranchAddress("BSdxdzErr", &BSdxdzErr);
  t->SetBranchAddress("BSdydzErr", &BSdydzErr);
  t->SetBranchAddress("BSWidthX", &BSWidthX);
  t->SetBranchAddress("BSWidthXErr", &BSWidthXErr);
  t->SetBranchAddress("BSWidthY", &BSWidthY);
  t->SetBranchAddress("BSWidthYErr", &BSWidthYErr);

  t->SetBranchAddress("Dindex", Dindex);
  t->SetBranchAddress("Dtype", Dtype);
  t->SetBranchAddress("Dmass", Dmass);
  t->SetBranchAddress("Dpt", Dpt);
  t->SetBranchAddress("Deta", Deta);
  t->SetBranchAddress("Dphi", Dphi);
  t->SetBranchAddress("Dy", Dy);
  t->SetBranchAddress("DvtxX", DvtxX);
  t->SetBranchAddress("DvtxY", DvtxY);
  t->SetBranchAddress("Dd0", Dd0);
  t->SetBranchAddress("Dd0Err", Dd0Err);
  t->SetBranchAddress("Ddxyz", Ddxyz);
  t->SetBranchAddress("DdxyzErr", DdxyzErr);
  t->SetBranchAddress("Dchi2ndf", Dchi2ndf);
  t->SetBranchAddress("Dchi2cl", Dchi2cl);
  t->SetBranchAddress("Ddtheta", Ddtheta);
  t->SetBranchAddress("Dlxy", Dlxy);
  t->SetBranchAddress("Dalpha", Dalpha);
  t->SetBranchAddress("DsvpvDistance", DsvpvDistance);
  t->SetBranchAddress("DsvpvDisErr", DsvpvDisErr);
  t->SetBranchAddress("DsvpvDistance_2D", DsvpvDistance_2D);
  t->SetBranchAddress("DsvpvDisErr_2D", DsvpvDisErr_2D);
  t->SetBranchAddress("DtktkRes_chi2ndf", DtktkRes_chi2ndf);
  t->SetBranchAddress("DtktkRes_chi2cl", DtktkRes_chi2cl);
  t->SetBranchAddress("DtktkRes_alpha", DtktkRes_alpha);
  t->SetBranchAddress("DtktkRes_svpvDistance", DtktkRes_svpvDistance);
  t->SetBranchAddress("DtktkRes_svpvDisErr", DtktkRes_svpvDisErr);
  t->SetBranchAddress("DlxyBS", DlxyBS);
  t->SetBranchAddress("DlxyBSErr", DlxyBSErr);
  t->SetBranchAddress("DMaxDoca", DMaxDoca);
  t->SetBranchAddress("Dtrk1Pt", Dtrk1Pt);
  t->SetBranchAddress("Dtrk2Pt", Dtrk2Pt);
  t->SetBranchAddress("Dtrk1Eta", Dtrk1Eta);
  t->SetBranchAddress("Dtrk2Eta", Dtrk2Eta);
  t->SetBranchAddress("Dtrk1Phi", Dtrk1Phi);
  t->SetBranchAddress("Dtrk2Phi", Dtrk2Phi);
  t->SetBranchAddress("Dtrk1PtErr", Dtrk1PtErr);
  t->SetBranchAddress("Dtrk2PtErr", Dtrk2PtErr);
  t->SetBranchAddress("Dtrk1Dxy", Dtrk1Dxy);
  t->SetBranchAddress("Dtrk2Dxy", Dtrk2Dxy);
  t->SetBranchAddress("Dtrk1PixelHit", Dtrk1PixelHit);
  t->SetBranchAddress("Dtrk2PixelHit", Dtrk2PixelHit);
  t->SetBranchAddress("Dtrk1StripHit", Dtrk1StripHit);
  t->SetBranchAddress("Dtrk2StripHit", Dtrk2StripHit);
  t->SetBranchAddress("Dtrk1nStripLayer", Dtrk1nStripLayer);
  t->SetBranchAddress("Dtrk2nStripLayer", Dtrk2nStripLayer);
  t->SetBranchAddress("Dtrk1nPixelLayer", Dtrk1nPixelLayer);
  t->SetBranchAddress("Dtrk2nPixelLayer", Dtrk2nPixelLayer);
  t->SetBranchAddress("Dtrk1Chi2ndf", Dtrk1Chi2ndf);
  t->SetBranchAddress("Dtrk2Chi2ndf", Dtrk2Chi2ndf);
  t->SetBranchAddress("Dtrk1MassHypo", Dtrk1MassHypo);
  t->SetBranchAddress("Dtrk2MassHypo", Dtrk2MassHypo);
  t->SetBranchAddress("Dtrk1Algo", Dtrk1Algo);
  t->SetBranchAddress("Dtrk2Algo", Dtrk2Algo);
  t->SetBranchAddress("Dtrk1originalAlgo", Dtrk1originalAlgo);
  t->SetBranchAddress("Dtrk2originalAlgo", Dtrk2originalAlgo);
  t->SetBranchAddress("Dtrk1highPurity", Dtrk1highPurity);
  t->SetBranchAddress("Dtrk2highPurity", Dtrk2highPurity);
  t->SetBranchAddress("Dtrk3Pt", Dtrk3Pt);
  t->SetBranchAddress("Dtrk4Pt", Dtrk4Pt);
  t->SetBranchAddress("Dtrk3Eta", Dtrk3Eta);
  t->SetBranchAddress("Dtrk4Eta", Dtrk4Eta);
  t->SetBranchAddress("Dtrk3Phi", Dtrk3Phi);
  t->SetBranchAddress("Dtrk4Phi", Dtrk4Phi);
  t->SetBranchAddress("Dtrk3PtErr", Dtrk3PtErr);
  t->SetBranchAddress("Dtrk4PtErr", Dtrk4PtErr);
  t->SetBranchAddress("Dtrk3Dxy", Dtrk3Dxy);
  t->SetBranchAddress("Dtrk4Dxy", Dtrk4Dxy);
  t->SetBranchAddress("Dtrk3PixelHit", Dtrk3PixelHit);
  t->SetBranchAddress("Dtrk4PixelHit", Dtrk4PixelHit);
  t->SetBranchAddress("Dtrk3StripHit", Dtrk3StripHit);
  t->SetBranchAddress("Dtrk4StripHit", Dtrk4StripHit);
  t->SetBranchAddress("Dtrk3nStripLayer", Dtrk3nStripLayer);
  t->SetBranchAddress("Dtrk4nStripLayer", Dtrk4nStripLayer);
  t->SetBranchAddress("Dtrk3nPixelLayer", Dtrk3nPixelLayer);
  t->SetBranchAddress("Dtrk4nPixelLayer", Dtrk4nPixelLayer);
  t->SetBranchAddress("Dtrk3Chi2ndf", Dtrk3Chi2ndf);
  t->SetBranchAddress("Dtrk4Chi2ndf", Dtrk4Chi2ndf);
  t->SetBranchAddress("Dtrk3MassHypo", Dtrk3MassHypo);
  t->SetBranchAddress("Dtrk4MassHypo", Dtrk4MassHypo);
  t->SetBranchAddress("Dtrk3Algo", Dtrk3Algo);
  t->SetBranchAddress("Dtrk4Algo", Dtrk4Algo);
  t->SetBranchAddress("Dtrk3originalAlgo", Dtrk3originalAlgo);
  t->SetBranchAddress("Dtrk4originalAlgo", Dtrk4originalAlgo);
  t->SetBranchAddress("Dtrk3highPurity", Dtrk3highPurity);
  t->SetBranchAddress("Dtrk4highPurity", Dtrk4highPurity);
  t->SetBranchAddress("Dtrk1Idx", Dtrk1Idx);
  t->SetBranchAddress("Dtrk2Idx", Dtrk2Idx);
  t->SetBranchAddress("Dtrk1EtaErr", Dtrk1EtaErr);
  t->SetBranchAddress("Dtrk2EtaErr", Dtrk2EtaErr);
  t->SetBranchAddress("Dtrk1PhiErr", Dtrk1PhiErr);
  t->SetBranchAddress("Dtrk2PhiErr", Dtrk2PhiErr);
  t->SetBranchAddress("Dtrk1Y", Dtrk1Y);
  t->SetBranchAddress("Dtrk2Y", Dtrk2Y);
  t->SetBranchAddress("Dtrk1D0Err", Dtrk1D0Err);
  t->SetBranchAddress("Dtrk2D0Err", Dtrk2D0Err);
  t->SetBranchAddress("Dtrk1MVAVal", Dtrk1MVAVal);
  t->SetBranchAddress("Dtrk2MVAVal", Dtrk2MVAVal);
  t->SetBranchAddress("Dtrk1Quality", Dtrk1Quality);
  t->SetBranchAddress("Dtrk2Quality", Dtrk2Quality);
  t->SetBranchAddress("Dtrk3Idx", Dtrk3Idx);
  t->SetBranchAddress("Dtrk4Idx", Dtrk4Idx);
  t->SetBranchAddress("Dtrk3EtaErr", Dtrk3EtaErr);
  t->SetBranchAddress("Dtrk4EtaErr", Dtrk4EtaErr);
  t->SetBranchAddress("Dtrk3PhiErr", Dtrk3PhiErr);
  t->SetBranchAddress("Dtrk4PhiErr", Dtrk4PhiErr);
  t->SetBranchAddress("Dtrk3Y", Dtrk3Y);
  t->SetBranchAddress("Dtrk4Y", Dtrk4Y);
  t->SetBranchAddress("Dtrk3D0Err", Dtrk3D0Err);
  t->SetBranchAddress("Dtrk4D0Err", Dtrk4D0Err);
  t->SetBranchAddress("Dtrk3MVAVal", Dtrk3MVAVal);
  t->SetBranchAddress("Dtrk4MVAVal", Dtrk4MVAVal);
  t->SetBranchAddress("Dtrk3Quality", Dtrk3Quality);
  t->SetBranchAddress("Dtrk4Quality", Dtrk4Quality);
  t->SetBranchAddress("DtktkResmass", DtktkResmass);
  t->SetBranchAddress("DtktkRespt", DtktkRespt);
  t->SetBranchAddress("DtktkReseta", DtktkReseta);
  t->SetBranchAddress("DtktkResphi", DtktkResphi);
  t->SetBranchAddress("DRestrk1Pt", DRestrk1Pt);
  t->SetBranchAddress("DRestrk1Eta", DRestrk1Eta);
  t->SetBranchAddress("DRestrk1Phi", DRestrk1Phi);
  t->SetBranchAddress("DRestrk1Y", DRestrk1Y);
  t->SetBranchAddress("DRestrk1Dxy", DRestrk1Dxy);
  t->SetBranchAddress("DRestrk1D0Err", DRestrk1D0Err);
  t->SetBranchAddress("DRestrk1originalAlgo", DRestrk1originalAlgo);
  t->SetBranchAddress("DRestrk2Pt", DRestrk2Pt);
  t->SetBranchAddress("DRestrk2Eta", DRestrk2Eta);
  t->SetBranchAddress("DRestrk2Phi", DRestrk2Phi);
  t->SetBranchAddress("DRestrk2Y", DRestrk2Y);
  t->SetBranchAddress("DRestrk2Dxy", DRestrk2Dxy);
  t->SetBranchAddress("DRestrk2D0Err", DRestrk2D0Err);
  t->SetBranchAddress("DRestrk2originalAlgo", DRestrk2originalAlgo);
  t->SetBranchAddress("DRestrk3Pt", DRestrk3Pt);
  t->SetBranchAddress("DRestrk3Eta", DRestrk3Eta);
  t->SetBranchAddress("DRestrk3Phi", DRestrk3Phi);
  t->SetBranchAddress("DRestrk3Y", DRestrk3Y);
  t->SetBranchAddress("DRestrk3Dxy", DRestrk3Dxy);
  t->SetBranchAddress("DRestrk3D0Err", DRestrk3D0Err);
  t->SetBranchAddress("DRestrk3originalAlgo", DRestrk3originalAlgo);
  t->SetBranchAddress("DRestrk4Pt", DRestrk4Pt);
  t->SetBranchAddress("DRestrk4Eta", DRestrk4Eta);
  t->SetBranchAddress("DRestrk4Phi", DRestrk4Phi);
  t->SetBranchAddress("DRestrk4Y", DRestrk4Y);
  t->SetBranchAddress("DRestrk4Dxy", DRestrk4Dxy);
  t->SetBranchAddress("DRestrk4D0Err", DRestrk4D0Err);
  t->SetBranchAddress("DRestrk4originalAlgo", DRestrk4originalAlgo);
  t->SetBranchAddress("Dgen", Dgen);
  t->SetBranchAddress("DgenIndex", DgenIndex);
  t->SetBranchAddress("DgennDa", DgennDa);
  t->SetBranchAddress("Dgenpt", Dgenpt);
  t->SetBranchAddress("Dgeneta", Dgeneta);
  t->SetBranchAddress("Dgenphi", Dgenphi);
  t->SetBranchAddress("Dgeny", Dgeny);
  t->SetBranchAddress("DgencollisionId", DgencollisionId);
  t->SetBranchAddress("DgenBAncestorpt", DgenBAncestorpt);
  t->SetBranchAddress("DgenBAncestorpdgId", DgenBAncestorpdgId);
}

#endif
