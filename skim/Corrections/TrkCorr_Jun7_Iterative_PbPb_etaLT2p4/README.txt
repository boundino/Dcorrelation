To use:

place all files in this folder into a subdirectory
include getTrkCorr.h header file (which should also grab the other header file)
make a trkCorr object at the start of your program somewhere:

TrkCorr* trkCorr = new TrkCorr("*Name of subdirectory*/"); //must have the '/' after the subdirectory name

To call get the correction, which is applied multiplicatively, use:

trkCorr->getTrkCorr(pt,eta,phi,hiBin,rmin);


********************************************************************
rmin is calculated thusly (assumes jet array is ordered by pt):
uses akPu4Calo jets
//find rmin parameters for the track
      float rmin = 999;
      for(int k = 0; k<nref; k++)
      {
        if(jtpt[k]<50) break;
        if((TMath::Abs(chargedSum[k]/rawpt[k])<0.01) || (TMath::Abs(jteta[k]>2))) continue;//jet quality cut
        float R = TMath::Power(jteta[k]-trkEta[j],2)+TMath::Power(TMath::ACos(TMath::Cos(jtphi[k]-trkPhi[j]))),2);
        if(rmin*rmin>R) rmin=TMath::Power(R,0.5);
      }



******************************************************************
Tracking selection to be used:

|eta|<2.4, 0.5<pt<400 GeV (0.5-0.7 has an experimental cut currently)
if(highPurity[j]!=1) continue;
if(trkPtError[j]/trkPt[j]>0.1 || TMath::Abs(trkDz1[j]/trkDzError1[j])>3 || TMath::Abs(trkDxy1[j]/trkDxyError1[j])>3) continue;
if(trkChi2[j]/(float)trkNdof[j]/(float)trkNlayer[j]>0.15) continue;
if(trkNHit[j]<11 && trkPt[j]>0.7) continue;

float Et = (pfHcal[j]+pfEcal[j])/TMath::CosH(trkEta[j]);
if(!(trkPt[j]<20 || (Et>0.5*trkPt[j]))) continue;


