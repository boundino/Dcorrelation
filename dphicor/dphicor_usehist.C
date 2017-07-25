#include "dphicor.h"

int dphicor_usehist(TString outfname, TString outplotname, TString collisionsyst, Int_t isMC, Float_t leading_ptmin, Float_t other_ptmin)
{
  std::cout<<std::endl;
  xjjrootuti::setgstyle();

  TFile* inf = new TFile(Form("%s.root",outfname.Data()));
  TH1D* hdphi = (TH1D*)inf->Get("hdphi");
  hdphi->Sumw2();
  xjjrootuti::setth1(hdphi);

  TCanvas* cdphi = new TCanvas("cdphi","",600,600);  
  cdphi->SetLogy();
  hdphi->Draw("pe");
  xjjrootuti::drawCMS(collisionsyst);
  cdphi->SaveAs(Form("plots/cdphi_%s.pdf",outplotname.Data()));

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc==7)
    {
      dphicor_usehist(argv[1], argv[2], argv[3], atoi(argv[4]), atof(argv[5]), atof(argv[6]));
      return 0;
    }
  else
    {
      return 1;
    }
}
