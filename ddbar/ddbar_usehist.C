#include "ddbar.h"

void ddbar_usehist(TString infhistname, TString outfname,
                   TString collisionsyst)
// Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
{
  int arguerr(TString collisionsyst);
  if(arguerr(collisionsyst)) return;

  if(initbinning()) return;

  TFile* infD = new TFile(Form("%s.root",infhistname.Data()));
  if(gethists(std::vector<TFile*>{infD}, "usehist")) return;

  //
  ahdphi_plot.push_back(ahdphiD[0]);
  ahdphi_plot.push_back(ahdphiD[1]);
  ahdphi_plot.push_back(ahdphiG[0]);
  ahdphi_plot.push_back(ahdphiG[1]);

  //
  TFile* outf = new TFile(Form("%s.root",outfname.Data()),"recreate");
  outf->cd();
  if(writehists("usehist")) return;
  outf->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==4)
    {
      ddbar_usehist(argv[1], argv[2], argv[3]);
      return 0;
    }
  else
    {
      return 1;
    }
}

int arguerr(TString collisionsyst)
{
  if(collsyst_list.find(collisionsyst)==collsyst_list.end())
    {
      std::cout<<"\033[1;31merror:\033[0m invalid \"collisionsyst\""<<std::endl;
      return 1;
    }
  return 0;
}
