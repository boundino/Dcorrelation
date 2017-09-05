#include "ptasymm.h"

void ptasymm_plothist(TString infname, TString outplotname, 
                      TString collisionsyst, Int_t isMC, 
                      Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
{
  int arguerr(TString collisionsyst, Int_t isMC);
  if(arguerr(collisionsyst, isMC)) return;

  xjjroot::setgstyle();
  TFile* inf = new TFile(Form("%s.root",infname.Data()));
  if(gethists(std::vector<TFile*>{inf}, "plothist")) return;

  std::vector<TString> canvdraw = {"data", "base", "bkgsub", "fitext", "final"};
  std::vector<std::vector<TString>> drawhist = 
    {
      std::vector<TString>{"hptasym_all_all_hist", "hptasym_all_all",         "hptasym_all_all_fit",        "hptasym_subtract_all_fit"},
      std::vector<TString>{"hptasym_all_all_hist", "hptasym_all_signal_hist", "hptasym_signal_signal_hist", "hptasym_all_all", "hptasym_all_signal", "hptasym_signal_signal"},
      std::vector<TString>{"hptasym_all_all_hist", "hptasym_all_signal_hist", "hptasym_signal_signal_hist", "hptasym_all_all", "hptasym_all_signal", "hptasym_signal_signal", "hptasym_subtract_signal"},
      std::vector<TString>{"hptasym_all_all_hist", "hptasym_all_signal_hist", "hptasym_signal_signal_hist", "hptasym_all_all", "hptasym_all_signal", "hptasym_signal_signal", "hptasym_all_all_fit"},
      std::vector<TString>{"hptasym_all_all_hist", "hptasym_all_signal_hist", "hptasym_signal_signal_hist", "hptasym_all_all", "hptasym_all_signal", "hptasym_signal_signal", "hptasym_subtract_all_fit"}
    };
  if(canvdraw.size()!=drawhist.size()) return;
  
  //  
  TH2F* hempty = new TH2F("hempty", ";#Delta#phi (rad);Entries (rad^{-1})", 10, ptasymBins[0], ptasymBins[nPtasymBins], 10, 1.e-1, 1.e+8);
  xjjroot::sethempty(hempty);

  //
  for(int i=0;i<canvdraw.size();i++)
    {
      if(!isMC && i) continue;
      TCanvas* cptasym = new TCanvas("cptasym", "", 600, 600);
      cptasym->SetLogy();
      hempty->Draw();
      Int_t nlegline = 0;
      for(int k=0;k<drawhist[i].size();k++)
        {
          if(histptr.find(drawhist[i].at(k))==histptr.end()) {std::cout<<Form("error: \"%s\" doesn't exist.", drawhist[i].at(k).Data())<<std::endl; return;}
          if(histstyle.find(drawhist[i].at(k))==histstyle.end()) {std::cout<<Form("error: plotting style of \"%s\" is not set.", drawhist[i].at(k).Data())<<std::endl; return;}

          xjjroot::setthgrstyle(histptr[drawhist[i].at(k)], histstyle[drawhist[i].at(k)]);
          histptr[drawhist[i].at(k)]->Draw(Form("%s same", histstyle[drawhist[i].at(k)].GetOption().at(0).Data()));
          if(histstyle[drawhist[i].at(k)].GetOption().at(0).Contains("hist")) continue;
          nlegline++;
        }
      TLegend* leg = new TLegend(0.55, 0.88-0.06*nlegline, 0.82, 0.88, NULL, "brNDC");
      for(int k=0;k<drawhist[i].size();k++) 
        {
          if(histstyle[drawhist[i].at(k)].GetOption().at(0).Contains("hist")) continue;
          leg->AddEntry(histptr[drawhist[i].at(k)], histstyle[drawhist[i].at(k)].GetOption().at(1).Data(), histstyle[drawhist[i].at(k)].GetOption().at(2).Data());
        }
      xjjroot::setlegndraw(leg);
      
      xjjroot::drawCMS(collisionsyst);
      Float_t texypos = 0.91, texxpos = 0.23, texdypos = 0.06;
      xjjroot::drawtex(texxpos, texypos=(texypos-texdypos), Form("|p_{T}^{trk}_{lead D}| > %s GeV/c",xjjc::number_remove_zero(leading_trkptmin).c_str()));
      xjjroot::drawtex(texxpos, texypos=(texypos-texdypos), "|y^{D}| < 1");
      xjjroot::drawtex(texxpos, texypos=(texypos-texdypos), Form("p_{T}^{D}_{lead} > %s GeV/c",xjjc::number_remove_zero(leading_ptmin).c_str()));
      xjjroot::drawtex(texxpos, texypos=(texypos-texdypos), Form("p_{T}^{D}_{sub} > %s GeV/c",xjjc::number_remove_zero(other_ptmin).c_str()));
      
      cptasym->RedrawAxis();
      cptasym->SaveAs(Form("plots/cptasym_%s_%s.pdf", outplotname.Data(), canvdraw[i].Data()));
      delete cptasym;
    }
}

int main(int argc, char* argv[])
{
  if(argc==8)
    {
      ptasymm_plothist(argv[1], argv[2], argv[3], atoi(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));
      return 0;
    }
  else
    {
      return 1;
    }
}

int arguerr(TString collisionsyst, Int_t isMC)
{
  if(collsyst_list.find(collisionsyst)==collsyst_list.end())
    {
      std::cout<<"\033[1;31merror:\033[0m invalid \"collisionsyst\""<<std::endl;
      return 1;
    }
  if(isMC!=0 && isMC!=1)
    {
      std::cout<<"\033[1;31merror:\033[0m invalid \"isMC\""<<std::endl;
      return 1;
    }
  return 0;
}
