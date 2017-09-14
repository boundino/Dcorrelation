#include "ddbar.h"

void ddbar_plothist(TString infname, TString outplotname, 
                    TString collisionsyst, 
                    Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
{
  int arguerr(TString collisionsyst);
  if(arguerr(collisionsyst)) return;

  xjjroot::setgstyle();
  TFile* inf = new TFile(Form("%s.root",infname.Data()));
  if(gethists(std::vector<TFile*>{inf}, "plothist")) return;

  std::vector<TString> canvdraw = {"reco", "gen", "recoNgen"};
  std::vector<std::vector<TString>> drawhist =
    {
      std::vector<TString>{"hdphi_DDbar_Reco", "hdphi_DD_Reco"},
      std::vector<TString>{"hdphi_DDbar_Gen", "hdphi_DD_Gen"},
      std::vector<TString>{"hdphi_DDbar_Reco", "hdphi_DD_Reco", "hdphi_DDbar_Gen", "hdphi_DD_Gen"}
    };
  if(canvdraw.size()!=drawhist.size()) return;

  //
  for(int i=0;i<canvdraw.size();i++)
    {
      TCanvas* cdphi = new TCanvas("cdphi", "", 600, 600);
      // cdphi->SetLogy();
      Float_t yaxismin = 0, yaxismax = 1.;
      TH2F* hempty = new TH2F("hempty", ";#Delta#phi (rad);Entries (rad^{-1})", 10, minDphi, maxDphi, 10, yaxismin, yaxismax);
      xjjroot::sethempty(hempty);
      hempty->Draw();
      TLegend* leg = new TLegend(0.55, 0.88-0.06*drawhist[i].size(), 0.82, 0.88, NULL, "brNDC");
      for(int k=0;k<drawhist[i].size();k++)
        {
          if(histptr.find(drawhist[i].at(k))==histptr.end()) {std::cout<<Form("error: \"%s\" doesn't exist.", drawhist[i].at(k).Data())<<std::endl; return;}
          if(histstyle.find(drawhist[i].at(k))==histstyle.end()) {std::cout<<Form("error: plotting style of \"%s\" is not set.", drawhist[i].at(k).Data())<<std::endl; return;}

          xjjroot::setthgrstyle(histptr[drawhist[i].at(k)], histstyle[drawhist[i].at(k)]);
          histptr[drawhist[i].at(k)]->Draw(Form("%s same", histstyle[drawhist[i].at(k)].GetOption().at(0).Data()));

          leg->AddEntry(histptr[drawhist[i].at(k)], histstyle[drawhist[i].at(k)].GetOption().at(1).Data(), histstyle[drawhist[i].at(k)].GetOption().at(2).Data());
        }
      xjjroot::setlegndraw(leg);

      xjjroot::drawCMS(collisionsyst);
      Float_t texypos = 0.91, texxpos = 0.23, texdypos = 0.06;
      xjjroot::drawtex(texxpos, texypos=(texypos-texdypos), Form("|p_{T}^{trk}_{lead D}| > %s GeV/c",xjjc::number_remove_zero(leading_trkptmin).c_str()));
      xjjroot::drawtex(texxpos, texypos=(texypos-texdypos), "|y^{D}| < 1");
      xjjroot::drawtex(texxpos, texypos=(texypos-texdypos), Form("p_{T}^{D}_{lead} > %s GeV/c",xjjc::number_remove_zero(leading_ptmin).c_str()));
      xjjroot::drawtex(texxpos, texypos=(texypos-texdypos), Form("p_{T}^{D} > %s GeV/c",xjjc::number_remove_zero(other_ptmin).c_str()));
      
      cdphi->RedrawAxis();
      cdphi->SaveAs(Form("plots/cdphi_%s_%s.pdf", outplotname.Data(), canvdraw[i].Data()));
      delete hempty;
      delete cdphi;
    }
}

int main(int argc, char* argv[])
{
  if(argc==7)
    {
      ddbar_plothist(argv[1], argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]));
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
