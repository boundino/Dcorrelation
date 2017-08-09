#include "stybkg.h"

int stybkg_usehist(TString outfDname, TString outffittpl, TString outplotname, TString collisionsyst, Int_t isMC, Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
{
  xjjroot::setgstyle();

  initbinning();

  TFile* infD = new TFile(Form("%s.root",outfDname.Data()));
  TFile* infS = new TFile(Form("%s.root",outffittpl.Data()));

  // leading D bkg subtraction
  TH1D* hmassSignalLD = (TH1D*)infS->Get("hmassSignalLD");
  TH1D* hmassSwappedLD = (TH1D*)infS->Get("hmassSwappedLD");
  TH1D* hmassLD = (TH1D*)infD->Get("hmassLD");

  xjjroot::dfitter* dftLD = new xjjroot::dfitter("3SC");
  dftLD->SetTexLinespc(0.009);
  dftLD->SetSidebandL(dmass_sideband_l);
  dftLD->SetSidebandH(dmass_sideband_h);
  std::vector<TString> vtexLD = {TString::Format("|p_{T}^{trk}_{lead D}| > %s GeV/c", xjjc::number_remove_zero(leading_trkptmin).c_str()),
                                 "|y^{D_{lead}}| < 1",
                                 TString::Format("p_{T}^{D}_{lead} > %s GeV/c",xjjc::number_remove_zero(leading_ptmin).c_str())};
  dftLD->fit(hmassLD, hmassSignalLD, hmassSwappedLD, collisionsyst, Form("plotfits/cmass_%s",outplotname.Data()), vtexLD);
  Float_t sidebandscale = dftLD->GetFun_not_mass()->Integral(dftLD->GetMassL(), dftLD->GetMassH()) / (dftLD->GetFun_f()->Integral(MASS_DZERO-dmass_sideband_h, MASS_DZERO-dmass_sideband_l) + dftLD->GetFun_f()->Integral(MASS_DZERO+dmass_sideband_l, MASS_DZERO+dmass_sideband_h));

  // preparing hists
  TH1D* ahdphi[nhist];
  for(int l=0;l<nhist;l++) ahdphi[l] = (TH1D*)infD->Get(Form("hdphi_%s",histname[l].Data()));

  TH1D* hdphi_incl_signal_signal_hist = (TH1D*)ahdphi[1]->Clone("hdphi_incl_signal_signal_hist");
  TH1D* hdphi_incl_signalNswap_signal_hist = (TH1D*)hdphi_incl_signal_signal_hist->Clone("hdphi_incl_signalNswap_signal_hist");
  hdphi_incl_signalNswap_signal_hist->Add(ahdphi[2]);
  TH1D* hdphi_incl_signalNswapNcomb_signal_hist = (TH1D*)hdphi_incl_signalNswap_signal_hist->Clone("hdphi_incl_signalNswapNcomb_signal_hist");
  hdphi_incl_signalNswapNcomb_signal_hist->Add(ahdphi[3]);  

  TH1D* hdphi_incl_notmass_signal_norm = (TH1D*)ahdphi[4]->Clone("hdphi_incl_notmass_signal_norm");
  hdphi_incl_notmass_signal_norm->Scale(1./hdphi_incl_notmass_signal_norm->Integral());
  TH1D* hdphi_sideband_all_signal_norm = (TH1D*)ahdphi[5]->Clone("hdphi_sideband_all_signal_norm");
  hdphi_sideband_all_signal_norm->Scale(1./hdphi_sideband_all_signal_norm->Integral());

  TH1D* hdphi_incl_signal_signal_norm = (TH1D*)ahdphi[1]->Clone("hdphi_incl_signal_signal_norm");
  hdphi_incl_signal_signal_norm->Scale(1./hdphi_incl_signal_signal_norm->Integral());
  TH1D* hdphi_incl_swap_signal_norm = (TH1D*)ahdphi[2]->Clone("hdphi_incl_swap_signal_norm");
  hdphi_incl_swap_signal_norm->Scale(1./hdphi_incl_swap_signal_norm->Integral());
  TH1D* hdphi_incl_comb_signal_norm = (TH1D*)ahdphi[3]->Clone("hdphi_incl_comb_signal_norm");
  hdphi_incl_comb_signal_norm->Scale(1./hdphi_incl_comb_signal_norm->Integral());

  TH1D* hdphi_sideband_all_signal_scale = (TH1D*)ahdphi[5]->Clone("hdphi_sideband_all_signal_scale");
  hdphi_sideband_all_signal_scale->Scale(sidebandscale);

  TH1D* hdphi_incl_all_signal_subtract = (TH1D*)ahdphi[0]->Clone("hdphi_incl_all_signal_subtract");
  hdphi_incl_all_signal_subtract->Add(hdphi_sideband_all_signal_scale, -1);

  // all hists are prepared
  std::vector<TH1D*> ahistdraw = 
    {
      hdphi_incl_signalNswapNcomb_signal_hist,  hdphi_incl_signalNswap_signal_hist,  hdphi_incl_signal_signal_hist,  
      ahdphi[0],                                hdphi_incl_notmass_signal_norm,      hdphi_sideband_all_signal_norm,
      hdphi_incl_comb_signal_norm,              hdphi_incl_swap_signal_norm,         hdphi_incl_signal_signal_norm,
      ahdphi[4],                                hdphi_sideband_all_signal_scale,     hdphi_incl_all_signal_subtract
    };
  const int nhistdraw = ahistdraw.size();
  Float_t yaxismin = ahistdraw[0]->GetMinimum(), yaxismax = ahistdraw[0]->GetMaximum();
  for(int k=0;k<nhistdraw;k++) 
    {
      if(histstyle.find(ahistdraw[k]->GetName())==histstyle.end()) return 1;
      xjjroot::setthgrstyle(ahistdraw[k], histstyle.at(ahistdraw[k]->GetName()));
      if(ahistdraw[k]->GetMinimum()<yaxismin) yaxismin = ahistdraw[k]->GetMinimum();
      if(ahistdraw[k]->GetMaximum()>yaxismax) yaxismax = ahistdraw[k]->GetMaximum();
    }

  //  
  // TH2F* hempty = new TH2F("hempty", ";#Delta#phi (rad);Entries (rad^{-1})", 10, minDphi, maxDphi, 10, (yaxismin>0?yaxismin:1)*1.e-1, yaxismax*5.e+1);
  TH2F* hempty = new TH2F("hempty", ";#Delta#phi (rad);Entries (rad^{-1})", 10, minDphi, maxDphi, 10, 1.e0, 5.e+6);
  xjjroot::sethempty(hempty);
  TH2F* hempty_norm = new TH2F("hempty_norm", ";#Delta#phi (rad);Probability (rad^{-1})", 10, minDphi, maxDphi, 10, 1.e-5, 5.e+2);
  xjjroot::sethempty(hempty_norm);

  TString  canvdraw[]  =  {"xcheck",  "directcompshape",  "components",  "directcomp",  "subtract"};
  bool     ifnorm[]    =  {false,     true,               true,          false,         false};
  bool ifdrawhist[][nhistdraw] = 
    {
      {true,   true,   true,   true,   false,  false,  false,  false,  false,  false,  false,  false},
      {false,  false,  false,  false,  true,   true,   false,  false,  false,  false,  false,  false},
      {false,  false,  false,  false,  false,  false,  true,   true,   true,   false,  false,  false},  
      {false,  false,  false,  false,  false,  false,  false,  false,  false,  true,   true,   false},
      {false,  false,  true,   true,   false,  false,  false,  false,  false,  false,  false,  true}
    };
  const Int_t ncanvdraw = sizeof(canvdraw)/sizeof(canvdraw[0]);
  
  //
  for(int i=0;i<ncanvdraw;i++)
    {
      TCanvas* cdphi = new TCanvas("cdphi","",600,600);
      cdphi->SetLogy();
      if(ifnorm[i]) hempty_norm->Draw();
      else hempty->Draw();
      int n_objects = 0;
      for(int k=0;k<nhistdraw;k++)
        {
          if(!ifdrawhist[i][k]) continue;
          ahistdraw[k]->Draw(Form("%s same", histstyle.at(ahistdraw[k]->GetName()).GetOption().Data()));
          n_objects++;
        }
      TLegend* leg = new TLegend(0.55, 0.88-0.06*n_objects, 0.82, 0.88, NULL, "brNDC");
      for(int k=0;k<nhistdraw;k++) 
        {
          if(!ifdrawhist[i][k]) continue;
          if(histleg.find(ahistdraw[k]->GetName())==histleg.end()) return 1;
          TString legopt = histstyle.at(ahistdraw[k]->GetName()).GetOption().Contains("hist")?"f":"p";
          leg->AddEntry(ahistdraw[k], histleg.at(ahistdraw[k]->GetName()), legopt.Data());
        }
      xjjroot::setlegndraw(leg);
      
      xjjroot::drawCMS(collisionsyst);
      Float_t texypos = xjjroot::y_tex_left_top, texxpos = xjjroot::x_tex_left_top;
      xjjroot::drawtex(texxpos, texypos=(texypos-xjjroot::dy_tex_left_top), Form("|p_{T}^{trk}_{lead D}| > %s GeV/c",xjjc::number_remove_zero(leading_trkptmin).c_str()));
      xjjroot::drawtex(texxpos, texypos=(texypos-xjjroot::dy_tex_left_top), "|y^{D}| < 1");
      xjjroot::drawtex(texxpos, texypos=(texypos-xjjroot::dy_tex_left_top), Form("p_{T}^{D}_{lead} > %s GeV/c",xjjc::number_remove_zero(leading_ptmin).c_str()));
      xjjroot::drawtex(texxpos, texypos=(texypos-xjjroot::dy_tex_left_top), Form("p_{T}^{D} > %s GeV/c",xjjc::number_remove_zero(other_ptmin).c_str()));

      cdphi->RedrawAxis();
      cdphi->SaveAs(Form("plots/cdphi_%s_%s.pdf", outplotname.Data(), canvdraw[i].Data()));
      delete cdphi;
    }

  //

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc==9)
    {
      stybkg_usehist(argv[1], argv[2], argv[3], argv[4], atoi(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8]));
      return 0;
    }
  else
    {
      return 1;
    }
}
