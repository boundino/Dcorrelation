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
  dftLD->fit(hmassLD, hmassSignalLD, hmassSwappedLD, collisionsyst, Form("plotfits/cmass_LD_%s",outplotname.Data()), vtexLD);
  Double_t N_signal_total = dftLD->GetFun_mass()->Integral(dftLD->GetMassL(), dftLD->GetMassH())/binwid_hist_dzero;
  Double_t N_swap_total = dftLD->GetFun_swap()->Integral(dftLD->GetMassL(), dftLD->GetMassH())/binwid_hist_dzero;
  Double_t N_comb_total = dftLD->GetFun_background()->Integral(dftLD->GetMassL(), dftLD->GetMassH())/binwid_hist_dzero;
  Double_t N_signal_sideband = dftLD->GetFun_mass()->Integral(MASS_DZERO-dmass_sideband_h, MASS_DZERO-dmass_sideband_l)/binwid_hist_dzero + dftLD->GetFun_mass()->Integral(MASS_DZERO+dmass_sideband_l, MASS_DZERO+dmass_sideband_h)/binwid_hist_dzero;
  Double_t N_swap_sideband = dftLD->GetFun_swap()->Integral(MASS_DZERO-dmass_sideband_h, MASS_DZERO-dmass_sideband_l)/binwid_hist_dzero + dftLD->GetFun_swap()->Integral(MASS_DZERO+dmass_sideband_l, MASS_DZERO+dmass_sideband_h)/binwid_hist_dzero;
  Double_t N_comb_sideband = dftLD->GetFun_background()->Integral(MASS_DZERO-dmass_sideband_h, MASS_DZERO-dmass_sideband_l)/binwid_hist_dzero + dftLD->GetFun_background()->Integral(MASS_DZERO+dmass_sideband_l, MASS_DZERO+dmass_sideband_h)/binwid_hist_dzero;

  Double_t N_s_total = N_signal_total + N_swap_total;
  Double_t N_s_sideband = N_signal_sideband + N_swap_sideband;
  Double_t N_b_total = N_comb_total;
  Double_t N_b_sideband = N_comb_sideband;

  Double_t scalefactor = 1./(1-(N_s_sideband/N_s_total)*(N_b_total/N_b_sideband));

  // preparing hists
  TH1D* ahdphi[nhist];
  for(int l=0;l<nhist;l++) ahdphi[l] = (TH1D*)infD->Get(Form("hdphi_%s",histname[l].Data()));

  TH1D* hdphi_incl_signal_signal_hist = (TH1D*)ahdphi[1]->Clone("hdphi_incl_signal_signal_hist");
  TH1D* hdphi_incl_signalNswap_signal_hist = (TH1D*)hdphi_incl_signal_signal_hist->Clone("hdphi_incl_signalNswap_signal_hist");
  hdphi_incl_signalNswap_signal_hist->Add(ahdphi[2]);
  TH1D* hdphi_incl_signalNswapNcomb_signal_hist = (TH1D*)hdphi_incl_signalNswap_signal_hist->Clone("hdphi_incl_signalNswapNcomb_signal_hist");
  hdphi_incl_signalNswapNcomb_signal_hist->Add(ahdphi[3]);  

  TH1D* hdphi_incl_signal_signal_norm = (TH1D*)ahdphi[1]->Clone("hdphi_incl_signal_signal_norm");
  hdphi_incl_signal_signal_norm->Scale(1./N_signal_total);
  TH1D* hdphi_incl_swap_signal_norm = (TH1D*)ahdphi[2]->Clone("hdphi_incl_swap_signal_norm");
  hdphi_incl_swap_signal_norm->Scale(1./N_swap_total);
  TH1D* hdphi_incl_comb_signal_norm = (TH1D*)ahdphi[3]->Clone("hdphi_incl_comb_signal_norm");
  hdphi_incl_comb_signal_norm->Scale(1./N_comb_total);

  TH1D* hdphi_sideband_signalNswap_signal_norm = (TH1D*)ahdphi[5]->Clone("hdphi_sideband_signalNswap_signal_norm");
  hdphi_sideband_signalNswap_signal_norm->Add(ahdphi[6]);
  hdphi_sideband_signalNswap_signal_norm->Scale(1./N_s_sideband);
  TH1D* hdphi_sideband_comb_signal_norm = (TH1D*)ahdphi[7]->Clone("hdphi_sideband_comb_signal_norm");
  hdphi_sideband_comb_signal_norm->Scale(1./N_b_sideband);
  TH1D* hdphi_incl_signalNswap_signal_norm_hist = (TH1D*)ahdphi[1]->Clone("hdphi_incl_signalNswap_signal_norm_hist");
  hdphi_incl_signalNswap_signal_norm_hist->Add(ahdphi[2]);
  hdphi_incl_signalNswap_signal_norm_hist->Scale(1./N_s_total);
  TH1D* hdphi_incl_comb_signal_norm_hist = (TH1D*)ahdphi[3]->Clone("hdphi_incl_comb_signal_norm_hist");
  hdphi_incl_comb_signal_norm_hist->Scale(1./N_b_total);

  TH1D* hdphi_incl_all_signal_subtract = (TH1D*)ahdphi[0]->Clone("hdphi_incl_all_signal_subtract");
  hdphi_incl_all_signal_subtract->Add(ahdphi[0], ahdphi[4], scalefactor, -scalefactor*N_b_total/N_b_sideband);

  // all hists are prepared
  std::vector<TH1D*> ahistdraw = 
    {
      hdphi_incl_signalNswapNcomb_signal_hist,  hdphi_incl_signalNswap_signal_hist,      hdphi_incl_signal_signal_hist,  
      ahdphi[0],                                
      hdphi_incl_comb_signal_norm,              hdphi_incl_swap_signal_norm,             hdphi_incl_signal_signal_norm,
      hdphi_incl_signalNswap_signal_norm_hist,  hdphi_sideband_signalNswap_signal_norm,  
      hdphi_incl_comb_signal_norm_hist,         hdphi_sideband_comb_signal_norm,         
      hdphi_incl_all_signal_subtract
    };
  const int nhistdraw = ahistdraw.size();
  for(int k=0;k<nhistdraw;k++) 
    {
      if(histstyle.find(ahistdraw[k]->GetName())==histstyle.end()) return 1;
      xjjroot::setthgrstyle(ahistdraw[k], histstyle.at(ahistdraw[k]->GetName()));
    }

  //
  TString  canvdraw[]  =  {"xcheck", "components", "inclvssideband_s", "inclvssideband_b", "subtract"};
  const Int_t ncanvdraw = sizeof(canvdraw)/sizeof(canvdraw[0]);
  bool ifdrawhist[ncanvdraw][nhistdraw] = 
    {
      {true,   true,   true,   true,   false,  false,  false,  false,  false,  false,  false,  false},
      {false,  false,  false,  false,  true,   true,   true,   false,  false,  false,  false,  false},
      {false,  false,  false,  false,  false,  false,  false,  true,   true,   false,  false,  false},  
      {false,  false,  false,  false,  false,  false,  false,  false,  false,  true,   true,   false},
      {false,  true,   true,   true,   false,  false,  false,  false,  false,  false,  false,  true}
    };
  
  //
  for(int i=0;i<ncanvdraw;i++)
    {
      TCanvas* cdphi = new TCanvas("cdphi","",600,600);
      cdphi->SetLogy();
      Float_t yaxismin = 1.e+30, yaxismax = -1;
      for(int k=0;k<nhistdraw;k++) 
        {
          if(!ifdrawhist[i][k]) continue;
          if(ahistdraw[k]->GetMinimum()<yaxismin) yaxismin = ahistdraw[k]->GetMinimum();
          if(ahistdraw[k]->GetMaximum()>yaxismax) yaxismax = ahistdraw[k]->GetMaximum();
        }
      TH2F* hempty = new TH2F("hempty", ";#Delta#phi (rad);Entries (rad^{-1})", 10, minDphi, maxDphi, 10, (yaxismin>0?yaxismin:1)*1.e-1, yaxismax*1.e+2);
      xjjroot::sethempty(hempty);
      hempty->Draw();
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
      delete hempty;
      delete cdphi;
    }

  std::cout<<hdphi_incl_all_signal_subtract->GetBinContent(1)<<" "<<hdphi_incl_all_signal_subtract->GetBinError(1)<<std::endl;  

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
