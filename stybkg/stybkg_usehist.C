#include "stybkg.h"

void stybkg_usehist(TString infhistname, TString inftplname, TString outfname, TString outplotname, 
                    TString collisionsyst, 
                    Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
{
  int arguerr(TString collisionsyst);
  if(arguerr(collisionsyst)) return;

  if(initbinning()) return;

  TFile* infS = new TFile(Form("%s.root",inftplname.Data()));
  TFile* infD = new TFile(Form("%s.root",infhistname.Data()));
  if(gethists(std::vector<TFile*>{infS, infD}, "usehist")) return;

  // 
  xjjroot::dfitter* dftLD = new xjjroot::dfitter("3SC");
  dftLD->SetTexLinespc(0.009);
  dftLD->SetSidebandL(dmass_sideband_l);
  dftLD->SetSidebandH(dmass_sideband_h);
  std::vector<TString> vtexLD = 
    {
      // TString::Format("|p_{T}^{trk}_{lead D}| > %s GeV/c", xjjc::number_remove_zero(leading_trkptmin).c_str()),
      "|y^{D_{lead}}| < 1",
      TString::Format("p_{T}^{D}_{lead} > %s GeV/c",xjjc::number_remove_zero(leading_ptmin).c_str())
    };
  dftLD->fit(hmassLD, hmassSignalLD, hmassSwappedLD, collisionsyst, Form("plotfits/cmass_LD_%s",outplotname.Data()), vtexLD);

  Double_t N_signal_total    = dftLD->GetFun_mass()->Integral(dftLD->GetMassL(), dftLD->GetMassH())/binwid_hist_dzero;
  Double_t N_swap_total      = dftLD->GetFun_swap()->Integral(dftLD->GetMassL(), dftLD->GetMassH())/binwid_hist_dzero;
  Double_t N_comb_total      = dftLD->GetFun_background()->Integral(dftLD->GetMassL(), dftLD->GetMassH())/binwid_hist_dzero;
  Double_t N_signal_sideband = dftLD->GetFun_mass()->Integral(MASS_DZERO-dmass_sideband_h, MASS_DZERO-dmass_sideband_l)/binwid_hist_dzero + dftLD->GetFun_mass()->Integral(MASS_DZERO+dmass_sideband_l, MASS_DZERO+dmass_sideband_h)/binwid_hist_dzero;
  Double_t N_swap_sideband   = dftLD->GetFun_swap()->Integral(MASS_DZERO-dmass_sideband_h, MASS_DZERO-dmass_sideband_l)/binwid_hist_dzero + dftLD->GetFun_swap()->Integral(MASS_DZERO+dmass_sideband_l, MASS_DZERO+dmass_sideband_h)/binwid_hist_dzero;
  Double_t N_comb_sideband   = dftLD->GetFun_background()->Integral(MASS_DZERO-dmass_sideband_h, MASS_DZERO-dmass_sideband_l)/binwid_hist_dzero + dftLD->GetFun_background()->Integral(MASS_DZERO+dmass_sideband_l, MASS_DZERO+dmass_sideband_h)/binwid_hist_dzero;

  Double_t N_s_total    = N_signal_total + N_swap_total;
  Double_t N_s_sideband = N_signal_sideband + N_swap_sideband;
  Double_t N_b_total    = N_comb_total;
  Double_t N_b_sideband = N_comb_sideband;
  Double_t scalefactor = 1./(1-(N_s_sideband/N_s_total)*(N_b_total/N_b_sideband));
  delete dftLD;

  //
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

  ahdphi_plot.push_back(hdphi_incl_signalNswapNcomb_signal_hist);
  ahdphi_plot.push_back(hdphi_incl_signalNswap_signal_hist);
  ahdphi_plot.push_back(hdphi_incl_signal_signal_hist);
  ahdphi_plot.push_back(ahdphi[0]);
  ahdphi_plot.push_back(hdphi_incl_comb_signal_norm);
  ahdphi_plot.push_back(hdphi_incl_swap_signal_norm);
  ahdphi_plot.push_back(hdphi_incl_signal_signal_norm);
  ahdphi_plot.push_back(hdphi_incl_signalNswap_signal_norm_hist);
  ahdphi_plot.push_back(hdphi_sideband_signalNswap_signal_norm);
  ahdphi_plot.push_back(hdphi_incl_comb_signal_norm_hist);
  ahdphi_plot.push_back(hdphi_sideband_comb_signal_norm);
  ahdphi_plot.push_back(hdphi_incl_all_signal_subtract);

  //
  TFile* outf = new TFile(Form("%s.root",outfname.Data()),"recreate");
  outf->cd();
  if(writehists("usehist")) return;
  outf->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==9)
    {
      stybkg_usehist(argv[1], argv[2], argv[3], argv[4], argv[5], atof(argv[6]), atof(argv[7]), atof(argv[8]));
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
