#include "dphicor.h"

void dphicor_usehist(TString infhistname, TString inftplname, TString outfname, TString outplotname, 
                     TString collisionsyst,  
                     Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
{
  int arguerr(TString collisionsyst);
  if(arguerr(collisionsyst)) return;

  if(initbinning()) return;
  if(createhists("usehist")) return;

  TFile* infS = new TFile(Form("%s.root",inftplname.Data()));
  TFile* infD = new TFile(Form("%s.root",infhistname.Data()));
  if(gethists(std::vector<TFile*>{infS, infD}, "usehist")) return;

  //
  xjjroot::dfitter* dftLD = new xjjroot::dfitter("3SC");
  dftLD->SetTexLinespc(0.009);
  dftLD->SetSidebandL(dmass_sideband_l);
  dftLD->SetSidebandH(dmass_sideband_h);
  std::vector<TString> vtexLD = {TString::Format("|p_{T}^{trk}_{lead D}| > %s GeV/c", xjjc::number_remove_zero(leading_trkptmin).c_str()),
                                 "|y^{D_{lead}}| < 1",
                                 TString::Format("p_{T}^{D}_{lead} > %s GeV/c",xjjc::number_remove_zero(leading_ptmin).c_str())};
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
  Double_t scalefactor  = 1./(1-(N_s_sideband/N_s_total)*(N_b_total/N_b_sideband));
  delete dftLD;

  // 
  for(int l=0;l<nhist;l++)
    {
      if(!histsave[l]) continue;
      for(int i=0;i<nDphiBins;i++) 
        {
          std::vector<TString> vtex = {TString::Format("%.2f < #Delta#phi < %.2f",dphiBins[i],dphiBins[i+1]), 
                                       TString::Format("|p_{T}^{trk}_{lead D}| > %s GeV/c", xjjc::number_remove_zero(leading_trkptmin).c_str()),
                                       "|y^{D}| < 1",
                                       TString::Format("p_{T}^{D}_{lead} > %s GeV/c",xjjc::number_remove_zero(leading_ptmin).c_str()),
                                       TString::Format("p_{T}^{D} > %s GeV/c",xjjc::number_remove_zero(other_ptmin).c_str())};
          xjjroot::dfitter* dft = new xjjroot::dfitter();
          dft->SetSidebandL(dmass_sideband_l);
          dft->SetSidebandH(dmass_sideband_h);
          dft->fit(ahmass[l].at(i), ahmassSignal.at(i), ahmassSwapped.at(i), collisionsyst, Form("plotfits/cmass_%s_%s_dphi_%d",outplotname.Data(),histname[l].Data(),i), vtex);
          ahdphi_fit[l]->SetBinContent(i+1, dft->GetY());
          ahdphi_fit[l]->SetBinError(i+1, dft->GetYE());
          delete dft;
        }
    }
  for(int l=0;l<nhist;l++) ahdphi_fit[l]->Scale(1, "width");

  //
  TH1D* hdphi_subtract_all_fit = (TH1D*)ahdphi_fit[0]->Clone("hdphi_subtract_all_fit");
  hdphi_subtract_all_fit->Add(ahdphi_fit[0], ahdphi_fit[4], scalefactor, -(N_b_total/N_b_sideband)*scalefactor);
  
  TH1D* hdphi_subtract_signal = (TH1D*)ahdphi[1]->Clone("hdphi_subtract_signal");
  hdphi_subtract_signal->Add(ahdphi[1], ahdphi[5], scalefactor, -(N_b_total/N_b_sideband)*scalefactor);  

  TH1D* hdphi_all_all_hist = (TH1D*)ahdphi[0]->Clone("hdphi_all_all_hist");
  TH1D* hdphi_all_signal_hist = (TH1D*)ahdphi[1]->Clone("hdphi_all_signal_hist");
  TH1D* hdphi_signal_signal_hist = (TH1D*)ahdphi[3]->Clone("hdphi_signal_signal_hist");

  //
  ahdphi_plot.push_back(ahdphi[0]);
  ahdphi_plot.push_back(ahdphi[1]);
  ahdphi_plot.push_back(ahdphi[3]);
  ahdphi_plot.push_back(hdphi_all_all_hist);
  ahdphi_plot.push_back(hdphi_all_signal_hist);
  ahdphi_plot.push_back(hdphi_signal_signal_hist);
  ahdphi_plot.push_back(ahdphi_fit[0]);
  ahdphi_plot.push_back(hdphi_subtract_all_fit);
  ahdphi_plot.push_back(hdphi_subtract_signal);

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
      dphicor_usehist(argv[1], argv[2], argv[3], argv[4], argv[5], atof(argv[6]), atof(argv[7]), atof(argv[8]));
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
