#include "dphicor.h"

int dphicor_usehist(TString outfDname, TString outffittpl, TString outplotname, TString collisionsyst, Int_t isMC, Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
{
  xjjroot::setgstyle();

  for(int i=0;i<=nDphiBins;i++) dphiBins[i] = minDphi+i*(maxDphi-minDphi)/nDphiBins;

  TFile* infD = new TFile(Form("%s.root",outfDname.Data()));
  TFile* infS = new TFile(Form("%s.root",outffittpl.Data()));

  TH1D* hdphi[nhist];
  for(int l=0;l<nhist;l++) hdphi[l] = (TH1D*)infD->Get(Form("hdphi_%s",histname[l].Data()));
  TH1D* hdphi_fit[nhist];
  for(int l=0;l<nhist;l++) hdphi_fit[l] = new TH1D(Form("hdphi_%s_fit",histname[l].Data()), ";#Delta#phi (rad);Entries (rad^{-1})", nDphiBins, dphiBins);
  TH1D* hmass[nhist][nDphiBins];
  TH1D** hmassSignal = new TH1D*[nDphiBins];
  TH1D** hmassSwapped = new TH1D*[nDphiBins];

  for(int i=0;i<nDphiBins;i++) 
    {
      hmassSignal[i] = (TH1D*)infS->Get(Form("hmassSignal_%d",i));
      hmassSwapped[i] = (TH1D*)infS->Get(Form("hmassSwapped_%d",i));
      for(int l=0;l<nhist;l++)
        {
          if(!histsave[l]) continue;
          hmass[l][i] = (TH1D*)infD->Get(Form("hmass_%s_%d",histname[l].Data(),i));  
          std::vector<TString> vtex = {TString::Format("%.2f < #Delta#phi < %.2f",dphiBins[i],dphiBins[i+1]), 
                                       TString::Format("|p_{T}^{trk}_{lead D}| > %s GeV/c", xjjc::number_remove_zero(leading_trkptmin).c_str()),
                                       "|y^{D}| < 1",
                                       TString::Format("p_{T}^{D}_{lead} > %s GeV/c",xjjc::number_remove_zero(leading_ptmin).c_str()),
                                       TString::Format("p_{T}^{D} > %s GeV/c",xjjc::number_remove_zero(other_ptmin).c_str())};
          xjjroot::dfitter* dft = new xjjroot::dfitter();
          dft->SetSidebandL(dmass_sideband_l);
          dft->SetSidebandH(dmass_sideband_h);
          dft->fit(hmass[l][i], hmassSignal[i], hmassSwapped[i], collisionsyst, Form("plotfits/cmass_%s_%s_dphi_%d",outplotname.Data(),histname[l].Data(),i), vtex);
          hdphi_fit[l]->SetBinContent(i+1, dft->GetY());
          hdphi_fit[l]->SetBinError(i+1, dft->GetYE());
          delete dft;
        }
    }
  for(int l=0;l<nhist;l++) xjjroot::dividebinwid(hdphi_fit[l]);

  TH1D* hmassSignalLD = (TH1D*)infS->Get("hmassSignalLD");
  TH1D* hmassSwappedLD = (TH1D*)infS->Get("hmassSwappedLD");
  TH1D* hmassLD = (TH1D*)infD->Get("hmassLD");

  xjjroot::dfitter* dftLD = new xjjroot::dfitter("YD");
  dftLD->SetSidebandL(dmass_sideband_l);
  dftLD->SetSidebandH(dmass_sideband_h);
  std::vector<TString> vtexLD = {TString::Format("|p_{T}^{trk}_{lead D}| > %s GeV/c", xjjc::number_remove_zero(leading_trkptmin).c_str()),
                                 "|y^{D}_{lead}| < 1",
                                 TString::Format("p_{T}^{D}_{lead} > %s GeV/c",xjjc::number_remove_zero(leading_ptmin).c_str())};
  dftLD->fit(hmassLD, hmassSignalLD, hmassSwappedLD, collisionsyst, Form("plotfits/cmass_%s",outplotname.Data()), vtexLD);
  Float_t sidebandscale = (dftLD->GetFun_f()->Integral(MASS_DZERO-dmass_sideband_h, MASS_DZERO-dmass_sideband_l) + dftLD->GetFun_f()->Integral(MASS_DZERO+dmass_sideband_l, MASS_DZERO-dmass_sideband_h)) / dftLD->GetFun_not_mass()->Integral(dftLD->GetMassL(), dftLD->GetMassH());
  hdphi[4]->Scale(1./sidebandscale);
  hdphi[5]->Scale(1./sidebandscale);
  hdphi_fit[4]->Scale(1./sidebandscale);

  TH1D* hdphi_subtract_all_fit = (TH1D*)hdphi_fit[0]->Clone("hdphi_subtract_all_fit");
  hdphi_subtract_all_fit->Add(hdphi_fit[4], -1.);


  //
  const int nhistdraw = 5;
  TH1D* histdraw[nhistdraw] = {hdphi[0], hdphi[3], hdphi_subtract_all_fit, hdphi[1], hdphi_fit[0]};
  for(int k=0;k<nhistdraw;k++) xjjroot::setthgrstyle(histdraw[k], histcolor.at(histdraw[k]->GetName()), 20, 1.1);


  //

  TH2F* hempty = new TH2F("hempty", ";#Delta#phi (rad);Entries (rad^{-1})", 10, minDphi, maxDphi, 10, (isMC?std::max(hdphi[3]->GetMinimum(),(double)1):std::max(hdphi_subtract_all_fit->GetMinimum(),(double)1))*0.1, hdphi[0]->GetMaximum()*5);
  xjjroot::sethempty(hempty);
  TCanvas* cdphi = new TCanvas("cdphi","",600,600);
  cdphi->SetLogy();
  hempty->Draw();
  int n_objects = 0;
  for(int k=0;k<nhistdraw;k++)
    {
      histdraw[k]->Draw("pe same");
      n_objects++;
    }
  
  TLegend* leg = new TLegend(0.55, 0.88-0.06*n_objects, 0.82, 0.88,NULL,"brNDC");
  for(int k=0;k<nhistdraw;k++) leg->AddEntry(histdraw[k], histleg.at(histdraw[k]->GetName()), "p");
  xjjroot::setlegndraw(leg);

  xjjroot::drawCMS(collisionsyst);
  Float_t texypos = xjjroot::y_tex_left_top, texxpos = xjjroot::x_tex_left_top;
  xjjroot::drawtex(texxpos, texypos=(texypos-xjjroot::dy_tex_left_top), Form("|p_{T}^{trk}_{lead D}| > %s GeV/c",xjjc::number_remove_zero(leading_trkptmin).c_str()));
  xjjroot::drawtex(texxpos, texypos=(texypos-xjjroot::dy_tex_left_top), "|y^{D}| < 1");
  xjjroot::drawtex(texxpos, texypos=(texypos-xjjroot::dy_tex_left_top), Form("p_{T}^{D}_{lead} > %s GeV/c",xjjc::number_remove_zero(leading_ptmin).c_str()));
  xjjroot::drawtex(texxpos, texypos=(texypos-xjjroot::dy_tex_left_top), Form("p_{T}^{D} > %s GeV/c",xjjc::number_remove_zero(other_ptmin).c_str()));
  cdphi->SaveAs(Form("plots/cdphi_%s.pdf",outplotname.Data()));

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc==9)
    {
      dphicor_usehist(argv[1], argv[2], argv[3], argv[4], atoi(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8]));
      return 0;
    }
  else
    {
      return 1;
    }
}
