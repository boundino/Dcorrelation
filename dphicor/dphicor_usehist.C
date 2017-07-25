#include "dphicor.h"

int dphicor_usehist(TString outfDname, TString outffittpl, TString outplotname, TString collisionsyst, Int_t isMC, Float_t leading_ptmin, Float_t other_ptmin)
{
  xjjrootuti::setgstyle();

  for(int i=0;i<=nDphiBins;i++) dphiBins[i] = minDphi+i*(maxDphi-minDphi)/nDphiBins;

  TFile* infD = new TFile(Form("%s.root",outfDname.Data()));
  TFile* infS = new TFile(Form("%s.root",outffittpl.Data()));
  TH1D** hdphi = new TH1D*[nhist];
  for(int l=0;l<nhist;l++) 
    {
      hdphi[l] = (TH1D*)infD->Get(histname[l]);
      hdphi[l]->Sumw2();
      xjjrootuti::setthgr(hdphi[l]);
      xjjrootuti::setthgrstyle(hdphi[l], hcolor[l]);
    }

  TH1D** hmass = new TH1D*[nDphiBins];
  TH1D** hmassSignal = new TH1D*[nDphiBins];
  TH1D** hmassSwapped = new TH1D*[nDphiBins];
  Float_t* fitresult = new Float_t[2];
  TH1D* hdphi_all_fit = new TH1D("hdphi_all_fit",";#Delta#phi;Entries", nDphiBins, minDphi, maxDphi);
  xjjrootuti::setthgr(hdphi_all_fit);
  xjjrootuti::setthgrstyle(hdphi_all_fit, kRed+1);

  for(int i=0;i<nDphiBins;i++) 
    {
      hmass[i] = (TH1D*)infD->Get(Form("hmass_%d",i));  
      hmassSignal[i] = (TH1D*)infS->Get(Form("hmassSignal_%d",i));
      hmassSwapped[i] = (TH1D*)infS->Get(Form("hmassSwapped_%d",i));
      Float_t dphimin = dphiBins[i], dphimax = dphiBins[i+1];
      TF1* fmass = fit(hmass[i], hmassSignal[i], hmassSwapped[i], 
                       Form("plotfits/cmass_%s_dphi_%d",outplotname.Data(),i), fitresult, collisionsyst,
                       leading_ptmin, other_ptmin, dphimin, dphimax,
                       false);
      hdphi_all_fit->SetBinContent(i+1, fmass->Integral(minhisto,maxhisto)/binwidthmass);
      hdphi_all_fit->SetBinError(i+1, fmass->Integral(minhisto,maxhisto)/binwidthmass*fmass->GetParError(0)/fmass->GetParameter(0));
      delete fmass;
    }
  hdphi_all_fit->Scale(nDphiBins*1./nDphiBins_fine);

  TH2F* hempty = new TH2F("hempty", ";#Delta#phi;Entries", nDphiBins, minDphi, maxDphi, 10, (isMC?std::max(hdphi[1]->GetMinimum(),(double)1):std::max(hdphi_all_fit->GetMinimum(),(double)1))*0.1, hdphi[0]->GetMaximum()*5);
  xjjrootuti::sethempty(hempty);
  TCanvas* cdphi = new TCanvas("cdphi","",600,600);
  cdphi->SetLogy();
  hempty->Draw();
  int n_objects = 0;
  for(int l=0;l<nhist;l++)
    {
      if(l && (l==3 || !isMC)) continue;
      hdphi[l]->Draw("pe same");
      n_objects++;
    }
  hdphi_all_fit->Draw("pe same");
  n_objects++;
  
  TLegend* leg = new TLegend(0.55, 0.88-0.06*n_objects, 0.82, 0.88,NULL,"brNDC");
  for(int l=0;l<nhist;l++)
    {
      if(l && (l==3 || !isMC)) continue;
      leg->AddEntry(hdphi[l], histleg[l], "p");
    }
  leg->AddEntry(hdphi_all_fit, "all D_{lead}, fit ext D", "p");
  xjjrootuti::setlegndraw(leg);

  xjjrootuti::drawCMS(collisionsyst);
  float texypos = xjjrootuti::y_tex_left_top, texxpos = xjjrootuti::x_tex_left_top;
  TLatex* texrap = new TLatex(texxpos,texypos=(texypos-xjjrootuti::dy_tex_left_top),"|y^{D}| < 1");
  xjjrootuti::settexndraw(texrap);
  TLatex* texleadpt = new TLatex(texxpos,texypos=(texypos-xjjrootuti::dy_tex_left_top),Form("p_{T}^{D}_{lead} > %s GeV/c",xjjuti::number_remove_zero(leading_ptmin).c_str()));
  xjjrootuti::settexndraw(texleadpt);
  TLatex* texpt = new TLatex(texxpos ,texypos=(texypos-xjjrootuti::dy_tex_left_top),Form("p_{T}^{D} > %s GeV/c",xjjuti::number_remove_zero(other_ptmin).c_str()));
  xjjrootuti::settexndraw(texpt);
  cdphi->SaveAs(Form("plots/cdphi_%s.pdf",outplotname.Data()));

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc==8)
    {
      dphicor_usehist(argv[1], argv[2], argv[3], argv[4], atoi(argv[5]), atof(argv[6]), atof(argv[7]));
      return 0;
    }
  else
    {
      return 1;
    }
}
