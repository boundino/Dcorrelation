#include "stybkg.h"

int stybkg_savehist(TString infname, TString outfname, TString collisionsyst, Int_t isMC, Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
{
  std::cout<<std::endl;
  initcutval(collisionsyst);
  
  initbinning();
  
  TFile* inf = new TFile(infname);
  TTree* ntDkpi = (TTree*)inf->Get("ntDkpi");
  TTree* ntGen = (TTree*)inf->Get("ntGen");
  TTree* ntHlt = (TTree*)inf->Get("ntHlt");
  TTree* ntHi = (TTree*)inf->Get("ntHi");
  TTree* ntSkim = (TTree*)inf->Get("ntSkim");
  
  readD dcand(MAX_XB, MAX_GEN, (bool)isMC);
  dcand.setbranchesaddress(ntDkpi, ntGen, ntHi);
  ntHlt->SetBranchStatus("*", 0);
  int val_hlt;
  xjjroot::setbranchaddress(ntHlt, cutval_hlt.Data(), &val_hlt);
  ntSkim->SetBranchStatus("*", 0);
  int* val_skim = new int[cutval_skim.size()];
  for(int k=0;k<cutval_skim.size();k++) xjjroot::setbranchaddress(ntSkim, cutval_skim[k].Data(), &val_skim[k]);
  
  TH1D* ahdphi[nhist];
  for(int l=0;l<nhist;l++) 
    {
      ahdphi[l] = new TH1D(Form("hdphi_%s",histname[l].Data()), ";#Delta#phi (rad);Entries (rad^{-1})", nDphiBins, dphiBins);
    }
  TH1D* hmassLD = new TH1D("hmassLD", ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
  
  int nentries = ntDkpi->GetEntries();
  for(int i=0;i<nentries;i++)
    {
      if(i%10000==0) xjjc::progressbar(i, nentries);
      ntDkpi->GetEntry(i);
      ntHlt->GetEntry(i);
      ntSkim->GetEntry(i);

      // if(!val_hlt && !isMC) continue;
      for(int k=0;k<cutval_skim.size();k++)
        {
          if(!val_skim[k]) continue;
        }
      if(TMath::Abs(dcand.PVz)>15) continue;
      
      // find leading D
      int jleading = -1;
      double ptleading = 0;
      std::map<int, double> dphi[nhist];
      for(int j=0;j<dcand.Dsize;j++)
        {
          int ipt = xjjc::findibin(&ptBins, (double)dcand.Dpt[j]);
          if(ipt<0) continue;
          int err_initcutval_ptdep = initcutval_ptdep(collisionsyst, ipt);
          if(err_initcutval_ptdep) return 1;
          dcand.settrkcut(leading_trkptmin, cutval_trkEta, cutval_trkPtErr);
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, leading_ptmin);
          if(dcand.isselected(j) && dcand.Dmass[j] > 1.7 && dcand.Dmass[j] < 2.0)
            {
              if(dcand.Dpt[j]>ptleading)
                {
                  jleading = j;
                  ptleading = dcand.Dpt[j];
                }
            }          
          dcand.settrkcut(cutval_trkPt, cutval_trkEta, cutval_trkPtErr);
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, other_ptmin);
          if(!dcand.isselected(j) || dcand.Dgen[j]!=23333) continue;
          for(int l=0;l<nhist;l++) dphi[l].insert(std::pair<int, double>(j, dcand.Dphi[j]));
        }

      // fill dphi
      if(jleading<0) continue;
      hmassLD->Fill(dcand.Dmass[jleading]);
      Bool_t leadingsel[nhist] = {
        true, 
        dcand.Dgen[jleading]==23333, 
        dcand.Dgen[jleading]==23344, 
        dcand.Dgen[jleading]!=23333 && dcand.Dgen[jleading]!=23344, 
        dcand.Dgen[jleading]!=23333, 
        TMath::Abs(dcand.Dmass[jleading]-MASS_DZERO)>dmass_sideband_l && TMath::Abs(dcand.Dmass[jleading]-MASS_DZERO)<dmass_sideband_h
      };
      for(int l=0;l<nhist;l++)
        {
          if(!leadingsel[l]) continue;
          if(dphi[l].empty()) continue;
          for(std::map<int, double>::iterator it=dphi[l].begin(); it!=dphi[l].end(); it++)
            {
              if(it->first==jleading) continue;
              double deltaphi = TMath::Abs(it->second - dcand.Dphi[jleading]);
              double filldeltaphi = deltaphi<M_PI?deltaphi:(2*M_PI-deltaphi);              
              ahdphi[l]->Fill(filldeltaphi);
            }
        }
    }
  xjjc::progressbar_summary(nentries);
  
  for(int l=0;l<nhist;l++) xjjroot::dividebinwid(ahdphi[l]);

  TFile* outf = new TFile(Form("%s.root",outfname.Data()),"recreate");
  outf->cd();
  hmassLD->Write();
  for(int l=0;l<nhist;l++)
    {
      ahdphi[l]->Write();
    }
  outf->Write();
  outf->Close();

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc==8)
    {
      stybkg_savehist(argv[1], argv[2], argv[3], atoi(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));
      return 0;
    }
  else
    {
      return 1;
    }
}
