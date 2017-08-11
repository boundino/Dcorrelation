#include "stybkg.h"

int stybkg_savefittpl(TString infname, TString outfname, TString collisionsyst, Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
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

  readD dcand(MAX_XB, MAX_GEN, true);
  dcand.setbranchesaddress(ntDkpi, ntGen, ntHi);
  ntHlt->SetBranchStatus("*", 0);
  int val_hlt;
  xjjroot::setbranchaddress(ntHlt, cutval_hlt.Data(), &val_hlt);
  ntSkim->SetBranchStatus("*", 0);
  int* val_skim = new int[cutval_skim.size()];
  for(int k=0;k<cutval_skim.size();k++) xjjroot::setbranchaddress(ntSkim, cutval_skim[k].Data(), &val_skim[k]);

  TH1D** hmassSignal = new TH1D*[nDphiBins];
  TH1D** hmassSwapped = new TH1D*[nDphiBins];
  TH1D* hmassSignalLD = new TH1D("hmassSignalLD", ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
  TH1D* hmassSwappedLD = new TH1D("hmassSwappedLD", ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
  for(int i=0;i<nDphiBins;i++)
    {
      hmassSignal[i] = new TH1D(Form("hmassSignal_%d",i), ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
      hmassSwapped[i] = new TH1D(Form("hmassSwapped_%d",i), ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);
    }
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
      
      int jleading = -1;
      double ptleading = 0;
      std::map<int, double> dphi;
      for(int j=0;j<dcand.Dsize;j++)
        {
          int ipt = xjjc::findibin(&ptBins, (double)dcand.Dpt[j]);
          if(ipt<0) continue;
          int err_initcutval_ptdep = initcutval_ptdep(collisionsyst, ipt);
          if(err_initcutval_ptdep) return 1;
          dcand.settrkcut(leading_trkptmin, cutval_trkEta, cutval_trkPtErr);
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, leading_ptmin);
          if(dcand.isselected(j))
            {
              if(dcand.Dpt[j]>ptleading)
                {
                  jleading = j;
                  ptleading = dcand.Dpt[j];
                }
            }
          dcand.settrkcut(cutval_trkPt, cutval_trkEta, cutval_trkPtErr);
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, other_ptmin);
          if(!dcand.isselected(j)) continue;
          dphi.insert(std::pair<int, double>(j, dcand.Dphi[j]));
        }
      if(jleading<0) continue;
      if(dcand.Dgen[jleading]==23333) hmassSignalLD->Fill(dcand.Dmass[jleading]);
      if(dcand.Dgen[jleading]==23344) hmassSwappedLD->Fill(dcand.Dmass[jleading]);
      for(std::map<int, double>::iterator it=dphi.begin(); it!=dphi.end(); it++)
        {
          double deltaphi = TMath::Abs(it->second - dphi.at(jleading));
          double filldeltaphi = deltaphi<M_PI?deltaphi:(2*M_PI-deltaphi);              
          if(it->first==jleading || dcand.Dpt[it->first]==ptleading) continue;
          int idphi = xjjc::findibin(&dphiBins, filldeltaphi);
          if(idphi<0) return 1;
          if(dcand.Dgen[it->first]==23333) 
            {
              hmassSignal[idphi]->Fill(dcand.Dmass[it->first]);
            }
          if(dcand.Dgen[it->first]==23344) hmassSwapped[idphi]->Fill(dcand.Dmass[it->first]);
        }
    }
  xjjc::progressbar_summary(nentries);
  
  TFile* outf = new TFile(Form("%s.root",outfname.Data()),"recreate");
  outf->cd();
  hmassSignalLD->Write();
  hmassSwappedLD->Write();
  for(int i=0;i<nDphiBins;i++) 
    {
      hmassSignal[i]->Write();
      hmassSwapped[i]->Write();
    }
  outf->Write();
  outf->Close();

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc==7)
    {
      stybkg_savefittpl(argv[1], argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]));
      return 0;
    }
  else
    {
      return 1;
    }
}
