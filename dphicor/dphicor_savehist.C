#include "dphicor.h"

int dphicor_savehist(TString infname, TString outfname, TString collisionsyst, Int_t isMC, Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
{
  std::cout<<std::endl;
  initcutval(collisionsyst);

  for(int i=0;i<=nDphiBins;i++) dphiBins[i] = minDphi+i*(maxDphi-minDphi)/nDphiBins;

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
  xjjrootuti::setbranchaddress(ntHlt, cutval_hlt.Data(), &val_hlt);
  ntSkim->SetBranchStatus("*", 0);
  int* val_skim = new int[cutval_skim.size()];
  for(int k=0;k<cutval_skim.size();k++) xjjrootuti::setbranchaddress(ntSkim, cutval_skim[k].Data(), &val_skim[k]);

  TH1D** hdphi = new TH1D*[nhist];
  TH1D** hmass = new TH1D*[nDphiBins];
  for(int l=0;l<nhist;l++) hdphi[l] = new TH1D(histname[l], ";#Delta#phi;Entries", nDphiBins_fine, 0, M_PI);
  for(int i=0;i<nDphiBins;i++) hmass[i] = new TH1D(Form("hmass_%d",i), ";m_{#piK} (GeV/c^{2});Entries / (5 MeV/c^{2})", 60, 1.7, 2.0);

  int nentries = ntDkpi->GetEntries();
  for(int i=0;i<nentries;i++)
    {
      if(i%10000==0) xjjuti::progressbar(i, nentries);
      ntDkpi->GetEntry(i);
      ntHlt->GetEntry(i);
      ntSkim->GetEntry(i);

      // if(!val_hlt && !isMC) continue;
      for(int k=0;k<cutval_skim.size();k++)
        {
          if(!val_skim[k]) continue;
        }
      if(TMath::Abs(dcand.PVz)>15) continue;
      
      int jleading[nhist];
      float ptleading[nhist], phileading[nhist];
      xjjuti::initarray(&jleading, -1);
      xjjuti::initarray(&ptleading, (float)0);
      xjjuti::initarray(&phileading, (float)0);
      std::map<int, float> dphi[nhist];
      for(int j=0;j<dcand.Dsize;j++)
        {
          int ipt = xjjuti::findibin(&ptBins, dcand.Dpt[j]);
          if(ipt<0) continue;
          int err_initcutval_ptdep = initcutval_ptdep(collisionsyst, ipt);
          if(err_initcutval_ptdep) return 1;
          dcand.settrkcut(leading_trkptmin, cutval_trkEta, cutval_trkPtErr);
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, leading_ptmin);
          if(dcand.isselected(j))
            {
              if(dcand.Dpt[j]>ptleading[0])
                {
                  jleading[0] = j;
                  ptleading[0] = dcand.Dpt[j];
                  phileading[0] = dcand.Dphi[j];
                  jleading[2] = j;
                  ptleading[2] = dcand.Dpt[j];
                  phileading[2] = dcand.Dphi[j];
                }
              if(dcand.Dgen[j]==23333 && dcand.Dpt[j]>ptleading[1])
                {
                  jleading[1] = j;
                  ptleading[1] = dcand.Dpt[j];              
                  phileading[1] = dcand.Dphi[j];
                  jleading[3] = j;
                  ptleading[3] = dcand.Dpt[j];              
                  phileading[3] = dcand.Dphi[j];
                }
            }
          dcand.settrkcut(cutval_trkPt, cutval_trkEta, cutval_trkPtErr);
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, other_ptmin);
          if(!dcand.isselected(j)) continue;
          dphi[0].insert(std::pair<int, float>(j, dcand.Dphi[j]));
          dphi[3].insert(std::pair<int, float>(j, dcand.Dphi[j]));
          if(dcand.Dgen[j]==23333) 
            {
              dphi[1].insert(std::pair<int, float>(j, dcand.Dphi[j]));
              dphi[2].insert(std::pair<int, float>(j, dcand.Dphi[j]));
            }
        }
      for(int l=0;l<nhist;l++)
        {
          if(jleading[l]<0 || dphi[l].size()==0) continue;
          for(std::map<int, float>::iterator it=dphi[l].begin(); it!=dphi[l].end(); it++)
            {
              if(it->first==jleading[l]) continue;
              float deltaphi = TMath::Abs(it->second - phileading[l]);
              float filldeltaphi = deltaphi<M_PI?deltaphi:(2*M_PI-deltaphi);              
              hdphi[l]->Fill(filldeltaphi);
              if(l) continue;
              int idphi = xjjuti::findibin(&dphiBins, filldeltaphi);
              if(idphi<0) return 1;
              hmass[idphi]->Fill(dcand.Dmass[it->first]);
            }
        }
    }
  xjjuti::progressbar_summary(nentries);
  
  TFile* outf = new TFile(Form("%s.root",outfname.Data()),"recreate");
  outf->cd();
  for(int l=0;l<nhist;l++) hdphi[l]->Write();
  for(int i=0;i<nDphiBins;i++) hmass[i]->Write();
  outf->Write();
  outf->Close();

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc==8)
    {
      dphicor_savehist(argv[1], argv[2], argv[3], atoi(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));
      return 0;
    }
  else
    {
      return 1;
    }
}