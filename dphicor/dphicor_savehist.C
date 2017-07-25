#include "dphicor.h"

int dphicor_savehist(TString infname, TString outfname, TString collisionsyst, Int_t isMC, Float_t leading_ptmin, Float_t other_ptmin)
{
  std::cout<<std::endl;
  initcutval(collisionsyst);

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

  int nentries = ntDkpi->GetEntries();
  TH1D** hdphi = new TH1D*[nhist];
  for(int l=0;l<nhist;l++) hdphi[l] = new TH1D(histname[l], ";#Delta#phi;Entries", 50, 0, 3.2);
  for(int i=0;i<nentries;i++)
    {
      if(i%10000==0) xjjuti::progressbar(i, nentries);
      ntDkpi->GetEntry(i);
      ntHlt->GetEntry(i);
      ntSkim->GetEntry(i);

      if(!val_hlt && !isMC) continue;
      for(int k=0;k<cutval_skim.size();k++)
        {
          if(!val_skim[k]) continue;
        }
      if(TMath::Abs(dcand.PVz)>15) continue;
      
      int jleading[nhist];
      float ptleading[nhist];
      for(int l=0;l<nhist;l++)
        {
          jleading[l] = -1;
          ptleading[l] = 0;
        }
      std::map<int, float> dphi[nhist];
      for(int j=0;j<dcand.Dsize;j++)
        {
          int ipt = xjjuti::findibin(&ptBins, dcand.Dpt[j]);
          if(ipt<0) continue;
          int err_initcutval_ptdep = initcutval_ptdep(collisionsyst, ipt);
          if(err_initcutval_ptdep) return 1;
          dcand.settrkcut(cutval_trkPt, cutval_trkEta, cutval_trkPtErr);
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, leading_ptmin);
          if(dcand.isselected(j))
            {
              if(dcand.Dpt[j]>ptleading[0])
                {
                  jleading[0] = j;
                  ptleading[0] = dcand.Dpt[j];
                }
              if(dcand.Dgen[j]==23333 && dcand.Dpt[j]>ptleading[1])
                {
                  jleading[1] = j;
                  ptleading[1] = dcand.Dpt[j];              
                }
            }
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, other_ptmin);
          if(!dcand.isselected(j)) continue;
          dphi[0].insert(std::pair<int, float>(j, dcand.Dphi[j]));
          if(dcand.Dgen[j]==23333) dphi[1].insert(std::pair<int, float>(j, dcand.Dphi[j]));
        }
      for(int l=0;l<nhist;l++)
        {
          if(jleading[l]<0) continue;
          for(std::map<int, float>::iterator it=dphi[l].begin(); it!=dphi[l].end(); it++)
            {
              if(it->first!=jleading[l]) hdphi[l]->Fill(TMath::Abs(it->second - dphi[l].at(jleading[l])));
            }
        }
    }
  std::cout<<std::endl;
  xjjuti::progressbar_summary(nentries);
  
  TFile* outf = new TFile(Form("%s.root",outfname.Data()),"recreate");
  outf->cd();
  for(int l=0;l<nhist;l++) hdphi[l]->Write();
  outf->Write();
  outf->Close();

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc==7)
    {
      dphicor_savehist(argv[1], argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]));
      return 0;
    }
  else
    {
      return 1;
    }
}
