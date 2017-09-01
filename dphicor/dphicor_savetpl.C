#include "dphicor.h"

void dphicor_savetpl(TString infname, TString outfname, 
                     TString collisionsyst, Int_t isMC,
                     Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
{
  int arguerr(TString collisionsyst, Int_t isMC);
  if(arguerr(collisionsyst, isMC)) return;

  // 
  if(initbinning()) return;
  if(initcutval(collisionsyst)) return;
  if(createhists("savetpl")) return;

  TFile* inf = new TFile(infname);
  TTree* ntDkpi = (TTree*)inf->Get("ntDkpi");
  TTree* ntGen = (TTree*)inf->Get("ntGen");
  TTree* ntHlt = (TTree*)inf->Get("ntHlt");
  TTree* ntHi = (TTree*)inf->Get("ntHi");
  TTree* ntSkim = (TTree*)inf->Get("ntSkim");

  readD dcand(MAX_XB, MAX_GEN, true);
  dcand.setbranchesaddress(ntDkpi, ntGen, ntHi);
  ntHlt->SetBranchStatus("*", 0);
  int* val_hlt = new int[cutval_hlt.size()];
  for(int k=0;k<cutval_hlt.size();k++) xjjroot::setbranchaddress(ntHlt, cutval_hlt[k].Data(), &(val_hlt[k]));
  ntSkim->SetBranchStatus("*", 0);
  int* val_skim = new int[cutval_skim.size()];
  for(int k=0;k<cutval_skim.size();k++) xjjroot::setbranchaddress(ntSkim, cutval_skim[k].Data(), &(val_skim[k]));

  // 
  int nentries = ntDkpi->GetEntries();
  for(int i=0;i<nentries;i++)
    {
      if(i%10000==0) std::cout<<std::setiosflags(std::ios::left)<<"  [ \033[1;36m"<<std::setw(10)<<i<<"\033[0m"<<" / "<<std::setw(10)<<nentries<<" ] "<<"\033[1;36m"<<std::setw(4)<<Form("%.0f%s",100.*i/nentries,"%")<<"\033[0m"<<"   >>   dphicor_savetpl("<<std::setw(5)<<Form("%s,",collisionsyst.Data())<<" "<<std::setw(20)<<Form("%s)",tMC[isMC].Data())<<"\r"<<std::flush;

      //
      ntDkpi->GetEntry(i);
      ntHlt->GetEntry(i);
      ntSkim->GetEntry(i);

      // add hlt sel...

      for(int k=0, npass=0;k<cutval_skim.size();k++)
        if(!val_skim[k]) {npass++; break;}
      if(npass) continue;
      if(TMath::Abs(dcand.PVz)>15) continue;

      // 
      int jleading = -1;
      double ptleading = 0;
      std::map<int, double> dphi;
      for(int j=0;j<dcand.Dsize;j++)
        {
          int ipt = xjjc::findibin(ptBins, (double)dcand.Dpt[j]);
          if(ipt<0) continue;
          if(initcutval_ptdep(collisionsyst, ipt)) return;
          dcand.settrkcut(leading_trkptmin, cutval_trkEta, cutval_trkPtErr);
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, leading_ptmin);
          if(dcand.isselected(j) && dcand.Dpt[j]>ptleading/* && dcand.Dgen[j]==23333*/)
            {
              jleading = j;
              ptleading = dcand.Dpt[j];              
            }
          dcand.settrkcut(cutval_trkPt, cutval_trkEta, cutval_trkPtErr);
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, other_ptmin);
          if(!dcand.isselected(j)) continue;
          dphi.insert(std::pair<int, double>(j, dcand.Dphi[j]));
        }

      // 
      if(jleading<0) continue;
      if(!(dcand.Dmass[jleading] > 1.7 && dcand.Dmass[jleading] < 2.0)) continue;
      if(dcand.Dgen[jleading]==23333) hmassSignalLD->Fill(dcand.Dmass[jleading]);
      if(dcand.Dgen[jleading]==23344) hmassSwappedLD->Fill(dcand.Dmass[jleading]);
      for(std::map<int, double>::iterator it=dphi.begin(); it!=dphi.end(); it++)
        {
          if(it->first==jleading || dcand.Dpt[it->first]==ptleading) continue; // skip leading D and swapped cand of leading D
          double deltaphi = TMath::Abs(it->second - dphi.at(jleading));
          double filldeltaphi = deltaphi<M_PI?deltaphi:(2*M_PI-deltaphi);              
          int idphi = xjjc::findibin(dphiBins, filldeltaphi);
          if(idphi<0) return;
          if(dcand.Dgen[it->first]==23333) ahmassSignal[idphi]->Fill(dcand.Dmass[it->first]);
          if(dcand.Dgen[it->first]==23344) ahmassSwapped[idphi]->Fill(dcand.Dmass[it->first]);
        }
    }
  std::cout<<std::setiosflags(std::ios::left)<<"  Processed "<<"\033[1;31m"<<nentries<<"\033[0m event(s)."<<"   >>   dphicor_savetpl("<<std::setw(5)<<Form("%s,",collisionsyst.Data())<<" "<<std::setw(30)<<Form("%s)",tMC[isMC].Data())<<std::endl;
  std::cout<<std::endl;

  // 
  TFile* outf = new TFile(Form("%s.root",outfname.Data()),"recreate");
  outf->cd();
  if(writehists("savetpl")) return;
  outf->Write();
  outf->Close();

}

int main(int argc, char* argv[])
{
  if(argc==8)
    {
      dphicor_savetpl(argv[1], argv[2], argv[3], atoi(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));
      return 0;
    }
  else
    {
      std::cout<<"  Error: invalid arguments number - dphicor_savetpl()"<<std::endl;
      return 1;
    }
}

int arguerr(TString collisionsyst, Int_t isMC)
{
  if(collsyst_list.find(collisionsyst)==collsyst_list.end())
    {
      std::cout<<"\033[1;31merror:\033[0m invalid \"collisionsyst\""<<std::endl;
      return 1;
    }
  if(isMC!=0 && isMC!=1)
    {
      std::cout<<"\033[1;31merror:\033[0m invalid \"isMC\""<<std::endl;
      return 1;
    }
  return 0;
}
