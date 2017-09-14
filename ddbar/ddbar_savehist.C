#include "ddbar.h"

void ddbar_savehist(TString infname, TString outfname, 
                    TString collisionsyst, 
                    Float_t leading_ptmin, Float_t other_ptmin, Float_t leading_trkptmin)
{
  int arguerr(TString collisionsyst);
  if(arguerr(collisionsyst)) return;

  //
  if(initbinning()) return;
  if(initcutval(collisionsyst)) return;
  if(createhists("savehist")) return;

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
  for(int k=0;k<cutval_skim.size();k++) xjjroot::setbranchaddress(ntSkim, cutval_skim[k].Data(), &val_skim[k]);
  
  //  
  int nentries = ntDkpi->GetEntries();
  for(int i=0;i<nentries;i++)
    {
      if(i%10000==0) std::cout<<std::setiosflags(std::ios::left)<<"  [ \033[1;36m"<<std::setw(10)<<i<<"\033[0m"<<" / "<<std::setw(10)<<nentries<<" ] "<<"\033[1;36m"<<std::setw(4)<<Form("%.0f%s",100.*i/nentries,"%")<<"\033[0m"<<"   >>   ddbar_savehist("<<std::setw(20)<<Form("%s)",collisionsyst.Data())<<"\r"<<std::flush;

      //
      ntHlt->GetEntry(i);
      ntSkim->GetEntry(i);
      ntDkpi->GetEntry(i);
      ntGen->GetEntry(i);

      // add hlt sel...
      
      for(int k=0, npass=0;k<cutval_skim.size();k++)
        if(!val_skim[k]) {npass++; break;}
      if(npass) continue;
      if(TMath::Abs(dcand.PVz)>15) continue;
      
      // 
      int jleadingD = -1;
      double ptleadingD = 0;
      std::vector<std::map<int, double>> dphiD(nhistD);
      for(int j=0;j<dcand.Dsize;j++)
        {
          int ipt = xjjc::findibin(ptBins, (double)dcand.Dpt[j]);
          if(ipt<0) continue;
          if(initcutval_ptdep(collisionsyst, ipt)) return;
          dcand.settrkcut(leading_trkptmin, cutval_trkEta, cutval_trkPtErr);
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, leading_ptmin);
          if(dcand.isselected(j) && dcand.Dpt[j]>ptleadingD)
            {
              jleadingD = j;
              ptleadingD = dcand.Dpt[j];
            }          
          dcand.settrkcut(cutval_trkPt, cutval_trkEta, cutval_trkPtErr);
          dcand.setDcut(cutval_Dy, cutval_Dsvpv, cutval_Dalpha, cutval_Dchi2cl, other_ptmin);
          if(!dcand.isselected(j) || dcand.Dgen[j]!=23333) continue;
          for(int l=0;l<nhistD;l++) 
            dphiD[l].insert(std::pair<int, double>(j, dcand.Dphi[j]));
        }
      // 
      if(jleadingD>=0)
        {
          // if(!(dcand.Dmass[jleadingD] > 1.7 && dcand.Dmass[jleadingD] < 2.0)) continue;
          std::vector<bool> leadingsel = // [nhist]
            {
              dcand.Dgen[jleadingD]==23333,
              dcand.Dgen[jleadingD]==23333,
            };
          for(int l=0;l<nhistD;l++)
            {
              if(!leadingsel[l]) continue;
              if(dphiD[l].empty()) continue;
              for(std::map<int, double>::iterator it=dphiD[l].begin(); it!=dphiD[l].end(); it++)
                {
                  if(it->first==jleadingD || dcand.Dpt[it->first]==ptleadingD) continue; // skip leading D and swapped cand of leading D
                  if(l==0 && dcand.Dtype[it->first]!=dcand.Dtype[jleadingD]) continue;
                  if(l==1 && dcand.Dtype[it->first]==dcand.Dtype[jleadingD]) continue;
                  double deltaphi = it->second - dcand.Dphi[jleadingD];
                  double filldeltaphi = deltaphi<-M_PI/2.?(deltaphi+2*M_PI):(deltaphi>3*M_PI/2.?(deltaphi-2*M_PI):deltaphi);
                  ahdphiD[l]->Fill(filldeltaphi);
                }
            }
        }

      // 
      int jleadingG = -1;
      double ptleadingG = 0;
      std::vector<std::map<int, double>> dphiG(nhistG);
      for(int j=0;j<dcand.Gsize;j++)
        {
          if(TMath::Abs(dcand.Gy[j])<1. && (dcand.GisSignal[j]==1 || dcand.GisSignal[j]==2) && dcand.Gpt[j]>ptleadingG)
            {
              jleadingG = j;
              ptleadingG = dcand.Gpt[j];
            }
          if(TMath::Abs(dcand.Gy[j])>1. && (dcand.GisSignal[j]==1 || dcand.GisSignal[j]==2)) continue;
          for(int l=0;l<nhistG;l++) 
            dphiG[l].insert(std::pair<int, double>(j, dcand.Gphi[j]));
        }
      // 
      if(jleadingG>=0)
        {
          for(int l=0;l<nhistG;l++)
            {
              if(dphiG[l].empty()) continue;
              for(std::map<int, double>::iterator it=dphiG[l].begin(); it!=dphiG[l].end(); it++)
                {
                  if(it->first==jleadingG) continue; // skip leading D
                  if(l==0 && dcand.GisSignal[it->first]!=dcand.GisSignal[jleadingG]) continue;
                  if(l==1 && dcand.GisSignal[it->first]==dcand.GisSignal[jleadingG]) continue;
                  double deltaphi = it->second - dcand.Gphi[jleadingG];
                  double filldeltaphi = deltaphi<-M_PI/2.?(deltaphi+2*M_PI):(deltaphi>3*M_PI/2.?(deltaphi-2*M_PI):deltaphi);
                  ahdphiG[l]->Fill(filldeltaphi);
                }
            }
        }
    }
  std::cout<<std::setiosflags(std::ios::left)<<"  Processed "<<"\033[1;31m"<<nentries<<"\033[0m event(s)."<<"   >>   ddbar_savehist("<<std::setw(30)<<Form("%s)",collisionsyst.Data())<<std::endl;
  std::cout<<std::endl;
  
  for(int l=0;l<nhistD;l++)
    {
      ahdphiD[l]->Sumw2();
      ahdphiD[l]->Scale(1./ahdphiD[l]->Integral(), "width");
    }
  for(int l=0;l<nhistG;l++)
    {
      ahdphiG[l]->Sumw2();
      ahdphiG[l]->Scale(1./ahdphiG[l]->Integral(), "width");
    }

  TFile* outf = new TFile(Form("%s.root",outfname.Data()),"recreate");
  outf->cd();
  if(writehists("savehist")) return;
  outf->Write();
  outf->Close();
}

int main(int argc, char* argv[])
{
  if(argc==7)
    {
      ddbar_savehist(argv[1], argv[2], argv[3], atof(argv[4]), atof(argv[5]), atof(argv[6]));
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
