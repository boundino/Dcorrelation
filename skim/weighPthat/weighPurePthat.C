using namespace std;
#include "uti.h"
#define MAX_GEN      6000

int weighPurePthat(TString ifname = "",
                   TString ofname = "")
{
  Bool_t isInsidebin(Float_t xpthat, Int_t i);
  cout<<endl;
  cout<<" -- Checking if input and output files are same"<<endl;
  if(ifname==ofname)
    {
      cout<<"    Error: Input file will be overwritten."<<endl;
      return 0;
    }
  cout<<" -- Opening unweighed sample"<<endl;
  TFile* inf = TFile::Open(ifname);
  TTree* dtrk = (TTree*)inf->Get("dtrk");
  Int_t Gsize; dtrk->SetBranchAddress("Gsize",&Gsize);
  std::vector<float>*     Gpt = 0;
  dtrk->SetBranchAddress("Gpt",&Gpt);
  std::vector<int>*       GisSignal = 0;
  dtrk->SetBranchAddress("GisSignal",&GisSignal);
  Float_t pthat; dtrk->SetBranchAddress("pthat",&pthat);

  Float_t weight[nBins],nweight[nBins];
  for(Int_t j=0;j<nBins;j++)
    {
      weight[j]=0;
      nweight[j]=0;
    }
  cout<<" -- Checking event number"<<endl;
  Int_t nentries = dtrk->GetEntries();
  cout<<" -- Calculating weights"<<endl;
  for(Int_t i=0;i<nentries;i++)
    {
      dtrk->GetEntry(i);
      if(i%100000==0) cout<<"    Processing event "<<setiosflags(ios::left)<<setw(7)<<i<<" / "<<nentries<<endl;
      for(Int_t j=0;j<nBins;j++)
        {
          if(isInsidebin(pthat,j)) nweight[j]++;
        }
    }
  cout<<" -- Weight results"<<endl;
  for(Int_t j=0;j<nBins;j++)
    {
      if(nweight[j]==0)
        {
          cout<<"    Error: Weight fails."<<endl;
          return 0;
        }
      weight[j] = (crosssec[j]-crosssec[j+1])/nweight[j];
      cout<<"    Pthat"<<setiosflags(ios::left)<<setw(3)<<pthatBin[j]<<": "<<weight[j]<<endl;
    }

  cout<<" -- Building weight branch"<<endl;
  TFile* otf = TFile::Open(ofname,"update");
  TTree* dtrknew = (TTree*)otf->Get("dtrk");
  Float_t pthatweight,maxDgenpt;
  TBranch* newBr_pthatweight = dtrknew->Branch("pthatweight", &pthatweight, "pthatweight/F");
  TBranch* newBr_maxDgenpt = dtrknew->Branch("maxDgenpt", &maxDgenpt, "maxDgenpt/F");
  cout<<" -- Filling weight branch"<<endl;
  for(Int_t i=0;i<nentries;i++)
    {
      dtrk->GetEntry(i);
      if(i%100000==0) cout<<"    Processing event "<<setiosflags(ios::left)<<setw(7)<<i<<" / "<<nentries<<endl;
      pthatweight=0;
      Float_t maxpt=0;
      for(Int_t k=0;k<Gsize;k++)
        {
          if(((*GisSignal)[k]==genSignal[0]||(*GisSignal)[k]==genSignal[1])&&(*Gpt)[k]>maxpt) maxpt=(*Gpt)[k];
        }
      maxDgenpt = maxpt;
      for(Int_t j=0;j<nBins;j++)
        {
          if(isInsidebin(pthat,j))
            {
              pthatweight = weight[j];
            }
        }
      newBr_pthatweight->Fill();
      newBr_maxDgenpt->Fill();
    }
  dtrknew->Write("", TObject::kOverwrite);

  cout<<" -- End"<<endl;
  cout<<endl;
  return 1;
}

Bool_t isInsidebin(Float_t xpthat, Int_t i)
{
  if(i>=nBins)
    {
      cout<<"    Error: invalid input"<<endl;
      return false;
    }
  if(i<(nBins-1)&&xpthat>pthatBin[i]&&xpthat<pthatBin[i+1]) return true;
  else if(i==(nBins-1)&&xpthat>pthatBin[i]) return true;
  else return false;
}

int main(int argc, char *argv[])
{
  if(argc==3)
    {
      weighPurePthat(argv[1], argv[2]);
      return 1;
    }
  else
    {
      std::cout<<"Invalid parameter number"<<std::endl;
      return 0;
    }
}
