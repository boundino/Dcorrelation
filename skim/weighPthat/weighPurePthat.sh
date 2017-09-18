#!/bin/bash

###
INPUTFILE="/export/d00/scratch/jwang/Dtrk/MC/DtrkFiles_20170918_pp_5TeV_TuneCUETP8M1_Dfinder_MC_20170906_noweight.root"
OUTPUTFILE="/export/d00/scratch/jwang/Dtrk/MC/DtrkFiles_20170918_pp_5TeV_TuneCUETP8M1_Dfinder_MC_20170906_pthatweight.root"

## pp

# Dzero prompt
sed '1iconst int nBins=9; float pthatBin[nBins]={0,5,10,15,30,50,80,120,170}; float crosssec[nBins+1]={6.885e+09,1.516e+09,2.113e+08,5.524e+07,4.411e+06,5.986e+05,8.142e+04,1.329e+04,2.404e+03,0.}; int genSignal[2]={1,2};' weighPurePthat.C > weighPurePthat_tmp.C

# Dzero nonprompt
#sed '1iconst int nBins=9; float pthatBin[nBins]={0,5,10,15,30,50,80,120,170}; float crosssec[nBins+1]={1.895e+08,1.393e+08,3.798e+07,1.206e+07,1.231e+06,1.853e+05,2.725e+04,4.566e+03,8.569e+02,0.}; int genSignal[2]={1,2};' weighPurePthat.C > weighPurePthat_tmp.C

## PbPb

# Dzero prompt
#sed '1iconst int nBins=9; float pthatBin[nBins]={0,5,10,15,30,50,80,120,170}; float crosssec[nBins+1]={6.885e+09,1.516e+09,2.113e+08,5.524e+07,4.411e+06,5.986e+05,8.142e+04,1.329e+04,2.404e+03,0.}; int genSignal[2]={1,2};' weighPurePthat.C > weighPurePthat_tmp.C

# Dzero nonprompt
#sed '1iconst int nBins=9; float pthatBin[nBins]={0,5,10,15,30,50,80,120,170}; float crosssec[nBins+1]={1.895e+08,1.393e+08,3.798e+07,1.206e+07,1.231e+06,1.853e+05,2.725e+04,4.566e+03,8.569e+02,0.}; int genSignal[2]={1,2};' weighPurePthat.C > weighPurePthat_tmp.C

###

cp "$INPUTFILE" "$OUTPUTFILE"
g++ weighPurePthat_tmp.C $(root-config --cflags --libs) -g -o weighPurePthat_tmp.exe 
./weighPurePthat_tmp.exe "$INPUTFILE" "$OUTPUTFILE"
rm weighPurePthat_tmp.exe

rm weighPurePthat_tmp.C
