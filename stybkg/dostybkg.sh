#!/bin/bash
# dostybkg.sh #

DO_STYBKG_SAVETPL=${1:-0}
DO_STYBKG_SAVEHIST=${2:-0}
DO_STYBKG_USEHIST=${3:-0}
DO_STYBKG_PLOTHIST=${3:-0}

# Select the systems the macros run on
jCOLSYST=(0 1)
lLEAD=(0)
kOTHER=(2)

##

# nCOL loop
COLSYST=('pp' 'PbPb')

# nLEAD loop
LEADING_TRKPTMIN=(1)
LEADING_PTMIN=(20)

# nOTHER loop
OTHER_PTMIN=(0 1 2)

# dataset[nCOL]
INPUTDNAME=(
    # "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20170427_DfinderMC_pp_20170427_D0_dPt0tkPt0p1_Pythia8_prompt_D0pt0p0_pp502_TuneCUETP8M1_pthatweight.root"
    "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20160513_DfinderMC_pp_20160502_dPt0tkPt0p5_D0Dstar_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
    "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
)
INPUTSNAME=(
    # "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20170427_DfinderMC_pp_20170427_D0_dPt0tkPt0p1_Pythia8_prompt_D0pt0p0_pp502_TuneCUETP8M1_pthatweight.root"
    "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20160513_DfinderMC_pp_20160502_dPt0tkPt0p5_D0Dstar_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
    "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
)

# Do not touch the macros below if you don't know what they mean #

[[ $DO_STYBKG_SAVETPL -eq 0 && $DO_STYBKG_SAVEHIST -eq 0 && $DO_STYBKG_USEHIST -eq 0 && $DO_STYBKG_PLOTHIST -eq 0 ]] && echo "./dostybkg.sh [DO_STYBKG_SAVETPL] [DO_STYBKG_SAVEHIST] [DO_STYBKG_USEHIST] [DO_STYBKG_PLOTHIST]"

#
nCOL=${#COLSYST[@]}
nLEAD=${#LEADING_TRKPTMIN[@]}
nOTHER=${#OTHER_PTMIN[@]}

#
NC='\033[0m'
FUNCOLOR='\033[1;33m'
ARGCOLOR='\033[1;32m'
ERRCOLOR='\033[1;31m'

#
FOLDERS=("rootfiles" "plots" "plotfits")
for i in ${FOLDERS[@]}
do
    if [[ ! -d $i ]]
    then
        mkdir -p $i
    fi
done

#
function float_to_string()
{
    if [[ $# -ne 1 ]]
    then
        echo -e "${ERRCOLOR}error:${NC} invalid argument number - float_to_string()"
        return 1
    fi
    part1=`echo $1 | awk -F "." '{print $1}'`
    part2=`echo $1 | awk -F "." '{print $2}'`
    rt_float_to_string=${part1:-0}p${part2:-0}
    echo $rt_float_to_string
}

##

# stybkg_savetpl.C + stybkg_savehist.C#
g++ stybkg_savetpl.C $(root-config --cflags --libs) -g -o stybkg_savetpl.exe || return 1
g++ stybkg_savehist.C $(root-config --cflags --libs) -g -o stybkg_savehist.exe || return 1
for j in ${jCOLSYST[@]}
do
    for l in ${lLEAD[@]}
    do
        for k in ${kOTHER[@]}
        do
            tPOSTFIX=STYBKG_${COLSYST[j]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[l]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[k]})_leadingtrkptmin_$(float_to_string ${LEADING_TRKPTMIN[l]})
            if [[ $DO_STYBKG_SAVETPL -eq 1 ]]
            then
                echo -e "-- Processing ${FUNCOLOR}stybkg_savetpl.C${NC} :: ${ARGCOLOR}${COLSYST[j]}${NC}"
                ./stybkg_savetpl.exe ${INPUTSNAME[j]} rootfiles/ftpl_${tPOSTFIX} ${COLSYST[j]} ${LEADING_PTMIN[l]} ${OTHER_PTMIN[k]} ${LEADING_TRKPTMIN[l]} &
                echo
            fi
            if [[ $DO_STYBKG_SAVEHIST -eq 1 ]]
            then
                echo -e "-- Processing ${FUNCOLOR}stybkg_savehist.C${NC} :: ${ARGCOLOR}${COLSYST[j]}${NC}"
                ./stybkg_savehist.exe ${INPUTDNAME[j]} rootfiles/fhist_${tPOSTFIX} ${COLSYST[j]} ${LEADING_PTMIN[l]} ${OTHER_PTMIN[k]} ${LEADING_TRKPTMIN[l]} &
                echo
            fi
        done
    done
done
wait
rm stybkg_savehist.exe
rm stybkg_savetpl.exe

# stybkg_usehist.C #
g++ stybkg_usehist.C $(root-config --cflags --libs) -g -o stybkg_usehist.exe || return 1
if [[ $DO_STYBKG_USEHIST -eq 1 ]]
then
    for j in ${jCOLSYST[@]}
    do
        for l in ${lLEAD[@]}
        do
            for k in ${kOTHER[@]}
            do
                tPOSTFIX=STYBKG_${COLSYST[j]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[l]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[k]})_leadingtrkptmin_$(float_to_string ${LEADING_TRKPTMIN[l]})
                echo -e "-- Processing ${FUNCOLOR}stybkg_usehist.C${NC} :: ${ARGCOLOR}${COLSYST[j]}${NC}"
                ./stybkg_usehist.exe rootfiles/fhist_${tPOSTFIX} rootfiles/ftpl_${tPOSTFIX} rootfiles/fdphi_${tPOSTFIX} ${tPOSTFIX} ${COLSYST[j]} ${LEADING_PTMIN[l]} ${OTHER_PTMIN[k]} ${LEADING_TRKPTMIN[l]}
                echo
            done
        done
    done
fi
rm stybkg_usehist.exe

# stybkg_plothist.C #
g++ stybkg_plothist.C $(root-config --cflags --libs) -g -o stybkg_plothist.exe || return 1
if [[ $DO_STYBKG_PLOTHIST -eq 1 ]]
then
    for j in ${jCOLSYST[@]}
    do
        for l in ${lLEAD[@]}
        do
            for k in ${kOTHER[@]}
            do
                tPOSTFIX=STYBKG_${COLSYST[j]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[l]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[k]})_leadingtrkptmin_$(float_to_string ${LEADING_TRKPTMIN[l]})
                echo -e "-- Processing ${FUNCOLOR}stybkg_plothist.C${NC} :: ${ARGCOLOR}${COLSYST[j]}${NC}"
                ./stybkg_plothist.exe rootfiles/fdphi_${tPOSTFIX} ${tPOSTFIX} ${COLSYST[j]} ${LEADING_PTMIN[l]} ${OTHER_PTMIN[k]} ${LEADING_TRKPTMIN[l]}
                echo
            done
        done
    done
fi
rm stybkg_plothist.exe

