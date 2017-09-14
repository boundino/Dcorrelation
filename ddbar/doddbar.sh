#!/bin/bash
# doddbar.sh #

DO_DDBAR_SAVETPL=${1:-0}
DO_DDBAR_SAVEHIST=${2:-0}
DO_DDBAR_USEHIST=${3:-0}
DO_DDBAR_PLOTHIST=${3:-0}

# Select the systems the macros run on
jCOLSYST=(0 1)
lLEAD=(0)
kOTHER=(0)

##

# nCOL loop
COLSYST=('pp' 'PbPb')

# nLEAD loop
LEADING_TRKPTMIN=(0.7)
LEADING_PTMIN=(20)

# nOTHER loop
OTHER_PTMIN=(2)

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

[[ $DO_DDBAR_SAVETPL -eq 0 && $DO_DDBAR_SAVEHIST -eq 0 && $DO_DDBAR_USEHIST -eq 0 && $DO_DDBAR_PLOTHIST -eq 0 ]] && echo "./doddbar.sh [DO_DDBAR_SAVETPL] [DO_DDBAR_SAVEHIST] [DO_DDBAR_USEHIST] [DO_DDBAR_PLOTHIST]"

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

# ddbar_savetpl.C + ddbar_savehist.C#
g++ ddbar_savehist.C $(root-config --cflags --libs) -g -o ddbar_savehist.exe || return 1
for j in ${jCOLSYST[@]}
do
    for l in ${lLEAD[@]}
    do
        for k in ${kOTHER[@]}
        do
            tPOSTFIX=DDBAR_${COLSYST[j]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[l]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[k]})_leadingtrkptmin_$(float_to_string ${LEADING_TRKPTMIN[l]})
            if [[ $DO_DDBAR_SAVEHIST -eq 1 ]]
            then
                echo -e "-- Processing ${FUNCOLOR}ddbar_savehist.C${NC} :: ${ARGCOLOR}${COLSYST[j]}${NC}"
                ./ddbar_savehist.exe ${INPUTDNAME[j]} rootfiles/fhist_${tPOSTFIX} ${COLSYST[j]} ${LEADING_PTMIN[l]} ${OTHER_PTMIN[k]} ${LEADING_TRKPTMIN[l]} &
                echo
            fi
        done
    done
done
wait
rm ddbar_savehist.exe

# ddbar_usehist.C #
g++ ddbar_usehist.C $(root-config --cflags --libs) -g -o ddbar_usehist.exe || return 1
if [[ $DO_DDBAR_USEHIST -eq 1 ]]
then
    for j in ${jCOLSYST[@]}
    do
        for l in ${lLEAD[@]}
        do
            for k in ${kOTHER[@]}
            do
                tPOSTFIX=DDBAR_${COLSYST[j]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[l]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[k]})_leadingtrkptmin_$(float_to_string ${LEADING_TRKPTMIN[l]})
                echo -e "-- Processing ${FUNCOLOR}ddbar_usehist.C${NC} :: ${ARGCOLOR}${COLSYST[j]}${NC}"
                ./ddbar_usehist.exe rootfiles/fhist_${tPOSTFIX} rootfiles/ftpl_${tPOSTFIX} rootfiles/fdphi_${tPOSTFIX} ${tPOSTFIX} ${COLSYST[j]} ${LEADING_PTMIN[l]} ${OTHER_PTMIN[k]} ${LEADING_TRKPTMIN[l]}
                echo
            done
        done
    done
fi
rm ddbar_usehist.exe

# ddbar_plothist.C #
g++ ddbar_plothist.C $(root-config --cflags --libs) -g -o ddbar_plothist.exe || return 1
if [[ $DO_DDBAR_PLOTHIST -eq 1 ]]
then
    for j in ${jCOLSYST[@]}
    do
        for l in ${lLEAD[@]}
        do
            for k in ${kOTHER[@]}
            do
                tPOSTFIX=DDBAR_${COLSYST[j]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[l]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[k]})_leadingtrkptmin_$(float_to_string ${LEADING_TRKPTMIN[l]})
                echo -e "-- Processing ${FUNCOLOR}ddbar_plothist.C${NC} :: ${ARGCOLOR}${COLSYST[j]}${NC}"
                ./ddbar_plothist.exe rootfiles/fdphi_${tPOSTFIX} ${tPOSTFIX} ${COLSYST[j]} ${LEADING_PTMIN[l]} ${OTHER_PTMIN[k]} ${LEADING_TRKPTMIN[l]}
                echo
            done
        done
    done
fi
rm ddbar_plothist.exe

