#!/bin/bash
# dophicor.sh #

# -1: loop all bins
jCOLSYST=-1
#
DO_DPHICOR_SAVEHIST=1
DO_DPHICOR_USEHIST=1

# nCOL loop
COLSYST=('pp' 'pp' 'PbPb' 'PbPb')
ISMC=(1 0 1 0)
LEADING_PTMIN=(20 20 20 20)
OTHER_PTMIN=(2 2 2 2)

INPUTNAME=("/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20160513_DfinderMC_pp_20160502_dPt0tkPt0p5_D0Dstar_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
    "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20160330_HeavyFlavor_DfinderData_pp_20160329_dPt0tkPt1_D0Dstar3p5p_goldenjson.root"
    "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
    "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20160405_HIHardProbes_DfinderData_PbPb_20160402_dPt0tkPt2p5_D0Dstar3p5p_FINALJSON.root")

# Do not touch the macros below if you don't know what they mean #
##
nCOL=${#COLSYST[@]}
TMC=("data" "MC")

#
NC='\033[0m'

#
function float_to_string()
{
    if [[ $# -ne 1 ]]
    then
        echo -e "\033[1;31merror:${NC} invalid argument number - float_to_string()"
        exit 1
    fi
    part1=`echo $1 | awk -F "." '{print $1}'`
    part2=`echo $1 | awk -F "." '{print $2}'`
    rt_float_to_string=${part1:-0}p${part2:-0}
    echo $rt_float_to_string
}

function run_this_bin()
{
    if [[ $# -ne 2 ]]
    then
        echo -e "\033[1;31merror:${NC} invalid argument number - run_this_bin()"
        exit 1
    fi
    rt_run_this_bin=0
    if [ $1 -eq $2 ] || [ $2 -eq -1 ]
    then
        rt_run_this_bin=1
    fi
    echo $rt_run_this_bin
}

#
FOLDERS=("rootfiles" "plots")
for i in ${FOLDERS[@]}
do
    if [ ! -d $i ]
    then
        mkdir -p $i
    fi
done

##

# dphicor_savehist.C #
if [ $DO_DPHICOR_SAVEHIST -eq 1 ]
then
    g++ dphicor_savehist.C $(root-config --cflags --libs) -g -o dphicor_savehist.exe
    j=0
    while ((j<$nCOL))
    do
        if [ $(run_this_bin $j $jCOLSYST) -eq 1 ]
        then
            TEND=DCOR_${COLSYST[j]}_${TMC[${ISMC[j]}]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[j]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[j]})
            echo -e "-- Processing \033[1;33mdphicor_savehist.C${NC}, \033[1;32m${COLSYST[j]}${NC}, \033[1;32m${TMC[${ISMC[j]}]}${MC} leading pT \033[1;32m> ${LEADING_PTMIN[j]} GeV${NC}, other pT \033[1;32m> ${OTHER_PTMIN[j]}${NC}"
            set -x
            ./dphicor_savehist.exe "${INPUTNAME[j]}" "rootfiles/fdphi_${TEND}" "${COLSYST[j]}" ${ISMC[j]} ${LEADING_PTMIN[j]} ${OTHER_PTMIN[j]}
            set +x
            echo
        fi
        j=$(($j+1))
    done
    rm dphicor_savehist.exe
fi

# dphicor_usehist.C #
if [ $DO_DPHICOR_USEHIST -eq 1 ]
then
    g++ dphicor_usehist.C $(root-config --cflags --libs) -g -o dphicor_usehist.exe
    j=0
    while ((j<$nCOL))
    do
        if [ $(run_this_bin $j $jCOLSYST) -eq 1 ]
        then
            TEND=DCOR_${COLSYST[j]}_${TMC[${ISMC[j]}]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[j]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[j]})
            echo -e "-- Processing \033[1;33mdphicor_usehist.C${NC}, \033[1;32m${COLSYST[j]}${NC}, leading pT \033[1;32m> ${LEADING_PTMIN[j]} GeV${NC}, other pT \033[1;32m> ${OTHER_PTMIN[j]}${NC}"
            set -x
            ./dphicor_usehist.exe "rootfiles/fdphi_${TEND}" "${TEND}" "${COLSYST[j]}" ${ISMC[j]} ${LEADING_PTMIN[j]} ${OTHER_PTMIN[j]}
            set +x
            echo
        fi
        j=$(($j+1))
    done
    rm dphicor_usehist.exe
fi


