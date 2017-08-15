#!/bin/bash
# dostybkg.sh #

# -1: loop all bins
jCOLSYST=-1
lLEAD=-1
kOTHER=2
#
DO_STYBKG_SAVEFITTPL=${1:-0}
DO_STYBKG_SAVEHIST=${2:-0}
DO_STYBKG_USEHIST=${3:-0}

# nCOL loop
COLSYST=('pp' 'PbPb')
ISMC=(1 1)
NEEDFITTPL=(1 1)
# nLEAD loop
LEADING_TRKPTMIN=(1)
LEADING_PTMIN=(20)

# nOTHER loop
OTHER_PTMIN=(0 1 2)

INPUTDNAME=("/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20170427_DfinderMC_pp_20170427_D0_dPt0tkPt0p1_Pythia8_prompt_D0pt0p0_pp502_TuneCUETP8M1_pthatweight.root"
    # "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20160513_DfinderMC_pp_20160502_dPt0tkPt0p5_D0Dstar_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
    "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
)
INPUTSNAME=("/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20170427_DfinderMC_pp_20170427_D0_dPt0tkPt0p1_Pythia8_prompt_D0pt0p0_pp502_TuneCUETP8M1_pthatweight.root"
    "/export/d00/scratch/jwang/DntupleRunII/ntD_EvtBase_20160513_DfinderMC_PbPb_20160502_dPt1tkPt0p5_D0_prompt_Dpt2Dy1p1tkPt0p7tkEta2Decay2p9Dalpha0p14Skim_pthatweight.root"
)

# Do not touch the macros below if you don't know what they mean #
##
nCOL=${#COLSYST[@]}
nLEAD=${#LEADING_TRKPTMIN[@]}
nOTHER=${#OTHER_PTMIN[@]}
TMC=("data" "MC")

#
FOLDERS=("rootfiles" "plots" "plotfits")
for i in ${FOLDERS[@]}
do
    if [ ! -d $i ]
    then
        mkdir -p $i
    fi
done

#
NC='\033[0m'

#
function float_to_string()
{
    if [[ $# -ne 1 ]]
    then
        echo -e "\033[1;31merror:${NC} invalid argument number - float_to_string()"
        return 1
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
        return 1
    fi
    rt_run_this_bin=0
    if [ $1 -eq $2 ] || [ $2 -eq -1 ]
    then
        rt_run_this_bin=1
    fi
    echo $rt_run_this_bin
}

##

# stybkg_savefittpl.C #
if [ $DO_STYBKG_SAVEFITTPL -eq 1 ]
then
    g++ stybkg_savefittpl.C $(root-config --cflags --libs) -g -o stybkg_savefittpl.exe || return 1
    j=0
    while ((j<$nCOL))
    do
        if [ $(run_this_bin $j $jCOLSYST) -eq 1 ] && [ ${NEEDFITTPL[j]} -eq 1 ]
        then
            l=0
            while ((l<$nLEAD))
            do
                if [ $(run_this_bin $l $lLEAD) -eq 1 ]
                then
                    k=0
                    while ((k<$nOTHER))
                    do
                        if [ $(run_this_bin $k $kOTHER) -eq 1 ]
                        then
                            TEND=STYBKG_${COLSYST[j]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[l]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[k]})_leadingtrkptmin_$(float_to_string ${LEADING_TRKPTMIN[l]})
                            echo -e "-- Processing \033[1;33mstybkg_savefittpl.C${NC}, \033[1;32m${COLSYST[j]}${NC}, leading pT \033[1;32m> ${LEADING_PTMIN[l]} GeV${NC}, other pT \033[1;32m> ${OTHER_PTMIN[k]} GeV${NC}"
                            set -x
                            ./stybkg_savefittpl.exe "${INPUTSNAME[j]}" "rootfiles/ffittpl_${TEND}" "${COLSYST[j]}" ${LEADING_PTMIN[l]} ${OTHER_PTMIN[k]} ${LEADING_TRKPTMIN[l]}
                            set +x
                            echo
                        fi
                        k=$(($k+1))
                    done
                fi
                l=$(($l+1))
            done
        fi
        j=$(($j+1))
    done
    rm stybkg_savefittpl.exe
fi

# stybkg_savehist.C #
if [ $DO_STYBKG_SAVEHIST -eq 1 ]
then
    g++ stybkg_savehist.C $(root-config --cflags --libs) -g -o stybkg_savehist.exe || return 1
    j=0
    while ((j<$nCOL))
    do
        if [ $(run_this_bin $j $jCOLSYST) -eq 1 ]
        then
            l=0
            while ((l<$nLEAD))
            do
                if [ $(run_this_bin $l $lLEAD) -eq 1 ]
                then
                    k=0
                    while ((k<$nOTHER))
                    do
                        if [ $(run_this_bin $k $kOTHER) -eq 1 ]
                        then
                            TEND=STYBKG_${COLSYST[j]}_${TMC[${ISMC[j]}]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[l]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[k]})_leadingtrkptmin_$(float_to_string ${LEADING_TRKPTMIN[l]})
                            echo -e "-- Processing \033[1;33mstybkg_savehist.C${NC}, \033[1;32m${COLSYST[j]}${NC}, \033[1;32m${TMC[${ISMC[j]}]}${NC}, leading pT \033[1;32m> ${LEADING_PTMIN[l]} GeV${NC}, other pT \033[1;32m> ${OTHER_PTMIN[k]} GeV${NC}"
                            set -x
                            ./stybkg_savehist.exe "${INPUTDNAME[j]}" "rootfiles/fdphi_${TEND}" "${COLSYST[j]}" ${ISMC[j]} ${LEADING_PTMIN[l]} ${OTHER_PTMIN[k]} ${LEADING_TRKPTMIN[l]}
                            set +x
                            echo
                        fi
                        k=$(($k+1))
                    done
                fi
                l=$(($l+1))
            done
        fi
        j=$(($j+1))
    done
    rm stybkg_savehist.exe
fi

# stybkg_usehist.C #
if [ $DO_STYBKG_USEHIST -eq 1 ]
then
    g++ stybkg_usehist.C $(root-config --cflags --libs) -g -o stybkg_usehist.exe || return 1
    j=0
    while ((j<$nCOL))
    do
        if [ $(run_this_bin $j $jCOLSYST) -eq 1 ]
        then
            l=0
            while ((l<$nLEAD))
            do
                if [ $(run_this_bin $l $lLEAD) -eq 1 ]
                then
                    k=0
                    while ((k<$nOTHER))
                    do
                        if [ $(run_this_bin $k $kOTHER) -eq 1 ]
                        then
                            TEND=STYBKG_${COLSYST[j]}_${TMC[${ISMC[j]}]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[l]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[k]})_leadingtrkptmin_$(float_to_string ${LEADING_TRKPTMIN[l]})
                            FITTEND=STYBKG_${COLSYST[j]}_leadingDptmin_$(float_to_string ${LEADING_PTMIN[l]})_otherDptmin_$(float_to_string ${OTHER_PTMIN[k]})_leadingtrkptmin_$(float_to_string ${LEADING_TRKPTMIN[l]})
                            echo -e "-- Processing \033[1;33mstybkg_usehist.C${NC}, \033[1;32m${COLSYST[j]}${NC}, \033[1;32m${TMC[${ISMC[j]}]}${NC}, leading pT \033[1;32m> ${LEADING_PTMIN[l]} GeV${NC}, other pT \033[1;32m> ${OTHER_PTMIN[k]} GeV${NC}"
                            set -x
                            ./stybkg_usehist.exe "rootfiles/fdphi_${TEND}" "rootfiles/ffittpl_${FITTEND}" "${TEND}" "${COLSYST[j]}" ${ISMC[j]} ${LEADING_PTMIN[l]} ${OTHER_PTMIN[k]} ${LEADING_TRKPTMIN[l]}
                            set +x
                            echo
                        fi
                        k=$(($k+1))
                    done
                fi
                l=$(($l+1))
            done
        fi
        j=$(($j+1))
    done
    rm stybkg_usehist.exe
fi


