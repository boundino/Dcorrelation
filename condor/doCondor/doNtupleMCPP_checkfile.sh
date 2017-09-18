#!/bin/bash

if [[ $0 != *.sh ]]
then
    echo -e "\e[31;1merror:\e[0m use \e[32;1m./script.sh\e[0m instead of \e[32;1msource script.sh\e[0m"
    return 1
fi

iPTHAT=0

#
isMC=1
isPP=1
MAXFILENO=100000000
ifCHECKEMPTY=0

#
movetosubmit=${1:-0}
runjobs=${2:-0}
merger=${3:-0}

#
PTHATS=(0 5 10 15 30 50 80 120 170)
DATES=(170907_172052 170908_155109 170908_155559 170908_155644 170908_155726 170908_155800 170908_155833 170908_155912 170908_155944)

INPUTDIR="/mnt/hadoop/cms/store/user/wangj/Pythia8_prompt_D0pt0p0_Pthat${PTHATS[iPTHAT]}_pp502_TuneCUETP8M1/crab_HiForestAOD_DfinderMC_pp_20170906_Pthat${PTHATS[iPTHAT]}_dPt0tkPt0p5_Dhadron/${DATES[iPTHAT]}/0000/"
OUTPUTPRIDIR="/mnt/hadoop/cms/store/user/wangj/condor/Dtrk/"
OUTPUTSUBDIR="DtrkFiles_20170918_pp_5TeV_TuneCUETP8M1_Dfinder_MC_Pthat${PTHATS[iPTHAT]}_20170906"

MERGEOUTPUTDIR="/export/d00/scratch/jwang/Dtrk/MC/"
WORKDIR="/work/$USER/Dtrk/"

##
OUTPUTDIR="${OUTPUTPRIDIR}/${OUTPUTSUBDIR}"
LOGDIR="logs/log_${OUTPUTSUBDIR}"
RESIDUALS="Corrections"

#
if [[ ! -d "$WORKDIR" ]]
then
    mkdir -p $WORKDIR
fi

if [ "$movetosubmit" -eq 1 ]
then 
    if [[ $(hostname) == "submit-hi2.mit.edu" || $(hostname) == "submit.mit.edu" ]]
    then
        cd ../skim
        g++ D_track_skim.C $(root-config --cflags --libs) -Werror -Wall -O2 -o D_track_skim.exe
        tar -czf ${RESIDUALS}.tar.gz ${RESIDUALS}/
        cd ../condor

        mv ../skim/D_track_skim.exe $WORKDIR/
        mv ../skim/${RESIDUALS}.tar.gz $WORKDIR/
        cp $0 $WORKDIR/
        cp skim_dtrk_checkfile.sh $WORKDIR/
        cp skim_condor_checkfile.sh $WORKDIR/
    else
        echo -e "\e[31;1merror:\e[0m compile macros on \e[32;1msubmit-hiX.mit.edu\e[0m or \e[32;1msubmit.mit.edu\e[0m."
    fi
fi

if [ "$runjobs" -eq 1 ]
then
    if [[ $(hostname) == "submit.mit.edu" ]]
    then
        ./skim_condor_checkfile.sh $INPUTDIR $OUTPUTDIR $MAXFILENO $LOGDIR $isPP $isMC $ifCHECKEMPTY ${RESIDUALS}.tar.gz
    else
        echo -e "\e[31;1merror:\e[0m submit jobs on \e[32;1msubmit.mit.edu\e[0m."
    fi
fi

if [ "$merger" -eq 1 ]
then
    if [[ $(hostname) == "submit-hi2.mit.edu" ]]
    then
        if [[ ! -d "$MERGEOUTPUTDIR" ]]
        then
            mkdir -p $MERGEOUTPUTDIR
        fi

        rm ${MERGEOUTPUTDIR}/${OUTPUTSUBDIR}.root
        cp -r ${OUTPUTPRIDIR}/${OUTPUTSUBDIR} ${MERGEOUTPUTDIR}/
        hadd ${MERGEOUTPUTDIR}/${OUTPUTSUBDIR}.root ${MERGEOUTPUTDIR}/${OUTPUTSUBDIR}/*.root
        rm -r ${MERGEOUTPUTDIR}/${OUTPUTSUBDIR}
    else
        echo -e "\e[31;1merror:\e[0m merge files on \e[32;1msubmit-hiX.mit.edu\e[0m."
    fi
fi
