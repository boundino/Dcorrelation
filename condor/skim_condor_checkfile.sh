#!/bin/bash

if [[ $# -ne 8 ]]; then
    echo "usage: ./skim_condor_checkfile.sh [input dir] [output dir] [max jobs] [log dir] [isPP] [isMC] [check input empty] [residuals]"
    exit 1
fi

DATASET=$1
DESTINATION=$2
MAXFILES=$3
LOGDIR=$4
isPP=$5
isMC=$6
ifCHECKEMPTY=$7
RESIDUALS=$8

PROXYFILE=$(ls /tmp/ -lt | grep $USER | grep -m 1 x509 | awk '{print $NF}')

OUTFILE="skim_dtrk"

rm filelist.txt
ls $DATASET | grep -v "/" | grep -v -e '^[[:space:]]*$' | awk '{print "" $0}' >> filelist.txt

SRM_PREFIX="/mnt/hadoop/"
SRM_PATH=${DESTINATION#${SRM_PREFIX}}

if [ ! -d $DESTINATION ]
then
    # gfal-mkdir -p srm://se01.cmsaf.mit.edu:8443/srm/v2/server?SFN=$DESTINATION
    gfal-mkdir -p gsiftp://se01.cmsaf.mit.edu:2811/${SRM_PATH}
fi

if [ ! -d $LOGDIR ]
then
    mkdir -p $LOGDIR
fi

counter=0
for i in `cat filelist.txt`
do
    if [ $counter -ge $MAXFILES ]
    then
        break
    fi
    #ifexist=`ls ${DESTINATION}/skim_djet_$i`
    #if [ -z $ifexist ]
    if [ ! -f ${DESTINATION}/${OUTFILE}_$i ] && [ -f ${DATASET}/$i ]
    then
        if [ -s ${DATASET}/$i ] || [ $ifCHECKEMPTY -eq 0 ]
        then
            echo -e "\033[38;5;242mSubmitting a job for output\033[0m ${DESTINATION}/${OUTFILE}_$i"
            infn=`echo $i | awk -F "." '{print $1}'`
            INFILE="${DATASET}/$i"
	    
	    cat > skim-dtrk.condor <<EOF

Universe     = vanilla
Initialdir   = $PWD/
Notification = Error
Executable   = $PWD/skim_dtrk_checkfile.sh
Arguments    = $INFILE $DESTINATION ${OUTFILE}_${infn}.root $isPP $isMC $PROXYFILE $RESIDUALS
GetEnv       = True
Output       = $LOGDIR/log-${infn}.out
Error        = $LOGDIR/log-${infn}.err
Log          = $LOGDIR/log-${infn}.log
Rank         = Mips
+AccountingGroup = "group_cmshi.$(whoami)"
requirements = GLIDEIN_Site == "MIT_CampusFactory" && BOSCOGroup == "bosco_cmshi" && HAS_CVMFS_cms_cern_ch && BOSCOCluster == "ce03.cmsaf.mit.edu"
job_lease_duration = 240
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = D_track_skim.exe,/tmp/$PROXYFILE,$RESIDUALS

Queue 
EOF

condor_submit skim-dtrk.condor -name submit.mit.edu
# condor_submit -pool submit.mit.edu:9615 -name submit.mit.edu -spool skim-djet.condor
mv skim-dtrk.condor $LOGDIR/log-${infn}.condor
counter=$(($counter+1))	
        fi
    fi
done

echo -e "Submitted \033[1;36m$counter\033[0m jobs to Condor."
