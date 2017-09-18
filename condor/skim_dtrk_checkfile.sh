#!/bin/bash

if [[ $# -ne 7 ]]; then
    echo "usage: ./skim-djet.sh [input file] [output dir] [output filename] [isPP] [isMC] [proxy] [residuals]"
    exit 1
fi

INFILE=$1
DESTINATION=$2
OUTFILE=$3
isPP=$4
isMC=$5
export X509_USER_PROXY=${PWD}/$6
RESIDUALS=$7

SRM_PREFIX="/mnt/hadoop/"
SRM_PATH=${DESTINATION#${SRM_PREFIX}}

tar -xzf $RESIDUALS

echo ./D_track_skim.exe $INFILE $OUTFILE $isPP $isMC
./D_track_skim.exe $INFILE $OUTFILE $isPP $isMC

if [[ $? -eq 0 ]]
then
    # gfal-copy file://${PWD}/${OUTFILE}  srm://se01.cmsaf.mit.edu:8443/srm/v2/server?SFN=${DESTINATION}/${OUTFILE}
    gfal-copy file://$PWD/${OUTFILE} gsiftp://se01.cmsaf.mit.edu:2811/${SRM_PATH}/${OUTFILE}
fi
