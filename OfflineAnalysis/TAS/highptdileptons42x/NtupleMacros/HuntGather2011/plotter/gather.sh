#!/bin/bash

#
# This script makes the plots
# it should be run as a cron job
#

cd /home/users/dlevans/gather/production/CMSSW_4_1_2_patch1/src/CMS2/NtupleMacros/HuntGather2011/plotter
export SCRAM_ARCH=slc5_amd64_gcc434
source /code/osgcode/cmssoft/cms/cmsset_default.sh
eval `scramv1 runtime -sh`

echo "getting new dcs run list"
curl http://dlevans.web.cern.ch/dlevans/dcs.txt > ../runlists/dcs.txt
python ../../Tools/convertGoodRunsList_JSON.py ../runlists/dcs.txt > ../runlists/dcs_jmu.txt
echo $?

root -q -b makeGatherPlots.C\(\"/nfs-3/userdata/cms2/gather/\"\)
cp /home/users/dlevans/gather/production/CMSSW_4_1_2_patch1/src/CMS2/NtupleMacros/HuntGather2011/output/run2011*.png /home/users/dlevans/public_html/gather2011/
cp /home/users/dlevans/gather/production/CMSSW_4_1_2_patch1/src/CMS2/NtupleMacros/HuntGather2011/output/max2011*.png /home/users/dlevans/public_html/gather2011/
touch /home/users/dlevans/public_html/gather2011/index.php

