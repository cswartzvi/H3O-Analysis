#! /bin/bash



PRO_DIR='/home/tua01780/MD_Simulations/AIMD/Analysis_Tools/OHminus_analysis/int-KS'
istart=$1
iend=$2
file_root='KS_wf_n'


echo " "
rm test.log

for ((i = $istart; i <= $iend; i++))
do

   echo "$i $($PRO_DIR'/KS' < $file_root$i.xsf)" >> test.log

   

done
