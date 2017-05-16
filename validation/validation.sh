#!/bin/bash

EWKfast='/Users/kazuki/Projects/EWKfast/EWKfast/EWKfast.py'
#EWKfast='/Users/kazuki/Projects/EWKfast/EWKfast_debug/EWKfast.py'

order=$1
rs=$2

if [[ $rs == '' ]]; then
    echo '[order] [rs]'
    exit
fi

tag=$order$rs

for i in {1..600}; do 
#for i in {1..2}; do 
    echo $i; 
    input_file='/Users/kazuki/Projects/EWKfast/validation/input_'$tag'/input_'$i'.dat'
    slha_file='/Users/kazuki/Projects/EWKfast/validation/slha_files/slha-'$i'.dat'
    outfile='/Users/kazuki/Projects/EWKfast/validation/EWKFresult_'$tag'/result-'$i'.dat'
    #outfile='/Users/kazuki/Projects/EWKfast/validation/EWKFresult_'$tag'_nearest/result-'$i'.dat'

    rm -f 'output.ewk'
    $EWKfast $slha_file $input_file > /dev/null    
    cp 'output.ewk' $outfile

done