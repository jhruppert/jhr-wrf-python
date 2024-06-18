#!/bin/bash

owd=`pwd`

dir='/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens/haiyan/memb_01/ctl/post/d02'

cd $dir

# All
for pclass in noncloud deepc congest shallowc strat anvil; do # Ensemble member
#for em in 03; do # Ensemble member

  infile="binned_isentrop_${pclass}_HiRes.nc"
  outfile="binned_isentrop_${pclass}_HiRes_out.nc"
  ls $infile

  # for varname in tmpk_mean,qv_mean,rho_mean,lw_mean,lwc_mean,sw_mean,swc_mean; do

  # echo $varname
  ncks -x -v tmpk_mean,qv_mean,rho_mean,lw_mean,lwc_mean,sw_mean,swc_mean,w_mean $infile $outfile
  mv $outfile $infile

  # done

done

cd $owd

exit
