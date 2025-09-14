#!/bin/bash

main_dir=`pwd`

# Specify the directory to work with
for em in 0{1..9} 10; do # Ensemble member
#for em in 01; do # Ensemble member

  memdir="memb_${em}"
  cd $main_dir/$memdir

  for test in ./*; do

    # Check if it's a file (not a directory)
    cd "$test/post/d02"
    # mv "rainrate_HiRes.nc" "rainrate.nc"

  done
done