
scdir=/ourdisk/hpc/radclouds/auto_archive_notyet/tape_2copies/tc_ens

storm="haiyan"
# storm="maria"

# itest="ctl"

for itest in "ncrf36h" "STRATANVIL_OFF" "STRATANVIL_ON" "STRAT_OFF"; do # Ensemble member
for em in 0{1..9} 10; do # Ensemble member
#for em in 0{2..9} 10; do # Ensemble member

  outputdir="$scdir/$storm/memb_${em}/$itest/post/d02"
  echo " "
  echo "Running: $outputdir"
  cd $outputdir
  if [ -f RTHRATLWCRF_HiRes.nc ]; then
    echo "File found; skipping!"
    continue
  fi

cat << EOF >> batch_cdo_${itest}.job
#!/bin/bash
#SBATCH -p radclouds
#SBATCH -N 1
#SBATCH -n 1 ##112 # max of N*56
#SBATCH --exclusive
#SBATCH --output=outbatch_out.%j.txt
#SBATCH --error=outbatch_err.%j.txt
#SBATCH -t 48:00:00

module purge
source ~/.bashrc

cdo sub -selname,RTHRATLW RTHRATLW_HiRes.nc -selname,RTHRATLWC RTHRATLWC_HiRes.nc RTHRATLWCRFv1.nc
cdo chname,RTHRATLW,RTHRATLWCRF RTHRATLWCRFv1.nc RTHRATLWCRF_HiRes.nc
rm RTHRATLWCRFv1.nc
chmod a+xr RTHRATLWCRF_HiRes.nc
EOF

  sbatch batch_cdo_${itest}.job > submit_cdo_out.txt

  # cdo sub -selname,RTHRATLW RTHRATLW_HiRes.nc -selname,RTHRATLWC RTHRATLWC_HiRes.nc RTHRATLWCRFv1.nc
  # cdo chname,RTHRATLW,RTHRATLWCRF RTHRATLWCRFv1.nc RTHRATLWCRF_HiRes.nc
  # rm RTHRATLWCRFv1.nc
  # chmod a+xr RTHRATLWCRF_HiRes.nc

done
done
