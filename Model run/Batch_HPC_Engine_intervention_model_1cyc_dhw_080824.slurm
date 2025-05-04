#!/bin/bash -l
#Edit this script to suit your purposes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=500G
#SBATCH --job-name=HPC_Engine_intervention_model_1cyc_dhw_080824
#SBATCH --time=14:00:00
#SBATCH --partition=general
#SBATCH --account=a_msel
#SBATCH --array=1-36

# 1) Get access to the Matlab module
#module load matlab/2022a
module load matlab

# 2) Check status of Matlab license: will retry up to 5 times with a wait of x=2 minutes each time (so 10 mins in total)
#echo Check license status

#for iwait in `seq 1 5`
#do
 # /sw/Matlab/R2022a/etc/glnxa64/lmutil lmstat -c /sw/Matlab/R2022a/licenses/network.lic > /dev/null
  #The special shell variable $? will equal 0 if good, not 0 if bad.
  #if [ $? -eq 0 ]; then
   # echo "License is OK" 
   # break
  #else
   # echo "License status failed on attempt $iwait" 
    #wait x minutes each time
    #sleep 2m
  #fi
#done

# 3) Specify the working directory for the job:
#echo /scratch/user/uqaadam5

cd $TMPDIR 
#cd /scratch/user/uqaadam5/Deployment/model_070723_temp/

echo Copying the current directory to TMPDIR

cp -rp /scratch/user/uqaadam5/Deployment/model_080124/* ./
# IMPORTANT The code above copies the rme v1.0.28 to temporary directory to run the model. Folder "model_Adapt_130624" contains folder "rme_rml" which is the linux version of rme v1.0.33


# 4) Run Matlab code
echo HPC_Engine_intervention_model_1cyc_dhw_080824

# YM: Note with singleCompThread we force execution to a single CPU, so ncpus >1 doesn't make sense?
matlab -singleCompThread -nodisplay -nosplash -r "HPC_Engine_intervention_model_1cyc_dhw_080824($SLURM_ARRAY_TASK_ID)"

# ./run_EXTRACTION_currents_GBR4.sh sw/Matlab/2018a $SLURM_ARRAY_TASK_ID

# 5) Take output files (.mat only!) to the working directory
echo Copying results back ...
cp *.mat /QRISdata/Q5785/Virtual_reef/version_080124/VR_1dist_1int_cyc_080824/
