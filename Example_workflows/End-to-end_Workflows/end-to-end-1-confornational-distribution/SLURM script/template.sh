#!/bin/bash

#SBATCH -J Strychnine
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH --export=NONE
#SBATCH --ntasks-per-node=24

##################################################
################    MODULES     ##################
##################################################

source $HOME/.bash_profile

g16root="/opt/apps/Gaussian/G16C"
export g16root
. $g16root/g16/bsd/g16.profile
export GAUSS_SCRDIR="/tmp"

source activate DL_CPU

##################################################
################    VARIABLES    #################
##################################################

childcpu=12 # Number of cpus per parallel job during batching
totalcpu=$SLURM_NTASKS_PER_NODE
jbatch=$(expr $totalcpu / $childcpu)
excdir=$(pwd)
jobname=QCALC

##################################################
################    FUNCTIONS    #################
##################################################

### Run QM calculation. In theory any program could be called here.
task () {

$g16root/g16/g16 < $1.com >> $2.log

}

### Wrapper for batching  and running jobs
gexc () {
#Variables; Use local here or otherwise it becomes global and overwrites original
local jobname=$1 # The current job
local jbatch=$2 # Number or parallel jobs.

if [ -d $jobname ]; then
        cd $jobname # We don't exit this directory within the function. Might cause issues.
else
        echo "No $jobname folder found, exiting the script"
        exit
fi
start=$(date +%s)
for i in $(ls -1vd *.com); do
        (
        file=$(echo $i | sed 's/.com//')
        echo "$file"
        task $file $excdir/$jobname/$file
        ) &
        # allow to execute up to $N jobs in parallel
        if [[ $(jobs -r -p | wc -l) -ge $jbatch ]]; then
        # now there are $N jobs already running, so wait here for any job
        # to be finished so there is a place to start next one.
                wait
        fi
done
wait
end=$(date +%s)
sec=$(expr $end - $start)
echo "Time $jobname: $sec seconds"

}


##################################################
################    MAIN SCRIPT     ##############
##################################################


### Run AQME-CSEARCH to generate conformers with RDKit from a list of SMILES in a CSV
python -m aqme --csearch --program rdkit --input smiles.csv --nprocs "$totalcpu"

### Run AQME-QPREP to generate input files for Gaussian
python -m aqme --qprep --program gaussian --files CSEARCH/rdkit/*sdf --qm_input 'wB97xd/def2TZVPP scrf=(solvent=chloroform,smd) opt freq' --mem '24GB' --nprocs "$childcpu"

### Run Gaussian with the generated COM files
echo "Starting opt-freq calculations"
gexc $jobname $jbatch
wait
touch JOB_DONE
cd $excdir

### Run AQME-QCORR to analyze the LOG files from Gaussian
python -m aqme --qcorr --w_dir_main "$excdir/$jobname" --program gaussian --files "*.log" --fullcheck --mem '16GB' --nprocs "$childcpu" 

### Run Gaussian with the COM files that were generated from unsuccessful calculations, and then run AQME-QCORR to analyze the new LOG files
fixed="unsuccessful_QM_outputs/run_1/fixed_QM_inputs"
if [ -d "$excdir/$jobname/$fixed" ]; then
    echo "Starting one rerun of failed calculations"
            gexc $excdir/$jobname/$fixed $jbatch
            wait
            cd $excdir
            python -m aqme --qcorr --w_dir_main "$excdir/$jobname/$fixed" --files "*.log" --fullcheck --mem '16GB' --nprocs "$childcpu" 
    else
            echo "No fixable errors found"
    fi
    
echo "Done with opt-freq calculations"
touch QCORR_DONE
cd $excdir

### Run GoodVibes with the processed LOG files
successful_QM_outputs='successful_QM_outputs'
GoodVibes_Analysis='GoodVibes_Analysis'

mkdir $excdir/$GoodVibes_Analysis

cp $excdir/$jobname/$successful_QM_outputs/*log $excdir/$GoodVibes_Analysis

python -m goodvibes --boltz --xyz $excdir/$GoodVibes_Analysis/*.log

cd $excdir
