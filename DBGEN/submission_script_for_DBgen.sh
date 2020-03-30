#!/bin/bash
# -*- coding: utf-8 -*-

"""
Submit script for DBgen.py to create conformers, gaussian input and submit jobs
@author: Shree Sowndarya S. V.
"""

#just for runnning the xtb jobs
#SBATCH --account=csu-general
#SBATCH -J xtb-and gaussian-submission
#SBATCH -p shas
#SBATCH --qos normal
#SBATCH -t 23:59:59
#SBATCH -N 1
#SBATCH --export=NONE
#SBATCH --ntasks-per-node 24

run xTB
get info

#pass the respective file as the first argument
#running the compute job for either csv, smi, sdf, cdx
python /DBcg/db_gen.py --compute --input $0

#now after the above step all jobs would have been submitted.
#check for any non normal termination using the output analysis
python /DBcg/db_gen.py --analysis --path $(PWD)/gaussian/

##resubmit the gaussian jobs if the didnt normally terminate
python /DBcg/db_gen.py --resubmit --path $(PWD)/gaussian/

##once resubmission is completed analyse the NEWLY CREATED Files and moves the normally terminated to the finished folder
python /DBcg/db_gen.py --analysis --secondrun --path $(PWD)/gaussian/

##creating the NMR input files from finished log files i.e., for doing NMR after DFT optimization
python /DBcg/db_gen.py --nmr --path $(PWD)/gaussian/

##carrying out the boltzmann from goodvibes
##note change the name of the log files accoring to your description in the reading of log files part
python /DBcg/db_gen.py --boltz --path $(PWD)/gaussian/

##combing all the energies for all molecules and printing the output in one csv_file
python /DBcg/db_gen.py --combine --path $(PWD)/gaussian/
