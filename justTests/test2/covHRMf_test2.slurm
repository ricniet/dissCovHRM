#!/bin/bash
#SBATCH -J covHRM_test		     # Job name
#SBATCH -o covHRM_test.o%j     # Name of stdout output file (%j expands to jobId)
#SBATCH -e covHRM_test.o%j     # Name of stderr output file(%j expands to jobId)
#SBATCH -n 48                 	# Total number of mpi tasks requested
#SBATCH -N 1
#SBATCH -t 48:00:00                # Run time (hh:mm:ss)
#SBATCH -p skx-normal
#SBATCH --mail-user=nietoric@utexas.edu # specify email address for notifications
#SBATCH --mail-type=all            # email when job ends

ml rm Rstats/3.4.0 RstatsPackages/3.4.0

# Get hostnames from SLURM
scontrol show hostnames $SLURM_NODELIST > nodelist.txt
scontrol show hostnames $SLURM_NODELIST >> nodelist.txt
scontrol show hostnames $SLURM_NODELIST >> nodelist.txt

# call R code

Rscript covHRMf_onlyRho_RSM_model_tacc_test2.R nodelist.txt 2 covHRMf_onlyRho 400N30R8J5K 3ch20Ki10Kb20t &> log
