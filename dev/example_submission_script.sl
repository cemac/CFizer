#!/bin/bash
#SBATCH --partition=test
#SBATCH -n 1
#SBATCH -c 4
##SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH --job-name=<cfize>
#SBATCH -o /<run_directory>/%j.out
#SBATCH -e /<run_directory>/%j.err
#SBATCH --time=60
# #SBATCH --account=short4hr
# #SBATCH --account=<account>
##SBATCH --mem=16384
#SBATCH --mem-per-cpu=16384

module add jaspy/3.10/r20230718
srun ./cfize.py -p 4 -v -i <directory_containing_nc_files>

