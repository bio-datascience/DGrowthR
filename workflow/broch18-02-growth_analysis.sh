#!/bin/bash

#------- Job Description -------

#SBATCH -J 'dgrowthr_broch18'
#SBATCH -A 'maccf'

#------- Parametrization -------

#SBATCH --requeue
#SBATCH --nodes=1 # Number of nodes
#SBATCH -n 20 # Number of cores
#SBATCH --ntasks=1 # Number of tasks (processes)
#SBATCH --cpus-per-task=20
#SBATCH --partition=bigmem
#SBATCH --time=3-0:0:0
#SBATCH --mem=300G # memory in MB
#SBATCH --mail-user=roberto.olayoalarcon@embl.de
#SBATCH --mail-type=ALL

#------- Input/Output -------

#SBATCH --output="/home/alarcon/DGrowthR/workflow/logdir/broch18/%x-%j.out"
#SBATCH --error="/home/alarcon/DGrowthR/workflow/logdir/broch18/%x-%j.err"

#------- Command -------
module load R
Rscript broch18-02-growth_analysis.R
