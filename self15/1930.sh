#!/bin/bash
#SBATCH --account=def-cfranc
#SBATCH --time=0-09:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3072M
#SBATCH --output=self15/1930.out

module load StdEnv/2020 julia/1.5.2
julia src/self.jl 15 48 2400 1930