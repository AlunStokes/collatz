import os
import sys
import shutil

#power of 10 being checked up to
k_str = sys.argv[1]
K = k_str.split(',')
K = [int(k) for k in K]
a = int(sys.argv[2])
#numer of cpu cores to use
n_cores = int(sys.argv[3])

s = """#!/bin/bash
#SBATCH --account=def-cfranc
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task={}
#SBATCH --mem-per-cpu=2048M
#SBATCH --output=neg_inv/k{}a{}.out

module load StdEnv/2020 julia/1.5.2
julia src/neg_inv.jl {} {} {}"""

if not os.path.exists('./neg_inv'):
    os.makedirs('./neg_inv')

for k in K:
    with open('./neg_inv/k{}a{}.sh'.format(k, a), 'w') as f:
        f.write(s.format(n_cores, k, a, k, a, n_cores))

with open('./neg_inv/run.sh', 'w') as f:
    for k in K:
        f.write('sbatch ./neg_inv/k{}a{}.sh\n'.format(k, a))
        f.write('sleep 1\n')
