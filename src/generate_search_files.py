import os
import sys
import shutil

#power of 10 being checked up to
p = int(sys.argv[1])
#numer of cpu cores to use
n_cores = int(sys.argv[2])
#number of machines to use
n_mach = int(sys.argv[3])

s = """#!/bin/bash
#SBATCH --account=def-cfranc
#SBATCH --time=28-00:00:00
#SBATCH --cpus-per-task={}
#SBATCH --mem-per-cpu=2048M
#SBATCH --output=self{}/{}.out

module load StdEnv/2020 julia/1.5.2
julia src/self.jl {} {} {} {}"""

if os.path.exists('./self{}'.format(p)):
    shutil.rmtree('./self{}'.format(p))
os.makedirs('./self{}'.format(p))

i = 0
while i < n_mach:
    with open('./self{}/{}.sh'.format(p, i), 'w') as f:
        f.write(s.format(n_cores, p, i, p, n_cores, n_mach, i))
    i += 1

with open('./self{}/run.sh'.format(p), 'w') as f:
    i = 0
    while i < n_mach:
        f.write('sbatch ./self{}/{}.sh\n'.format(p, i))
        f.write('sleep 1\n')
        i += 1
