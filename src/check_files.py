import os
import sys
import shutil

#power of 10 being checked up to
p = int(sys.argv[1])

s = """#!/bin/bash
#SBATCH --account=def-cfranc
#SBATCH --time=28-00:00:00
#SBATCH --cpus-per-task={}
#SBATCH --mem-per-cpu=2048M
#SBATCH --output=self.out

module load StdEnv/2020 julia/1.5.2
julia src/self.jl"""

files = os.listdir('./self{}'.format(p))

for file in files:
    if '.sh' not in file:
        continue
    if file == 'run.sh':
        continue
    with open('./self{}/{}'.format(p, file), 'r') as f:
        r = f.read()
        if r.count('\n') != s.count('\n'):
            print("Error with file {}".format(file))
