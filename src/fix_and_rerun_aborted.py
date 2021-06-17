import os
import sys
import shutil

def mean(L):
    m = 0
    for i in L:
        m += i
    m /= len(L)
    return m

#power of 10 being checked up to
p = int(sys.argv[1])
#new timecode
new_time = sys.argv[2]

s = """#!/bin/bash
#SBATCH --account=def-cfranc
#SBATCH --time=0-09:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3072M
#SBATCH --output=self15/0.out

module load StdEnv/2020 julia/1.5.2
julia src/self.jl 15 48 2400 0"""

files = os.listdir('./self{}'.format(p))

A = []

for file in files:
    if '.out' not in file:
        continue
    with open('./self{}/{}'.format(p, file), 'r') as f:
        r = f.read()
    if 'TIME LIMIT' in r:
        A.append(file.split('.')[0])

for a in A:
    with open('./self{}/{}.sh'.format(p, a), 'r') as f:
        r = f.read()
    R = r.split('\n')
    if R[-1] == '':
        del R[-1]
    for i, s in enumerate(R):
        if not '--time' in s:
            continue
        index = s.index('=')
        s = s[:index + 1]
        s += new_time
        R[i] = s
        break
    new_text = ""
    for s in R:
        new_text += s + '\n'
    with open('./self{}/{}.sh'.format(p, a), 'w') as f:
        f.write(new_text)
    os.remove('./self{}/{}.out'.format(p, a))

run_text = ""
for a in A:
    run_text += 'sbatch ./self{}/{}.sh\n'.format(p, a)
    run_text += 'sleep 2\n'
with open('./self{}/rerun.sh'.format(p), 'w') as f:
    f.write(run_text)
