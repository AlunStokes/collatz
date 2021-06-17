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

s = """#!/bin/bash
#SBATCH --account=def-cfranc
#SBATCH --time=28-00:00:00
#SBATCH --cpus-per-task={}
#SBATCH --mem-per-cpu=2048M
#SBATCH --output=self.out

module load StdEnv/2020 julia/1.5.2
julia src/self.jl"""

files = os.listdir('./self{}'.format(p))

T = []
N = []
n = 0

for file in files:
    if '.out' not in file:
        continue
    if file == 'run.sh':
        continue
    with open('./self{}/{}'.format(p, file), 'r') as f:
        r = f.read()
    if 'TIME LIMIT' in r:
        print('{} ran out of time'.format(file))
        continue
    R = r.split('\n')
    index = 0
    while index < len(R):
        if R[index] == "TEST PRINT":
            index += 2
            break
        index += 1
    if len(R) <= index:
        print('{} is not finished yet'.format(file))
        continue
    n += 1
    s_time = R[index]
    time = s_time.split('.')[0]
    time = int(time)
    T.append(time)
    index += 1
    while index < len(R) - 1:
        s_num = R[index]
        s_num = int(s_num)
        N.append(s_num)
        index += 1

print('Took {} computers an average of {}s (max {}s ; min {}s) each'.format(n, int(mean(T)), min(T), max(T)))
if len(N) > 0:
    print('Found the following self-contained numbers:')
    for n in N:
        print('\t{}'.format(n))
else:
    print('No new self-contained numbers found.')
