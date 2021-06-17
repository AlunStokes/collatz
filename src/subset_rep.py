import time

from math import floor
from multiprocessing.pool import Pool
from functools import partial
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from itertools import combinations, permutations, product

from utils import *
from scaling_factors import get_wab
from garb import nu
from scaling_factors import *
from test_r import get_all_r

def prime_fac(n):
    if n < 0:
        n = -n
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors

def is_power_of_p(n, p):
    a = 0
    while n % p**a == 0:
        a += 1
    a -= 1
    if n == p**a:
        return True
    return False

def tuple_diff(a, b):
    return tuple([i - j for i, j in zip(a, b)])


def seq_diff(a):
    b = [a[0]] + [a[i] - a[i - 1] for i in range(1, len(a))]
    return tuple(b)

#Uses 0 indexing
def inv_r(r, n):
    if n == 0:
        if r == 0:
            return tuple()
        else:
            return (-1,)
    if n == 1:
        v = nu(r, 2)
        if 2**v == r:
            return (v,)
        else:
            return (-1,)
    L = []
    i = 0
    while i < n:
        s = 0
        j = 0
        while j < i:
            s += 3**(n-(j+1)) * 2**(L[j] - 1)
            j += 1
        L.append(nu(r - s) + 1)
        i += 1
    return tuple([i - 1 for i in L])

#Uses 0 indexing
def r_func(A, right_offset=0):
    r = 0
    for i, v in enumerate(A[::-1]):
        r += 2**(v) * 3**(i + right_offset)
        #r += 3**(len(A) - 1 - i - offset) * 2**(A[i])
    return r

#Considers which parts of original rep are given by comb
def r_func_comb(A, comb, total_len=-1):
    if total_len == -1:
        total_len = len(A)
    r = 0
    for c, v in zip(comb[::-1], A[::-1]):
        r += 2**(v) * 3**(total_len - c - 1)
        #r += 3**(len(A) - 1 - i - offset) * 2**(A[i])
    return r

def r_func_comb(A, comb):
    total_len = len(A)
    A = [A[i] for i in comb]
    r = 0
    for c, v in zip(comb[::-1], A[::-1]):
        r += 2**(v) * 3**(total_len - c - 1)
        #r += 3**(len(A) - 1 - i - offset) * 2**(A[i])
    return r

def get_all_r(a, b, include_even = False, starting_with=tuple()):
    R = []
    Is = []
    if len(starting_with) > 0:
        comb = [(*starting_with, *c) for c in list(combinations(list(range(starting_with[-1] + 1, b + 1)), a - len(starting_with)))]
    elif include_even:
        comb = list(combinations(list(range(0, b + 1)), a))
    else:
        comb = [(0, *c) for c in list(combinations(list(range(1, b + 1)), a - 1))]
    if len(starting_with) > 0:
        starting_length = len(starting_with)
        comb = [c for c in comb if c[:starting_length] == starting_with]
    for I in comb:
        R.append(r_func(I))
        Is.append(I)
    return R, Is

def get_all_comb_sum(v1, total, max_len):
    C = []
    for l in range(1, max_len):
        nums = list(range(1, total - v1 - l + 2))
        result = [(v1, *seq) for seq in product(nums, repeat=l) if sum(seq) == total - v1]
        C += result
    C = sorted(list(set(C)))
    D = {}
    for c in C:
        if not len(c) in D:
            D[len(c)] = [c]
        else:
            D[len(c)].append(c)
    return C, D

def get_matching_patterns(I, Maps):
    Patterns = []
    for pattern in Maps:
        match = True
        sub_index = pattern.index('_')
        v = I[sub_index]
        if v not in Maps[pattern][1]:
            match = False
        i = sub_index + 1
        while i < len(pattern):
            if pattern[i] == '*':
                i += 1
                continue
            v = v + pattern[i]
            if I[i] != v:
                match = False
                break
            i += 1
        if match:
            Patterns.append(pattern)
    return Patterns

def apply_pattern(I, pattern, Maps):
    #Assume pattern applies, no checking
    Is = []
    sub_index = pattern.index('_')
    v = I[sub_index]
    if not v in Maps[pattern][1]:
        return []
    Res_patterns = Maps[pattern][0]
    for res_pattern in Res_patterns:
        J = [v]
        for j in res_pattern[1:]:
            J.append(J[-1] + j)
        J = tuple(J)
        Is.append(J)
    return Is

def pattern_to_comb(pattern):
    c = []
    for i, j in enumerate(pattern):
        if j != '*':
            c.append(i)
    return c

def get_Idw(check_d, w):
    a_offset = 0
    while True:
        a = max(5, check_d - 2) + a_offset
        R, Is = get_all_r(a, 2*a + 1, include_even=False)
        L = []
        for r, I in zip(R, Is):
            r_ = r + 3**a * w
            d = nu(r, 2)
            d_ = nu(r_, 2)
            #d = nu(r, 4)
            #d_ = nu(r_, 4)
            if d_ == check_d:
                L.append(I)
        if len(L) in [0, 1]:
            a_offset += 1
            continue
        a_offset = 0
        matches_up_to = -1
        for i in range(a):
            match = True
            v = L[0][i]
            for I in L:
                if I[i] != v:
                    match = False
                    break
            if match:
                matches_up_to += 1
            else:
                break
        break
    return L[0][:matches_up_to + 1]

if __name__ == '__main__':

    #Check for single d
    '''w = 9
    a = 8
    R, Is = get_all_r(a, 2*a + 5, include_even=False)
    for r, I in zip(R, Is):
        r_ = r + 3**a * w
        d = nu(r, 2)
        d_ = nu(r_, 2)
        #if not d == min(I):
            #print('{}: {} | {}'.format(r, nu(r, 2), I))
        if d_ == 9:
            print(I)
        #print('{} ({}) {} -> {} ({})'.format(r, d, I, r_, d_))
    exit()'''

    #BIG GUY - KEEP ME
    #This finds the mappings for subsequences
    a = 5
    b = 2 * a + 5

    '''t = 15
    for v1 in range(0, 6):
        C, D = get_all_comb_sum(v1, v1 + t, a)
        print('{}: {}'.format(v1, C))
    exit()'''

    R, Is = get_all_r(a, b)

    combs = []
    for k in range(1, a + 1):
        combs += list(combinations(range(a), k))

    #Finds ALL subsets
    D = {}
    for r, I in zip(R, Is):
        #print('{}: {}'.format(r, I))
        for comb in combs:
            s = tuple([I[c] for c in comb])
            rs = r_func_comb(I, comb=comb)
            #print('  {}, {} = {}'.format(s, right_offset, rs))
            for m in range(1, a+3):
                if in_Rn(rs, m):
                    Im = inv_r(rs, m)
                    if Im == s:
                        continue
                    k = (s, tuple(comb))
                    if k in D:
                        if Im not in D[k]:
                            D[k].append(Im)
                    else:
                        D[k] = [Im]
                    #print('    m = {}: {}'.format(m, Im))
        #print('')

    print('found all subsets')
    #P[starting digit of original seq, sum of remaining differences, pattern of remaining differences] = [difference sequences mapped to by this]
    P = {}
    #Read through mappings
    for k in D:
        #if maps to multiple others
        #if len(D[k]) > 1:
        diff = seq_diff(k[0])
        t = diff[0]
        s = sum(diff[1:])
        if t not in P:
            P[t] = {}
        if not s in P[t]:
            P[t][s] = {}
        if not diff[1:] in P[t][s]:
            L = []
            for d in D[k]:
                d_diff = seq_diff(d)
                if d_diff[0] == t:
                    d_diff = ('_', *d_diff[1:])
                L.append(d_diff)
            P[t][s][diff[1:]] = (L, k[1])
        #print('k = {} | {}'.format(k, diff))
        for d in D[k]:
            #print('  {} | {}'.format(d, seq_diff(d)))
            pass
        #print()
    print('partitioning completed')
    F = {}
    #T tracks where each diff seq gets mapped and by which starting num
    T = {}
    #C Tracks which starting num have a seq of each pattern getting mapped
    C = {}
    for t in P.keys():
        print('Ic starting with {}'.format(t))
        for s in P[t].keys():
            if s not in F:
                F[s] = []
            F[s].append(len(P[t][s].keys()))
            print('  D(I) has sum {} ({} entrie(s))'.format(s, len(list(P[t][s].keys()))))
            maxlentup = 0
            maxlenlist = 0
            for diff in P[t][s]:
                #tupstr = str(('_', *diff))
                tupstr = []
                j = 0
                for i in range(a):
                    if i == P[t][s][diff][1][0]:
                        tupstr.append('_')
                    elif i not in P[t][s][diff][1]:
                        tupstr.append('*')
                    else:
                        tupstr.append(diff[j])
                        j += 1
                tupstr = str(tuple(tupstr))
                liststr = str(P[t][s][diff][0])
                if len(tupstr) > maxlentup:
                    maxlentup = len(tupstr)
                if len(liststr) > maxlenlist:
                    maxlenlist = len(liststr)
            for diff in P[t][s]:
                #D(Ic)[1:] -> D(p)
                #if len(P[t][s][diff][1]) == 1:
                tupstr = []
                j = 0
                for i in range(a):
                    if i == P[t][s][diff][1][0]:
                        tupstr.append('_')
                    elif i not in P[t][s][diff][1]:
                        tupstr.append('*')
                    else:
                        tupstr.append(diff[j])
                        j += 1
                tup = tuple(tupstr)
                tupstr = str(tuple(tupstr))
                liststr = str(P[t][s][diff][0])
                print('    {}{} -> {}{} ({}) | {}'.format(tupstr, ' ' * (maxlentup - len(tupstr)), liststr, ' ' * (maxlenlist - len(liststr)), [sum(x[1:]) for x in P[t][s][diff][0]], P[t][s][diff][1]))
                if not tup in T:
                    T[tup] = []
                    C[tup] = []
                C[tup].append(t)
                for pattern in P[t][s][diff][0]:
                    if not pattern in T[tup]:
                        T[tup].append(pattern)
        print('------------------------------------------------')
        print()

    '''for f in F:
        print('{}: {}'.format(f, F[f]))
    print()'''

    '''for t in T:
        print('{}: {} | {}'.format(t, T[t], C[t]))'''

    Maps = {}
    for t in T:
        Maps[t] = (T[t], C[t])
    '''for pattern in Maps:
        print('{}: {} | {}'.format(pattern, *Maps[pattern]))'''

    print()
    #I = (0,4,6,8,11)
    #I = (1,3,5,7,9)
    I = inv_r(5 * (4**a - 3**a), a)
    print('{} matches the patterns:'.format(I))
    Patterns = get_matching_patterns(I, Maps)
    for p in Patterns:
        c = pattern_to_comb(p)
        Js = apply_pattern(I, p, Maps)
        r = r_func_comb(I, c)
        print('  {} : {}'.format(p, Js))
        for J in Js:
            rJ = r_func(J)
            if r == rJ:
                print('    {} matches'.format(J))
            else:
                print('    {} DOES NOT match ({} vs. {})'.format(J, r, rJ))
