import time

from math import floor
from multiprocessing.pool import Pool
from functools import partial
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

from utils import *
from garb import nu
from scaling_factors import *
from test_r import *
from self_ref import get_self_cont_multi, get_pseudo_self_cont_multi, get_neg_inv_multi

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

#ensure sequential makes sure that this is also valid for a, a + 1, ..., a + (ensure_sequential)
#if ensure_sequential == 0, we don't check.
def get_wa(a, w_bound=False, ensure_sequential=0, include_even=False, verbose=False):
    wab = []
    b = 2*a
    if verbose:
        print((a, b))

    r = 2**b - 3**a
    skip = 2
    if include_even:
        skip = 1
    W = range(1, 3**a, skip)
    if w_bound:
        W = range(1, min(3**a, w_bound), skip)
    for w in W:
        if in_Rn(w*r, a):
            if verbose:
                print('{} | {}: {}'.format(w, w * r, inv_r(w*r, a)))
            in_all = True
            for k in range(ensure_sequential):
                bp = 2*(a + k + 1)
                rp = 2**bp - 3**(a + k + 1)
                if not in_Rn(w*rp, a + k + 1):
                    in_all = False
                    break
            if in_all:
                wab.append(w)

    return wab

def seq_diff(a):
    b = [a[0]] + [a[i] - a[i - 1] for i in range(1, len(a))]
    return tuple(b)

def most_common(lst):
    return max(set(lst), key=lst.count)

def index_of_longest_seq_of(L, s):
    longest_seq = -1
    longest_index = -1

    for i, l in enumerate(L):
        if l == s:
            c = 1
            for k in range(1, len(L) - i):
                if L[i + k] == s:
                    c += 1
                else:
                    break
            if c > longest_seq:
                longest_seq = c
                longest_index = i

    return longest_index, longest_seq


if __name__ == '__main__':

    #Check how to decompose into r_1 + 3^t*2^(s - 2*t)
    '''num_rep = 5

    w = 85
    l, s, t = 1,6,3
    #w = 53
    #l, s, t = 2,5,2
    #w = 4835
    #l, s, t = 3,11,5
    #w = 4277
    #l, s, t = 5,10,5

    #w = 1099
    #l, s, t = 8, 7, 2
    #l, s, t = 8, 7, 3

    a = l + t + num_rep
    I = tuple([i - 1 for i in inv_r(w * (2**(2*a) - 3**a), a)])
    Idiff = [i - 2*j for j, i in enumerate(I)]
    leading = tuple(Idiff[:l])
    trailing = tuple(Idiff[l + num_rep:])
    r1 = w - 3**t * 2**(s - 2*t)
    print('w = {}'.format(w))
    print('For a = {}: {} | {}'.format(a, I, Idiff))
    print('(l, s, t) = {}'.format((l, s, t)))
    print('{}, {}, {}'.format(leading, s, trailing))
    print('r_1 = {}'.format(3**t * 2**(s - 2*t), r1))
    Ir1 = tuple([i - 1 for i in inv_r(r1, t)])
    Ir1p = tuple([s - 2*t,] + list(Ir1))
    print('I_(t)(r_1) = {}'.format(Ir1))
    print('I_(t+1)(r_1\') = {}'.format(Ir1p))
    if in_Rn(w, t + 1):
        Iw = tuple([i - 1 for i in inv_r(w, t+1)])
        print('I_(t+1)(w) = {}'.format(Iw))
    else:
        print('w does NOT have a (t+1)-special representation')
    print()
    exit()'''

    '''w = 341
    l = 1
    s = 8
    t = 4
    r1 = w - 3**t * 2**(s - 2 * t)
    r2 = 2**(s + 2*l) - 3**l * w
    a = 9
    print(tuple_diff([i - 1 for i in inv_r(w*(4**a - 3**a), a)], [2*i for i in range(a)]))
    print('{}'.format(r1 + 3**t*2**(s - 2*t)))
    Ir1 = tuple([i - 1 for i in inv_r(r1, t)])

    print(in_Rn(r1, t))
    print(in_Rn(r2, l))
    if in_Rn(r1, t):
        print('{}: {}'.format(r1, tuple([i - 1 for i in inv_r(r1, t)])))
    if in_Rn(r2, l):
        print('{}: {}'.format(r2, tuple([i - 1 for i in inv_r(r2, l)])))
    exit()'''

    '''a = 20
    W = get_wa(a, 100000)[1:]
    for w in W:
        wa = a
        r = w * (4**wa - 3**wa)
        while in_Rn(r, wa):
            wa -= 1
            r = w * (4**wa - 3**wa)
        wa += 1
        r = w * (4**wa - 3**wa)
        I = tuple([i - 1 for i in inv_r(r, wa)])
        L = []
        T = []
        for i in I:
            if i < 2*wa:
                L.append(i)
            else:
                T.append(i - 2*wa)
        r1 = r_func([i + 1 for i in T])
        r2 = r_func([i + 1 for i in L])
        print('{} (a = {}): ({}) {} {} | ({}) {} {}'.format(w, wa, r2, False, L, r1, is_power_of_p((w -r1)//3**(len(T)), 2), T))
    exit()'''


    '''a = 7
    b = 2*a
    m = 5
    W = get_wa(a, 1000)[1:]
    for w in W:
        S = []
        s = 1
        while s < a:
            if in_Rn(w, s):
                S.append(s)
            s += 1
        print('w  = {}'.format(w))
        print('  has a-special rep.:')
        for s in S:
            Is = tuple([i - 1 for i in inv_r(w, s)])
            print('    {}: {}'.format(s, Is))

        r = w * (2**b - 3**a)
        I = tuple([i - 1 for i in inv_r(r, a)])
        Id = tuple([i - j for i, j in zip(I, [2*i for i in range(a)])])
        L = []
        T = []
        for i in I:
            if i < 2*a:
                L.append(i)
            else:
                T.append(i - 2*a)
        L = tuple([i + 1 for i in L])
        T = tuple([i + 1 for i in T])
        rL = r_func(L)
        rT = r_func(T)
        print('{}: {}'.format(a, I))
        print('{}: {}'.format(' ' * len(str(a)), Id))
        print('{}: {} | {}'.format(' ' * len(str(a)), prime_fac(rL), prime_fac(rT)))
        k = 1
        while k <= 5:
            ak = a + k
            bk = 2 * ak
            rk = w*(2**bk - 3**ak)
            if in_Rn(rk, ak):
                Ik = tuple([i - 1 for i in inv_r(rk, ak)])
                Ikd = tuple([i - j for i, j in zip(Ik, [2*i for i in range(ak)])])
                print('{}: {}'.format(ak, Ik))
                print('{}: {}'.format(' ' * len(str(ak)), Ikd))
            else:
                print('{}: ------------------'.format(ak))
                print('{}: ------------------'.format(' ' * len(str(ak))))
            k += 1
        print()
    exit()'''



    '''w = 821
    a = 7
    b = 2*a
    r = w * (2**b - 3**a)
    I = [i - 1 for i in inv_r(r, a)]
    I1 = [2*i for i in range(a)]
    D = tuple_diff(I, I1)
    s = 0
    for j, i in enumerate(I[::-1]):
        s += 3**j * 2**i
    print(I)
    print(D)
    print(s)
    print(r)
    exit()'''

    '''for w in range(821, 1000, 2):
        Sol1 = []
        Sol2 = []
        for t in range(0,60):
            for s in range(-20, 20):
                r = w - 3**t * 2**(s - 2*t)
                if in_Rn(r, t):
                    #print('(_, {}, {}): {}'.format(s, t, [i - 1 for i in inv_r(r, t)]))
                    Sol1.append((-1, s, t))

        for l in range(0,60):
            for s in range(-20, 20):
                r = 2**(s + 2*l) - 3**l * w
                if in_Rn(r, l):
                    #print('({}, {}, _): {}'.format(l, s, [i - 1 for i in inv_r(r, l)]))
                    Sol2.append((l, s, -1))

        S1 = set([x[1] for x in Sol1])
        S2 = set([x[1] for x in Sol2])

        S = S1.intersection(S2)
        if len(S) == 0:
            continue
        Sol = []
        for s in S:
            for s1 in Sol1:
                for s2 in Sol2:
                    if s1[1] == s and s2[1] == s:
                        sol = (s2[0], s, s1[2])
                        Sol.append(sol)


        min_a = -1
        for a in range(1, 100):
            if in_Rn(w * (4**a - 3**a), a):
                min_a = a
                break

        if min_a == -1:
            continue

        print('{}: {}'.format(w, Sol[:5]))
        sol = Sol[0]
        l, s, t = sol
        r1 = w - 3**t * 2**(s - 2*t)
        r2 = 2**(s + 2*l) - 3**l * w
        print('V_w: {}'.format(tuple([i - 1 for i in inv_r(w*(4**min_a - 3**min_a), min_a)])))
        print('V_r1: {}'.format(tuple([i - 1 for i in inv_r(r1, t)])))
        print('V_r2: {}'.format(tuple([i - 1 for i in inv_r(r2, l)])))
        print('Minimum a = {}'.format(a))
        print('')

    exit()'''

    #KEEP ME
    #find smallest m so that p_2 + 3^m*w divides 4^a
    w = 5
    a = 4
    for m in range(a+1):
        for p2 in get_all_r(m, 2*a - m + 15):
        #for p2 in get_all_r(m, 2*a - m + 2):
            r = p2 + 3**m * w
            if r % 4**a == 0:
                I = tuple([i - 1 for i in inv_r(p2, m)])
                print('m = {}'.format(m))
                print('r = {} + 3^{} * w = {}'.format(p2, m, r))
                print('R^(-1)(p2) = {}'.format(I))
                print('r = {}'.format(prime_fac(r)))
                print()
    exit()

    #Check min m so that 4^a divides p2 + 3^m * w
    '''a = 7
    W = get_wab(a, 2*a, 10000)
    #W = range(1, 100, 2)
    #W = get_self_cont_multi(10**5)
    for w in W:
        found = False
        M = {}
        for m in range(a + 1):
            c = 0
            for p2 in get_all_r(m, 2*a - m + 12):
            #for p2 in get_all_r(m, 2*a - m + 2):
                r = p2 + 3**m * w
                if r % 4**a == 0:
                    c += 1
                    #print((w,m))
                    #found = True
                    #break
            if found:
                break
            if c > 0:
                M[m] = c
        print('{}: {}'.format(w, M))


    exit()'''


    #Check how prime fac of r + 3^a changes
    '''w = 821
    a = 7
    m = 5
    b = 2*a + 8
    R = get_all_r(m, b)
    R = [r for r in R if r % 2 == 1]
    Is = [tuple([i - 1 for i in inv_r(r, m)]) for r in R]
    for r, I in zip (R, Is):
        rp = r + 3**m * w
        rp_fac = prime_fac(rp)
        r_fac = prime_fac(r)
        if rp_fac.count(2) >= 2*a:
            print('{} ({}) = {} -> {} = {} ({})'.format(r, I, r_fac, rp, rp_fac, rp_fac.count(2)))
    exit()'''

    ###KEEP ME
    #Find solutions corresponding to minimal m for (I)
    a = 8
    W = get_wab(a, 2*a, 10000)
    for w in W:
        print('w = {} ({})'.format(w, prime_fac(w)))
        found = False
        for m in range(a + 1):
            for p2 in get_all_r(m, 2*a - m + 12):
                r = p2 + 3**m * w
                if r % 4**a == 0:
                    r_fac = prime_fac(r)
                    I = tuple([i - 1 for i in inv_r(p2, m)])
                    print('m = {}'.format(m))
                    print('r = {} + 3^{} * w = {}'.format(p2, m, r))
                    print('R^(-1)(p2) = {}'.format(I))
                    print('r = {} ({})'.format(r_fac, r_fac.count(2)))
                    print()

                    found = True
                    break
            if found:
                break
    exit()

    #Check num of m so that 3^(m - a) divides w - p1
    '''a = 7
    W = get_wab(a, 2*a, 10000)
    #W = range(1, 100, 2)
    #W = get_self_cont_multi(10**5)
    for w in W:
        found = False
        M = {}
        for m in range(a + 1):
            c = 0
            for p1 in get_all_r(a - m, 12):
            #for p2 in get_all_r(m, 2*a - m + 2):
                r = w - p1
                if r % 3**(a - m) == 0:
                    c += 1
                    #print((w,m))
                    #found = True
                    #break
            if found:
                break
            if c > 0:
                M[m] = c
        print('{}: {}'.format(w, M))
    exit()'''

    w = 821
    a = 7
    b = 2*a + 10
    for m in range(1, a + 1):
        R = get_all_r(m, b)
        R = [r for r in R if r % 2 == 1]
        Is = [tuple([i - 1 for i in inv_r(r, m)]) for r in R]
        R, Is = (list(t) for t in zip(*sorted(zip(R, Is), key = lambda x: nu(x[0] + 3**m * w, 2))))
        for r, I in zip(R, Is):
            rp = r + 3**m * w
            #print('{} ({}) = {} -> {} = {}'.format(r, I, prime_fac(r), rp, prime_fac(rp)))
            if nu(r + 3**m * w) >= 2*a:
                print('{} ({}) = {} -> {} = {}'.format(r, I, prime_fac(r), rp, prime_fac(rp)))
        '''for r, I in zip(R, Is):
            rp = r + 3**a * w
            print('{} ({}) = {} -> {} = {}'.format(r, I, prime_fac(r), rp, prime_fac(rp)))'''

    exit()

    #KEEP THESE TWO
    #Run on one
    #print(tuple([i - 1 for i in inv_r(85*(4**7 - 3**7), 7)]))
    w = 325
    a = 7
    for m in range(0, a + 1):
        print('m = {}'.format(m))
        for p1 in get_all_r(a - m, 11):
            if p1 == 0:
                p2bound = 2*a + 10
            else:
                minp1 = min([i - 1 for i in inv_r(p1, a - m)])
                p2bound = 2*a + minp1 + 1
            for p2 in get_all_r(m, p2bound):
                lhs = 3**(a - m) * (p2 + 3**m*w)
                rhs = 4**a * (w - p1)
                if lhs == rhs:
                    Ip1 = tuple([i - 1 for i in inv_r(p1, a-m)])
                    Ip2 = tuple([i - 1 for i in inv_r(p2, m)])
                    print('w = {}: ({}, {}, {})'.format(w, m, Ip1, Ip2))
                    print('  p1 = {}'.format(prime_fac(p1)))
                    print('  p2 = {}'.format(prime_fac(p2)))
                    print('  w - p1 = {}'.format(prime_fac(w - p1)))
                    print('  p2 + 3^m*w = {}'.format(prime_fac(p2 + 3**m  * w)))

    exit()

    #run over a bunch
    a = 8
    U = 100
    W_known = get_wab(a, 2*a, U)
    W_guess = []

    print(W_known)

    #w = 83
    #for w in range(1, U, 2):
    for w in W_known:
        found = False
        for m in range(0, a + 1):
            #print('m = {}'.format(m))
            for p1 in get_all_r(a - m, 11):
                if p1 == 0:
                    p2bound = 2*a + 10
                else:
                    minp1 = min([i - 1 for i in inv_r(p1, a - m)])
                    p2bound = 2*a + minp1 + 1
                for p2 in get_all_r(m, p2bound):
                    lhs = 3**(a - m) * (p2 + 3**m*w)
                    rhs = 4**a * (w - p1)
                    if lhs == rhs:
                        Ip1 = tuple([i - 1 for i in inv_r(p1, a-m)])
                        Ip2 = tuple([i - 1 for i in inv_r(p2, m)])
                        W_guess.append((w, (m, Ip1, Ip2)))
                        print('w = {}: ({}, {}, {})'.format(w, m, Ip1, Ip2))
                        found = True
                        '''print((m, (p1, p2)))
                        Ip1 = tuple([i - 1 for i in inv_r(p1, a-m)])
                        Ip2 = tuple([i - 1 for i in inv_r(p2, m)])
                        print('  p1: {} = {}'.format(p1, Ip1))
                        print('  p2: {} = {}'.format(p2, Ip2))'''
                        break
                if found:
                    pass
                    break
            if found:
                pass
                break
    print(W_guess)
