import time

from math import floor
from multiprocessing.pool import Pool
from functools import partial

import matplotlib.pyplot as plt
import numpy as np

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

def vp(x, p):
    if x == 0:
        return 0
    if x == 1:
        return 0
    n = 0
    while x % p == 0:
        x //= p
        n += 1
    return n

def collatz_accel(x):
    code = ''
    n = 0

    l = [x]
    while x != 1:
        n += 1
        if x % 2 == 1:
            x = (3*x + 1) // 2**(vp(3*x + 1, 2))
            code += '1'
        else:
            x = x // 2**(vp(x, 2))
            code += '0'
        l.append(x)

    res = {
    'iter': n,
    'code': code,
    'values': l
    }

    return res

def collatz(x):
    code = ''
    n = 0

    l = [x]
    while x != 1:
        n += 1
        if x % 2 == 1:
            x = (3*x + 1)
            code += '1'
        else:
            x = x // 2
            code += '0'
        l.append(x)

    res = {
    'iter': n,
    'code': code,
    'values': l
    }

    return res

#Utility function for multiprocessing self_cont
def _get_self_cont(bound, workers, offset, verbose=True):
    L = []
    n = 2 * offset + 1
    if offset == 1 and verbose:
        t0 = time.time()
        last = 1
        while n <= bound:
            if 100 * n / bound > last:
                last += 1

                print('{:.1f}% - {:.2f}s'.format(100 * n / bound, time.time() - t0))
                t0 = time.time()
            V = collatz_accel(n)['values'][1:]
            V = [v % n for v in V]
            if 0 in V:
                L.append(n)
            n += 2 * workers
    else:
        while n <= bound:
            V = collatz_accel(n)['values'][1:]
            V = [v % n for v in V]
            if 0 in V:
                L.append(n)
            n += 2 * workers
    return L

#Returns all self-contained numbers less than bound. Workers is number of workers to use in multiprocessing
def get_self_cont_multi(bound, workers=8):
    f = partial(_get_self_cont, bound, workers)

    params = list(range(workers))
    pool = Pool()
    res = pool.map(f, params)
    res = [r for s in res for r in s]

    R = [r for r in res]
    # add in even numbers
    for r in res:
        p = 1
        while True:
            V = collatz(r * 2**p)['values'][1:]
            V = [v % (r * 2**p) for v in V]
            if 0 not in V:
                break
            R.append(r*2**p)
            p += 1

    return sorted(R)

#gcd by euclidean algorithm
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

#modular inverse
def mod_inv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        return None
    else:
        return x % m

#Returns all nuumbers with some element equivalent to -3^-1 in their P orbit
def get_inv_P_orbit(bound, include_even=False):
    L = [1]
    n = 1
    step = 2
    if include_even:
        step = 1
    while n < bound:
        if n % 3 == 0:
            n += step
            continue
        V = collatz_accel(n)['values']
        V = [v % n for v in V]
        r = n - mod_inv(3, n)
        if r in V:
            L.append(n)
        n += step
    return L

#print same length arrays in unison
def print_arrs(A, B):
    sa = "["
    sb = "["
    for a, b in zip(A, B):
        l = max([len(str(a)), len(str(b))])
        sa += str(a).ljust(l)
        sa += ', '
        sb += str(b).ljust(l)
        sb += ', '
    sa = sa[:-1]
    sb = sb[:-1]
    sa += "]"
    sb += "]"
    print(sa)
    print(sb)

def isCSelfContained(n):
    d = n
    while n > 1:
        if n % 2 == 1:
            n = 3*n + 1
        else:
            n //= 2
        if n % d == 0:
            return True
    return False

#w is starting point
#depth is number of steps backwards
#height is highest power of 2 to multiply by
def backward_collatz(w, depth, height):
    L = [[w]]
    i = 0
    while i < depth:
        l = []
        for p in L[-1]:
            j = 1
            pref = p
            while j <= height:
                p = 2**j * pref - 1
                if p % 3 != 0:
                    j += 1
                    continue
                p //= 3
                l.append(p)
                j += 1
        L.append(l)
        i += 1
    return L


if __name__ == '__main__':

    #looks at obrits of powers of 2 multiplied by odd self-contianed numbers
    #Not necessarily useful, just something being played around with
    '''W = get_self_cont_multi(10**4)
    W = [w for w in W if w % 2 == 1]
    for w in W:
        V = collatz(w)['values'][1:]
        Vw = [v % w for v in V]
        k = V[Vw.index(0)] // w
        print('{} (k = {})'.format(w, k))
        P = range(0, 5)
        for p in P:
            wp = w * 2**p
            r = w - mod_inv(3, w)
            rp = wp - mod_inv(3, wp)
            V = collatz_accel(wp)['values']
            Vw = [v % w for v in V]
            Vwp = [v % wp for v in V]
            if rp in Vwp:
                print('2^{} ({}, {}) **'.format(p, r, rp))
            else:
                print('2^{} ({}, {})'.format(p, r, rp))
            #print(V)
            #print(Vw)
            #print(Vwp)
            print_arrs(Vw, Vwp)
            print('')
        print('--------')'''

    '''L = backward_collatz(379837, 16, 6)
    L_flat = [i for s in L for i in s]
    #for l in L:
        #print(l)

    print(max(L_flat))

    for w in L_flat:
        if isCSelfContained(w):
            print(w)'''

    W = [31, 62, 83, 166, 293, 347, 586, 671, 694, 1342, 2684, 19151, 38302, 2025797, 4051594]
    #W = [w for w in W if w % 2 == 1]

    for w in W:
        d = 3
        if w % 3 == 2:
            d = 6
        r = (4**1*w - d//3)//d
        O = collatz(w)['values'][1:]
        Ow = [o % w for o in O]
        v = O[Ow.index(0) - 1]
        p = (v - r) // w
        q = O[Ow.index(0)] // w
        print('{} ({}): {}'.format(w, q, p))

    '''W = range(1, 10**2)
    for w in W:
        if w % 3 == 0:
            continue
        d = 3
        if w % 3 == 2:
            d = 6
        O = collatz_accel(w)['values']
        Sw = [(2**1*w - d//3) // d]
        m = max(O)
        k = 2
        while Sw[-1] < m:
            Sw.append((2**k*w - d//3) // d)
            k += 1
        if len(set(Sw).intersection(set(O))) != 0:
            print('{}: {}'.format(w, Sw))
        print('{}: {}'.format(w, Sw))'''


    '''#Get self contained numbers up to bound, multi refers to multiprocessing
    #W = get_self_cont_multi(10**7, workers=8)
    W = [31, 62, 83, 166, 293, 347, 586, 671, 694, 1342, 2684, 19151, 38302, 2025797, 4051594]
    W = [w for w in W if w % 2 == 1]
    Worb = [collatz_accel(w)['values'] for w in W]
    for w, orb in zip(W, Worb):
        print('{}: {}'.format(w, orb))
        print('')
    #Get numbers with -3^-1 in their P orbit
    #P = get_inv_P_orbit(10**9, True)
    #print(P)'''

    '''W = [31, 62, 83, 166, 293, 347, 586, 671, 694, 1342, 2684, 19151, 38302, 2025797, 4051594]
    W = [w for w in W if w % 2 == 1]
    for w in W:
        r = w - mod_inv(3, w)
        print('{}: {}'.format(w, r))

    plt.plot(np.log(W))
    plt.show()'''
