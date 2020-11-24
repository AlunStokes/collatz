from math import sqrt, log2, log
import operator as op
from functools import reduce
from itertools import combinations
from scipy.special import comb
from utils import *

import matplotlib.pyplot as plt

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom

#smallest n such that d divides (an + b)
def smallest_n(d, a, b):
    if b == 0:
        return 0
    n = 0
    while True:
        if (a*n + b) > 0:
            if (a*n + b) % d == 0:
                return n
        n += 1

def tern(n, pad_to=None):
    if n == 0:
        s = '0'
        if pad_to:
            s = '0' * (pad_to - len(s)) + s
        return s
    nums = []
    while n:
        n, r = divmod(n, 3)
        nums.append(str(r))
    s = ''.join(reversed(nums))
    if pad_to:
        s = '0' * (pad_to - len(s)) + s
    return s

#n is length of ternary strings
def get_nondec_ternary(a, b):
    L = []
    num = ncr(b, a)
    bound = int(sqrt(num)) + 1

    n = 0
    while n < bound:
        m = 0
        while m <= n:
            L.append((3**m + 3**n) // 2)
            m += 1
        n += 1

    L = [tern(l - 1, a) for l in L[:num]]
    return L

#a is |A|, b is |b|
def get_all_r(a ,b, include_even=True):
    '''L = []
    M = [list(range(i + 1, i + 1 + a)) for i in range(b - a + 1)]
    print(M)
    T = get_nondec_ternary(a, b)
    T = [t[::-1] for t in T]
    for t in T:
        r = 0
        for i, c in enumerate(t):
            r += 3**i * 2**(M[len(M)-1-i][int(c)] - 1)
        L.append(r)
    print(L)'''

    L = []
    M = list(combinations(list(range(1, b + 1)), a))
    if not include_even:
        M = [m for m in M if m[0] == 1]
    for m in M:
        A = []
        r = 0
        for i, p in enumerate(m[::-1]):
            A.append(p)
            r += 3**(i) * 2**(p - 1)
        #L.append((tuple(A[::-1]), r))
        L.append(r)
    return L

def is_prime(n):
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    i = 3
    while i < int(sqrt(n)) + 1:
        if n % i == 0:
            return False
        i += 2
    return True

def phi(x):
    if x == 0:
        return 0
    i = 1
    while True:
        if x % 2**i != 0:
            return i - 1
        i += 1

def r_func(A):
    r = 0
    i = 1
    while i <= len(A):
        r += 3**(i-1) * 2**(A[len(A) - i] - 1)
        i += 1
    return r

def r_func_seq(A):
    r = []
    i = 1
    while i <= len(A):
        r.append(3**(i-1) * 2**(A[len(A) - i] - 1))
        i += 1
    return r[::-1]

def inv_r(r, n):
    L = []
    i = 0
    while i < n:
        s = 0
        j = 0
        while j < i:
            s += 3**(n-(j+1)) * 2**(L[j] - 1)
            j += 1
        L.append(phi(r - s) + 1)
        i += 1
    return tuple(L)

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

def check_hyp(R, n):
    for r in R:
        for i in range(n-1, 0, -1):
            r = r // 2**phi(r)
            r -= 3**i

        if 2**phi(r) != r:
            return False
    return True

def in_Rn(r, n):
    if n == 0 and r == 0:
        return True
    elif n == 0 and r != 0:
        return False
    for i in range(n-1, 0, -1):
        r = r // 2**phi(r)
        r -= 3**i

    if 2**phi(r) != r:
        return False
    return True

def recurse_image(R, n):
    P = [prime_fac(r) for r in R]
    WR = [r for r, p in zip(R,P) if 2 not in p]
    WRp = []
    for r in WR:
        k = 1
        while k < 7:
            WRp.append(2**k * r + 3**n)
            k += 1
    return sorted(list(set(WRp)))

def is_power_of(x, n):
    if n < 0:
        n = -n
    i = 0
    while n**i <= x:
        if n**i == x:
            return True
        i += 1
    return False

#Checks if number of form 2^b - 3^a
def of_form(x):
    A = range(0, 20)
    for a in A:
        if is_power_of(x + 3**a, 2):
            return True
        if is_power_of(x + 2**a, 3):
            return True
    return False

def n_odd(S):
    n = 0
    for s in S:
        if s % 2 == 1:
            n += 1
    return n

def n_even(S):
    n = 0
    for s in S:
        if s % 2 == 0:
            n += 1
    return n

def is_subset(A, B):
    for a in A:
        if a not in B:
            return False
    for a in A:
        if B.count(a) < A.count(a):
            return False
    return True

def divisor_in(r, P):
    for p in P:
        if p == r:
            continue
        if r % p == 0:
            return (True, p)
    return (False, -1)

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

def merge_lists(P):
    Q = []
    for p in P:
        Q += p
    return Q

def num_and_mult(P):
    S = []
    for p in P:
        for q in P:
            if p == q:
                continue
            if q % p == 0:
                S.append((p, q))
    return S



if __name__ == '__main__':

    a = 5
    b = 2*a + 8

    w = 5
    #R = sorted(get_all_r(a, b))
    R = get_all_r(a, b)
    A = [inv_r(r,a) for r in R]
    #P = [prime_fac(r) for r in R]

    ref_A = tuple(range(1, 2*a, 2))
    ref_r = r_func(ref_A)
    ref_fac = prime_fac(ref_r)

    Q = [[r for a, r in zip(A, R) if max(a) == i] for i in range(a, b + 1)]

    for i, q in enumerate(Q):
        print('{} [{} - {}]'.format(i + a, 2**(i  + a - 1) + 3 * (3**(a - 1) - 2**(a - 1)), 2**(i) * (3**a - 2**a) ))
        M = []
        for p in q:
            if p % 2 == 1:
                M.append(p)
                if p % ref_r == 0:
                    A = inv_r(p, a)
                    print('{}: {} | {} - {}'.format(p, p // ref_r, A, len(collatz_accel(p // ref_r)['code'])))
        #print(num_and_mult(M))
        M = [M[i] - M[0] for i in range(1, len(M))]
        #print(M)
        #print(sorted([(m // 6) for m in M]))
        #print([prime_fac(m) for m in M])
        print('')
    print(ref_A)
    print(ref_r)
    print(ref_fac)
    N = [M[i] for i in range(1, len(M)) if M[i] < M[i - 1]]

    #plt.plot(M)
    #plt.show()
    #print(N)
    #print([2 * (n + 30) for n in N])



    '''a = 4
    b = 9
    R = sorted(get_all_r(a, b))
    A = [inv_r(r, a) for r in R]

    P = []
    for i in range(a, b):
        P.append([(r, a) for r, a in zip(R, A) if a[-1] == i])

    print('2^b - 3^a = {}'.format(2**b - 3**a))

    for i, p in enumerate(P):
        print(i+a, '({})'.format(len([x for x in p if x[0] % 2 == 1])))
        for j in p:
            #if 2 not in prime_fac(j[0]):
            print('\t{} : {}'.format(j, prime_fac(j[0])))'''
