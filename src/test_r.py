from math import sqrt, log2, log
import operator as op
from functools import reduce
from itertools import combinations
from itertools import permutations
from itertools import product
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
    if n == 0:
        if r == 0:
            return tuple()
        else:
            return (-1,)
    if n == 1:
        v = vp(r, 2)
        if 2**v == r:
            return (v+1,)
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
    if int(r) != r:
        return False
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


def factor(n):
    L = []
    i = 2
    while i < int(np.sqrt(n)) + 1:
        if n % i == 0:
            L.append(i)
            L.append(n // i)
        i += 1
    #if len(L) == 0:
        #L.append(n)
    L = list(sorted(L))
    return L

def primes_less_than(N):
    L = [2]
    n = 3
    while n < N:
        if is_prime(n):
            L.append(n)
        n += 1
    return L

def f(a, b):
    return 2**b - 3**a


if __name__ == '__main__':

    w = 4835
    L = list(range(12))
    Ls = [L for _ in range(6)]
    perms = product(*Ls)
    for perm in perms:
        s = 0
        for i in range(len(perm)):
            s += 3**(len(perm) - i -1) * 2**perm[i]
        if s == w:
            print(perm)
    exit()

    #Old guyg who goes through combiuations for which w works but largest element is too big
    #***
    '''a = 4
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
    #print([2 * (n + 30) for n in N])'''

    #Check what happens upon multiplication of numbers in image R_a
    '''n = 4
    R = get_all_r(n, 8)
    R = list(sorted(R))
    print(R)
    T = []
    for r in R:
        for q in R:
            if in_Rn(r * q, 2*n - 1):
                if (r,q) not in T and (q,r) not in T:
                    T.append((r, q))
    for t in T:
        print('{}: {} * {} = {}'.format(t, t[0], t[1], t[0] * t[1]))
        print('\t{}: {}'.format(t[0], inv_r(t[0], n)))
        print('\t{}: {}'.format(t[1], inv_r(t[1], n)))
        print('\t{}: {}'.format(t[0] * t[1], inv_r(t[0] * t[1], 2*n - 1)))
        print('')'''

    #Check if all numbers in image(R_a) can be factored into other such numbers
    '''n = 3
    R = get_all_r(n, 10)
    R = list(sorted(R))
    print(R)
    T = []
    for r in R:
        I = inv_r(r, n)
        F = factor(r)
        F = F[:len(F)//2]
        fac = False
        for f in F:
            p = r // f
            q = f
            i = 1
            p_fac = False
            while i < n:
                if in_Rn(p, i):
                    p_fac = True
                    p_I = inv_r(p, i)
                    break
                i += 1
            i = 1
            q_fac = False
            while i < n:
                if in_Rn(q, i):
                    q_fac = True
                    q_I = inv_r(q, i)
                    break
                i += 1
            if p_fac and q_fac:
                fac = True
                break
        if fac:
            print('{} can be factored into {} = {} * {}'.format(r, I, p_I, q_I))
        else:
            print('{} is not factorable'.format(r))'''


    #checking some numbers in image(R_a)
    M = 10
    #a = 5
    #b = 10
    T = []
    B = range(1, M + 1)
    for b in B:
        A = range(1, b + 1)
        for a in A:
            W = range(1, 100)
            for w in W:
                r = w * (2**b - 3**a)
                if r > 0:
                    T.append(((a, b, w), r))

    T.sort(key=lambda x:x[1])
    #print(T)
    for t in T:
        if in_Rn(t[1], t[0][0]):
            print('{} {} in R_{}'.format(t[1], t[0], t[0][0]))

    #Check which (a,b) allow scaling
    '''U = 100
    D = {}
    for b in range(1, U):
        print(b)
        for a in range(1, b + 1):
            D[(a, b)] = []
            for w in range(1, 100, 2):
                if in_Rn(w * (2**b - 3**a), a):
                    D[(a, b)].append(w)

    for d in D:
        if len(D[d]) > 0:
            print('{}: {}'.format(d, D[d]))'''

    #checking for specific linear relation
    '''p = 0
    U = 30
    D = {}
    for a in range(1, U):
        b = 2 * a + p
        D[(a, b)] = []
        for w in range(1, 2000, 2):
            if in_Rn(w * (2**b - 3**a), a):
                D[(a, b)].append(w)

    for d in D:
        if len(D[d]) > 0:
            print('{}: {}'.format(d, D[d]))

    plt.plot(D[(U - 1, 2 * (U - 1) + p)])
    plt.show()'''
