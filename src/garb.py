from math import sqrt, log, log2
from test_r import *
from scaling_factors import *
from utils import *
import itertools

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

def smallest_n(d, a, b):
    if b == 0:
        return 0
    n = 0
    while True:
        if (a*n + b) > 0:
            if (a*n + b) % d == 0:
                 return n
        n += 1

def phi(n):

    result = 1
    for i in range(2, n):
        if (gcd(i, n) == 1):
            result+=1
    return result

def mod_inv(a, m):
    return a**(phi(m)-1) % m

def check_r(r, a):
    A = []
    i = 0
    while i < a:
        s = []
        j = a - (i + 1)
        while 3**i * 2**j < r:
            s.append(3**i * 2**j)
            j += 1
        A.append(s)
        i += 1

    print(A)
    L = []

    i = 0
    while i < len(A[0]):
        j = 0
        while j < len(A[1]):
            k = 0
            while k < len(A[2]):
                if A[0][i] + A[1][j] + A[2][k] == r:
                    L.append((A[0][i], A[1][j], A[2][k]))
                k += 1
            j += 1
        i += 1
    return L

def check_bounds(a):
    lower = a * log(3) / log(2)
    upper = (a * log(2) + 2*a * log(3) - log(3**a*2**a - (3**a - 2**a))) / log(2)

    d = upper - lower

    return (lower, upper, d)

def is_subset(A, B):
    for a in A:
        if a not in B:
            return False
    for a in A:
        if A.count(a) > B.count(a):
            return False
    return True

def num_and_mult(P):
    S = []
    for p in P:
        for q in P:
            if p == q:
                continue
            if q % p == 0:
                S.append((p, q))
    return S

def all_equal(L):
    i = 0
    while i < len(L) - 1:
        if L[i] != L[i+1]:
            return False
        i += 1
    return True

def find_conv(w_bound, n_rep=4):
    L = []
    a = 2
    while len(L) < n_rep or not all_equal(L[-n_rep:]):
        L.append(len(get_wa(a, w_bound)))
        a += 1
    return a - n_rep

def r_func_mat(X, Y):
    Z = []
    for x, y in zip(X, Y):
        Z.append(r_func((x, y)))
    return np.array(Z)

def has_multiple(S):
    for s in S:
        if S.count(s) > 1:
            return True
    return False

def get_gaps(S):
    S = sorted(S)
    L = []
    i = 1
    while i < len(S):
        L.append(S[i] - S[i - 1])
        i += 1
    return L

def all_squares_up_to(U):
    L = []
    n = 1
    while n**2 <= U:
        L.append(2**n)
        n += 1
    return L

if __name__ == '__main__':
    u = 2
    while u < 60 :
        U = 2**u
        S = all_squares_up_to(U)
        L = []
        for s in S:
            n = 0
            while s * 2**n <= U:
                L.append((s + 3)*2**n)
                n += 1
        #print(L)
        r = len(L) / U
        c = (sqrt(2**(1-sqrt(U)) * U)-2*sqrt(U)) / ((sqrt(2) - 2) * U)
        print('{} ({}): {}, {} - {}'.format(U, len(L), r, c, r < c))
        u += 1
