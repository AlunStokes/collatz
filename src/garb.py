from math import sqrt, log, log2, ceil
from test_r import *
from scaling_factors import *
from utils import *
import itertools
from self_ref import get_self_cont_multi

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

def is_subset(A, B):
    A = set(A)
    B = set(B)
    for a in A:
        if a not in B:
            return False
    return True

def num_diff(a, b):
    if len(a) != len(b):
        raise Exception('Must be same size')
    n = 0
    for i, j in zip(a, b):
        if i != j:
            n += 1
    return n

#Returns the w for which w(2^b-3^a) in image(R_a)
def get_w_up_to(a, b, w_bound, stop_at_first=False):
    L = []
    W = range(1, w_bound + 1, 2)
    r = 2**b - 3**a
    for w in W:
        if in_Rn(w*r, a):
            L.append(w)
            if stop_at_first:
                break
    return L

def has_shared_element(A, B):
    for a in A:
        if a in B:
            return True
    return False

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def mod_inv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m

def all_equal(A, n):
    for a in A:
        if a != n:
            return False
    return True

if __name__ == '__main__':

    w = 10
    K = range(1, 12)
    S = []
    for k in K:
        d = 3
        if w % 3 == 2:
            d = 6
        S.append((4**k * w - d // 3) // d)
    Sw = [3 * s + 1 for s in S]
    r = w - mod_inv(3, w)
    Sr = [(s - r) // w for s in S]
    print(S)
    print(Sw)
    print(Sr)



    '''f = lambda x: 2*x + 4
    a = 7

    R = get_all_r(a, f(a))
    R = [r for r in R if r % 2 == 0]
    T = num_and_mult(R)
    A = [(inv_r(t[0], a), inv_r(t[1], a)) for t in T]
    for a in A:
        if num_diff(*a) == 1:
            print('{}: {}'.format(r_func(a[0]), a[0]))
            print('{}: {}'.format(r_func(a[1]), a[1]))
            print('')'''


    #KEEP ME
    #Check for stabilizatyion in linear pattern in b of a for mult factors
    '''f = lambda x : 2*x + 4

    N = 20
    D = {}

    A = range(1, N)
    for a in A:
        b = f(a)
        W = range(1, min(1000, 3**a), 2)
        L = []
        for w in W:
            r = w * (2**b - 3**a)
            if in_Rn(r, a):
                L.append(w)
                #print('{}: {}'.format(w, r))
            if len(L) > 0:
                D[(a,b)] = L

    for d in D:
        a = d[0] + 1
        b = f(a)
        if (a, b) in D:
            print('{}: {} - {}'.format(d, D[d], is_subset(D[d], D[(a, b)])))
        else:
            print('{}: {}'.format(d, D[d]))'''

    #check width for which w(2^b-3^a) in image(R_a) for b=2x+k has solutions as a grows
    '''P = []
    A = range(1, 20)
    for a in A:
        print(a)
        K = range(-a//2 - 1, a//2 + 2)
        L = []
        for k in K:
            f = lambda x: 2*x + k
            b = f(a)
            W = get_w_up_to(a, b, 100000, stop_at_first=True)
            if len(W) > 0:
                L.append(k)
        P.append(L)
    P_min = [min(p) for p in P]
    P_max = [max(p) for p in P]
    plt.plot(list(range(1, len(P_min) + 1)), P_min)
    plt.plot(list(range(1, len(P_max) + 1)), P_max)
    plt.show()'''


    '''#Something about self referential sequences
    D = []
    S = get_self_cont_multi(10**6)
    for s in S:
        L = []
        P = range(100)
        for p in P:
            w = s * 2**p
            if w % 3 == 0:
                continue
            V = collatz_phi(w)['values']
            d = 3
            if w % 3 == 2:
                d = 6
            S = []
            k = 1
            while True:
                r = (w * 4**k - d // 3) // d
                if r > max(V):
                    S.append(r)
                    break
                S.append(r)
                k += 1
            r = w - mod_inv(3, w)
            Vw = [v % w for v in V]
            if r in Vw:
                L.append(p)
                D.append(w)
                #print('{}'.format(V[Vw.index(r)]))
                #print('{}'.format(S))
        print('{}: {}'.format(s, L))
    print(D)
    print('')

    for w in D:
        V = collatz_phi(w)['values'][1:]
        Vm = [v % (w // 2**v2(w)) for v in V]
        d = 3
        if w % 3 == 2:
            d = 6
        r = w - mod_inv(3, w)
        print('{} ({}): {}'.format(w, r, Vm))
        print('')'''
