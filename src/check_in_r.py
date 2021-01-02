import time
import matplotlib.pyplot as plt

from math import log2, log
from test_r import *
from test import *

def v2(x):
    if x == 0:
        return 0
    i = 1
    while True:
        if x % 2**i != 0:
            return i - 1
        i += 1

def seq_diff(A):
    D = []
    i = 1
    while i < len(A):
        D.append(A[i] - A[i - 1])
        i += 1
    return D

def z(n):
    if n == 0:
        return 0
    if n == 1:
        return 1
    return int(-1/3 * (-1 + (-1)**(n)*(2)**(n)))

def psi(x):
    n = 0
    while True:
        #print((x + z(n)))
        if (x + z(n)) % (2**(n + 1)) == 0:
            return n
        n += 1

def v2(x):
    n = 1
    while x % 2**n == 0:
        n += 1
    return n - 1

def tuple_diff(A, B):
    return tuple([b - a for a, b in zip(A, B)])

def P(ds, de, p, r):
    return ds + r * (p, ) + de

def add_base(S):
    L = []
    i = 0
    while i < len(S):
        L.append(S[i] + i* 2 + 1)
        i += 1
    return tuple(L)


if __name__ == '__main__':
    '''a = 3
    P = []
    b = 2*a
    #print((a, b))
    r = 2**b - 3**a
    R = inv_r(r, a)
    #W = range(1, 1000)
    W = range(1, 3**a)
    for w in W:
        if in_Rn(w*r, a):
            if w % 2 == 1:
                A = inv_r(w*r, a)
                P.append((w, A))
                t = tuple_diff(R, A)
                #if w == 53:
                #print('{} | {}: {} - {}'.format(w, w * r, A, t))

    P = sorted(P, key = lambda x: x[1])
    for p in P:
        print(p)

    X = []
    Y = []
    Xp = []
    Yp = []

    Xt = []
    Yt = []

    i = 0
    while i < len(P):
        if P[i][0] == 1:
            for j in P[i][1]:
                Xp.append(i)
                Yp.append(j)
        else:
            for j in P[i][1]:
                X.append(i)
                Y.append(j)
        Xt.append(i)
        Yt.append(P[i][1][-1])
        i += 1

    plt.plot(X, Y, 'b.')
    plt.plot(Xp, Yp, 'r.')
    plt.plot(Xt, Yt)
    plt.show()'''

    print(r_func((1,3,5)))
    exit()

    B = range(1, 20)
    for b in B:
        A = range(1, b + 1)
        for a in A:
            #b = 2*a
            print((a, b))
            r = 2**b - 3**a
            R = inv_r(r, a)
            #W = range(1, 1000)
            W = range(1, min(3**a,1000))
            for w in W:
                if in_Rn(w*r, a):
                    if w % 2 == 1:
                        A = inv_r(w*r, a)
                        t = tuple_diff(R, A)
                        #if w == 53:
                        print('{} | {}: {} - {}'.format(w, w * r, A, w % 6))
            print('')


    '''for s in S:
        print('{}: {}'.format(s, S[s]))
    print('')

    for s in S:
        p = P(*S[s], 0)
        a = len(S[s][0]) + len(S[s][1]) + 0
        print(s, p, r_func(add_base(p)), r_func(add_base(p)) == s*(4**(a) - 3**a))'''




    '''D = []

    T = range(1, 250)
    for t in T:
        L = []
        b = 1
        while len(L) < 2:
            A = range(1, b + 1)
            for a in A:
                if in_Rn(2**b - t*3**a, a):
                    L.append((a, b))
                    #print('{}'.format((a, b)))
                    if len(L) >= 2:
                        break
            b += 1

        d1 = L[0][1] - 2 * L[0][0]
        d2 = L[1][1] - 2 * L[1][0]
        if d1 == d2:
            D.append((t, (2, d1), L[0][0]))

    for d in D:
        print(d)

    print(max([d[1][1] for d in D]))
    print(min([d[1][1] for d in D]))'''



    '''B = range(1, 9)
    for b in B:
        A = range(1, b + 1)
        for a in A:
            r = 2**b - 3**a
            W = range(1, 3**a)
            for w in W:
                if in_Rn(w*r, a):
                    print((a, b, w, r))
                    A = inv_r(w*r, a)
                    print('\t{}'.format(A))'''
