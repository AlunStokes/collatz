#from math import log, ceil
from mpmath import *
from utils import gcd

import matplotlib.pyplot as plt

mp.dps = 100
mp.pretty = True

def lower_bound(a):
    return a * ln(3) / ln(2)

def upper_bound(a):
    return (a * ln(2) + 2 * a * ln(3) - ln(3**a * 2**a - 3**a + 2**a)) / ln(2)

def log_ceil_dist(a):
    return ceil(a * ln(3) / ln(2)) - a * ln(3) / ln(2)

def log_floor_dist(a):
    return a * ln(3) / ln(2) - floor(a * ln(3) / ln(2))

def log_dist(a):
    return min([log_ceil_dist(a), log_floor_dist(a)])

def graph_decreasing_log_ceil_dist(A):
    L = []
    for a in A:
        if len(L) == 0:
            L.append((a, log_ceil_dist(a)))
            continue
        r = log_ceil_dist(a)
        if r <= L[-1][1]:
            L.append((a,r))

    X = [l[0] for l in L]
    Y = [l[1] for l in L]

    #print(L)

    plt.plot(X, Y, 'r')
    #plt.show()

def graph_decreasing_interval_width(A):
    L = []
    for a in A:
        L.append((a, upper_bound(a) - lower_bound(a)))

    X = [l[0] for l in L]
    Y = [l[1] for l in L]

    plt.plot(X, Y, 'b')

def graph_denom_order(A):
    L = []
    for a in A:
        if len(L) == 0:
            L.append((a, log_ceil_dist(a)))
            continue
        r = log_ceil_dist(a)
        if r <= L[-1][1]:
            L.append((a,r))

    X = [log(l[0]) for l in L]
    Y = [log(l[1]) for l in L]

    print(L[-1])

    plt.plot(X,Y, 'r')

def graph_reference(A):
    X = [log(a) for a in A]
    Y = [-log(a)**(4/3) for a in A]
    plt.plot(X,Y, 'b')

def reduce_frac(num, den):
    while gcd(num, den) != 1:
        i = 2
        while True:
            if num % i == 0 and den % i == 0:
                num //= i
                den //= i
                break
            i += 1
    return num, den

def graph_order_approx(acc=30):
    temp = mp.dps
    mp.dps = acc
    X = []
    Y = []
    r = log(3) / log(2)
    rs = str(r).replace('.', '')
    print(rs)
    i = 1
    while i < len(rs):
        num = int(rs[:i])
        den = 10**(i-1)

        num, den = reduce_frac(num, den)
        v = num / den
        print('{}/{}'.format(num, den))
        X.append(den)
        Y.append(abs(r - v))
        i += 1

    X = [log(x) for x in X]
    Y = [log(y) for y in Y]

    #print(X)
    #print(Y)

    plt.plot(X,Y)
    plt.show()

    mp.dps = temp

if __name__ == '__main__':
    '''A = range(1, 500000)
    #graph_decreasing_log_ceil_dist(A)
    #graph_decreasing_interval_width(A)
    graph_denom_order(A)
    graph_reference(A)
    plt.show()'''

    A = range(1, 500000)
    graph_order_approx(120)
