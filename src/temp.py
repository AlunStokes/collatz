from test_r import *
from utils import *

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    '''m = 7
    n = 1
    while n < 10:
        A = tuple(range(1, m*n, m))
        print(A)
        print(r_func(A))
        print(((2**m)**(n) - 3**n) // (2**m - 3))
        print('')
        n += 1'''

    '''plt.figure(num=None, figsize=(16, 6), dpi=140, facecolor='w', edgecolor='k')

    S = [4**n - 3**n for n in range(1, 15)]
    R = get_all_r(6, 18)
    R = [r for r in R if r % 2 == 1]
    #R = sorted(R)
    #print(sorted(R))
    M = num_and_mult(R)
    #print(M)
    #print(R)
    U = sorted(list(set([m[0] for m in M])))

    #U = [u for u in U if u % 2 == 1]
    #print(U)
    i = 0
    for u in R:
        if u not in U:
            continue
        plt.plot(R, 'g', markersize=1)
        T = '{}: ('.format(u)
        if u in S:
            T = '**' + T
        added = False
        for m in M:
            i1 = R.index(m[0])
            i2 = R.index(m[1])
            plt.plot([i1], [m[0]], 'r.', markersize=2)
            plt.plot([i2], [m[1]], 'k.', markersize=2)
            if i1 == R.index(u):
                if  (m[1] // u) % 2 == 0:
                    pass
                    #plt.plot([i1, i2], [m[0], m[1]], 'r')
                else:
                    added = True
                    plt.plot([i1, i2], [m[0], m[1]], 'r', linewidth=1)
                    T += '{}[{}], '.format(m[1], m[1] // u)
        if not added:
            plt.clf()
            continue
        T = T[:-2]
        T += ')'
        plt.title(T)
        plt.savefig('./images/graph_pointers/{}.png'.format(i))
        #plt.show()

        plt.clf()
        i += 1'''

    u = 3
    while u < 20:
        #print('{}^{}'.format(2, u))
        N = 0
        d = 3
        U = 2**u

        n = 1
        while n < U:
            if n % (U // 2**d) == 0:
                #print('{}/{}'.format(n // (U//2**d), 2**d))
                pass
            if n % 3 == 0:
                # print('')
                n += 1
                continue
            a = 1
            found = False
            while (3**a - 2**a) <= n:
                if in_Rn(n, a):
                    #print('{}: {}'.format(n, inv_r(n, a)))
                    N += 1
                    found = True
                    break
                a += 1
            #if not found:
                #print('{}: BAD'.format(n))
            n += 1
        print('{},{}'.format(2**u, N / U))
        u += 1

    '''a = 9
    b = 18

    R = get_all_r(a, b)
    M = num_and_mult(R)

    #print(M)

    for m in M:
        A1 = inv_r(m[0], a)
        A2 = inv_r(m[1], a)

        if A1[-1] >= A2[-1] and (m[1] // m[0]) % 2 == 1:
            print('{}: {}'.format(m, m[1] // m[0]))
            print('{} - {}'.format(A1, sum(A1)))
            print('{} - {}'.format(A2, sum(A2)))'''
