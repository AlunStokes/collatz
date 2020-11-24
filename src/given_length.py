import time

from sympy import *
from sympy.utilities.lambdify import lambdify, implemented_function

from utils import *

init_printing(use_unicode=True)

def collatz(x, code):
    #code = code[::-1]
    for c in code:
        if c == '0':
            x = x / 2
        else:
            x = (3*x + 1) / 2
    return x

def tokenize_exp(exp):
    exp = str(exp)
    d = int(exp.split('/')[-1])
    if '*' in exp:
        m = int(exp[exp.index('(') + 1 : exp.index('*')])
        n = int(exp[exp.index('+') + 2: exp.index(')')])
    else:
        m = 1
        n = 0

    return m, n, d

def tokenize_relation(rel):
    rel = str(rel)
    if len(rel) == 1:
        return 1, 0

    if '*' in rel:
        a = int(rel[:rel.index('*')])
    else:
        a = 1
    if '+' in rel or '-' in rel:
        if '+' in rel:
            b = int(rel[rel.index('+') + 2:])
        else:
            b = -int(rel[rel.index('-') + 2:])
    else:
        b = 0

    return a,b

def set_linear_relation(a, b, n):
    k = 0
    l = []
    while a * k + b <= n:
        l.append(a * k + b)
        k += 1
    return l

def all_bin_len(n):
    return [bin(i)[2:].rjust(n, '0') for i in range(2**n)]

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

def get_const_r(b):#
    #b = b[::-1]
    m = 0
    A = []
    for j, i in enumerate(b):
        if i == '1':
            #Extra for 1-indexing
            A.append(j + 1)
    i = 0
    while i < len(A):
        m += 3**(i)*2**(A[len(A) - i - 1] - 1)
        i += 1
    return m

def get_x_generator(B):
    r = factor(collatz(x, B[:-1]))
    #print(r)

    m,n,d = tokenize_exp(r)
    if B[-1] == '0':
        g, p, q = eea(2*d, m)
        a = int(m / g)
        b = int(n * p / g)
        b = b % a
        l = a * k + b
        v = expand((2 * d * l - n) / m)
    else:
        g, p, q = eea(2*d, m)
        a = int(m / g)
        b = int((n - d) * p / g)
        b = b % a
        l = a * k + b
        v = expand((2 * d * l - (n-d)) / m)
    return v

def get_Tx_generator(b):
    xg = get_x_generator(b)
    return expand(collatz(xg, b))

if __name__ == '__main__':
    x, k = symbols('x k')

    t0 = time.time()

    L = []

    m = 11
    B = all_bin_len(m)
    i = 0
    for b in B:
        #A =  2**(m - 1)*(collatz(x, b[:-1])-int('0b'+b[0],2))
        #xg = get_x_generator(b)
        #txg = get_Tx_generator(b)
        s = 2**len(b)
        r = get_const_r(b)
        v = 3**b.count('1')
        c = smallest_n(v, s, -r)
        t = (s * c - r) // v
        w = c

        #if (w == 1) or (v == 1 and w == 0):
            #i += 1
        if w == 1:
            print(b)
            print(r)
            print(w)
            p = (s - v)**(2*3**(b.count('1') - 1) - 1) % v
            #print(p)
            print(w * (s - v) % v)
            print(r * p % v)
            print(s-v)
            #print(A)
            #print(int('0b' + b, 2), '->', r)
            #print('x = {}'.format(str(xg).replace('*', '')))
            print('x = {}k + {}'.format(s, t))
            #print('T(x) = {}'.format(str(txg).replace('*', '')))
            print('T(x) = {}k + {}'.format(v, w))
            print('')
            L.append(t)

            #print('${}$ & ${}$ & ${}$ \\\\'.format(b, '{}k + {}'.format(s, t), '{}k + {}'.format(v, w)))
    #print(i)

    print('Took {:.4f}s'.format(time.time() - t0))
    L = sorted(list(set(L)))
    print(L)
