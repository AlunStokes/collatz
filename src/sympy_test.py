from sympy import *
from sympy.utilities.lambdify import lambdify, implemented_function

from utils import *

init_printing(use_unicode=True)

'''
def code_to_lambda(code, x):
    if code[0] == '0':
        fi = implemented_function('f', lambda x: x / 2)
        f = lambdify(x, fi(x))
    else:
        fi = implemented_function('f', lambda x: (3*x + 1) / 2)
        f = lambdify(x, fi(x))
    for c in code[1:]:
        if c == '0':
            fi = implemented_function('f', lambda x: f(x) / 2)
            f = lambdify(x, fi(x))
        else:
            fi = implemented_function('f', lambda x: (3*f(x) + 1) / 2)
            f = lambdify(x, fi(x))

    return f
'''
def collatz(x, code):
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

def get_const(b):
    m = 0
    A = []
    for j, i in enumerate(b[1:]):
        if i == '1':
            #Extra for skipping first index, second extra for 1-indexing
            A.append(j + 2)
    i = 0
    while i < len(A):
        m += 3**(i)*2**(len(b)-A[i])
        i += 1
    if b[0] == '1':
        m -= 2**(len(b) - 1)
    return m

if __name__ == '__main__':
    '''x, y = symbols('x y')

    V = [8*x, 8*x+1, 8*x + 2, 8*x + 3, 8*x+4, 8*x+5, 8*x + 6, 8*x + 7]
    C = ['000', '101', '010', '110', '001', '100', '011', '111']

    for v,c in zip(V, C):
        print(factor(collatz(v, (c))))'''

    x, k = symbols('x k')

    V = []
    C = []

    bins = all_bin_len(3)

    for B in bins:
        C.append(B[::-1])
        r = factor(collatz(x, B[::-1][:-1]))
        #print(r)
        m,n,d = tokenize_exp(r)
        if B[0] == '0':
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
        V.append(v)

    R = [expand(collatz(v, (c))) for v,c in zip (V,C)]
    #print(len([x for x in R if R.count(x) == 1]))

    s1 = max([len(str(s)) for s in V])
    s2 = max([len(str(s)) for s in R])


    n = 0
    for i, (v,c,r) in enumerate(zip(V, C, R)):
        #r = expand(collatz(v, (c)))
        matches = int('0b' + c[::-1], 2) == tokenize_relation(v)[1]
        #print(matches)
        if matches:
            print('{}: {}{} ->  {}{}*'.format(c, v, ' ' * (s1 - len(str(v)) + 1), r, ' ' * (s2 - len(str(r)) + 1)))
        else:
            n += 1
            print('{}: {}{} ->  {}{}'.format(c, v, ' ' * (s1 - len(str(v)) + 1), r, ' ' * (s2 - len(str(r)) + 1)))
