from multiprocessing.pool import Pool
from functools import partial
import os

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def is_prime(x):
    if x < 2:
        return False
    if x in [2,3,5]:
        return True
    if x % 2 == 0:
        return False
    i = 3
    u = int(np.sqrt(x))
    while i < u + 1:
        if x % i == 0:
            return False
        i += 2
    return True

def get_primes_up_to(n):
    l = [2]
    i = 3
    while i <= n:
        if is_prime(i):
            l.append(i)
        i += 2
    return l

def get_factor_vector(x, pad_to=False):
    n = int(np.sqrt(x)) + 1

    p = []
    i = 2
    if pad_to:
        while i <= pad_to:
            if is_prime(i):
                p.append(i)
            i += 1
    else:
        while i <= x:
            if is_prime(i):
                p.append(i)
            i += 1

    f = {}

    i = 0
    while i < len(p):
        c = 0
        while x % p[i] == 0 and x != 0:
            c += 1
            x //= p[i]
        f[p[i]] = c
        i += 1

    return f

def factor_prod(f):
    p = 1
    for i, v in zip(f.keys(), f.values()):
        p *= i**v
    return p

def dict_to_vec(f):
    return list(f.values())

#F(x) = ax + b if odd, x / 2 if even
def collatz_mod(x, mod=17, return_list=False):
    code = ''
    failed = False
    n = 0
    l = [x]
    while x != 1:
        if x % 2 == 1:
            x = (3*x + 1) % mod
            code += '1'
        else:
            x = x // 2
            #x = (x * (mod + 1) // 2) % mod
            code += '0'
        if x in l:
            failed = True
            break

        l.append(x)
        n += 1

    res = {
    'iter': n,
    'code': code,
    'values': l,
    'failed': failed
    }

    return res

#0 is x/2, 1 is 3x+1
def collatz(x):
    code = ''
    n = 0

    l = [x]
    while x != 1:
        n += 1
        if x % 2 == 1:
            x = 3*x + 1
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

def collatz_T(x):
    code = ''
    n = 0

    l = [x]
    while x != 1:
        n += 1
        if x % 2 == 1:
            x = (3*x + 1) // 2
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

def collatz_phi(x):
    code = ''
    iter = 0

    l = [x]
    while x != 1:
        iter += 1
        if x % 2 == 1:
            n = phi(x)
            x = (3*x + 1) // 2**n
            code += '1'
        else:
            n = v2(x)
            x = x // 2**n
            code += '0'
        l.append(x)

    res = {
    'iter': iter,
    'code': code,
    'values': l
    }

    return res

def collatz_psi(x):
    code = ''
    n = 0

    l = [x]
    while x != 1:
        n += 1
        if x % 2 == 1:
            n = psi(x)
            x = (3*x + 1) // 2**n
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

def reverse_colaltz(code):
    x = 1
    uses_exact = True
    for c in code:
        if c == '1' and ((x - 1) % 3 != 0 or (x - 1) == 0):
            return -1
        if c == '0':
            x *= 2
        else:
            x = (x - 1) // 3
    if x == 2 and len(code) > 1 or x == 4 and len(code) > 2:
        return -1
    return x

def find_all_num_of_len(length):
    #codes = [bin(x)[2:].rjust(length, '0') for x in range(2**(length + 1))]

    nums = set()
    i = 0
    while i < 2**(length + 1):
        code = bin(i)[2:].rjust(length, '0')
        r = reverse_colaltz(code)
        if r != -1:
            l = collatz(r)['iter']
            if l == length:
                nums.add((r, code))
        i += 1

    return nums

def _find_all_num_of_len(length, skip, offset):
    #codes = [bin(x)[2:].rjust(length, '0') for x in range(2**(length + 1))]

    nums = set()
    i = offset
    while i < 2**(length + 1):
        code = bin(i)[2:].rjust(length, '0')
        r = reverse_colaltz(code)
        if r != -1:
            l = collatz(r)['iter']
            if l == length:
                nums.add((r, code))
        i += skip

    return nums

def get_coord(x):
    code = collatz(x)['code']
    x_coord = np.sum(np.array(list(code)) == '1')
    y_coord = np.sum(np.array(list(code)) == '0')

    return (x_coord, y_coord)

def find_all_num_of_len_multi(length, n_threads=8):
    func = partial(_find_all_num_of_len, length, n_threads)

    params = [x for x in range(n_threads)]

    pool = Pool()
    res = pool.map(func, params)
    pool.close()
    pool.join()

    nums = res[0]
    for r in res[1:]:
        nums.update(r)


    return nums

def build_graph(radius):
    G = {}
    G[1] = []

    r = 1
    while r <= radius:
        for k in list(G.keys()):
            nk = 2 * k
            if nk in G:
                if not k in G[nk]:
                    G[nk].append(k)
                    G[k].append(nk)
            else:
                G[nk] = [k]
                G[k].append(nk)

            if (k - 1) % 3 == 0 and (k - 1) != 0:
                nk = (k - 1) // 3
                if nk in G:
                    if not k in G[nk]:
                        G[nk].append(k)
                        G[k].append(nk)
                else:
                    G[nk] = [k]
                    G[k].append(nk)
        r += 1

    nx.draw(nx.Graph(G), node_size=10)
    plt.show()

def code_complement(code):
    comp = ''.join([str(1 - int(x)) for x in list(code)])
    return comp

def plot_path(n, show=True):
    V = collatz(n)['values']
    coords = []
    for v in V[::-1]:
        coords.append(get_coord(v))
    X = [x[0] for x in coords]
    Y = [x[1] for x in coords]

    plt.plot(X, Y, label=n)
    if show:
        plt.show()

def plot_paths(N):
    for n in N:
        plot_path(n, show=False)
    plt.legend()
    plt.show()

def v2(x):
    if x == 0:
        return 0
    i = 1
    while True:
        if x % 2**i != 0:
            return i - 1
        i += 1

def plot_iteration(n, show=True):
    V = collatz(n)['values']
    coords = []
    for i, v in enumerate(V):
        coords.append((v, i))
    X = [x[0] for x in coords]
    Y = [x[1] for x in coords]
    plt.plot(X, Y)
    if show:
        plt.show()

def plot_iterations(N):
    for n in N:
        plot_iteration(n, show=False)
    plt.grid(True, axis='y', alpha=1)
    plt.show()

def plot_ordered_iterations(N, subfolder=None, save=False):
    if subfolder is None:
        #width and skip
        subfolder = 'w' + str(max(N) - min(N) + 1) + '-s' + str(N[1] - N[0])
    if not os.path.exists('./images/ordered-iteration/{}'.format(subfolder)):
        os.makedirs('./images/ordered-iteration/{}'.format(subfolder))

    V = []
    for n in N:
        V.append(collatz(n)['values'])

    m = min([len(x) for x in V])

    X = []
    Y = []
    i = 0
    while i < len(V):
        X.append([])
        Y.append([])
        i += 1
    i = 0
    while i < m:
        p = [x[i] for x in V if len(x) > i]
        S = np.argsort(np.array(p))

        for j, s in enumerate(S):
            X[s].append(j)
            Y[s].append(i)
        i += 1

    plt.figure(figsize=(7,11), dpi=160)
    for n, x, y in zip(N, X, Y):
        if True:
            plt.plot(x, y)
    if not save:
        plt.show()
    else:
        plt.savefig('./images/ordered-iteration/{}/{}-{}.png'.format(subfolder, min(N), max(N)))
    plt.close()

#phi(x) = max{y : 2^y | 3x + 1}
def phi(x):
    if x % 2 == 0:
        return 0
    x = 3 * x + 1
    i = 0
    while x % 2 == 0:
        x //= 2
        i += 1
    return i

def psi(x):
    x = 3*x + 1
    n = 0
    while True:
        if (x - 2**n) % 2**(n + 1) == 0:
            return n
        n += 1

def num_to_even_odd_coord(y):
    n = 0
    while True:
        if (y - 2**n) % 2**(n + 1) == 0:
            k = (y // 2**n - 1) // 2
            return n, k
        n += 1

def even_odd_coord_to_num(n, k):
    return 2**n * (2 * k + 1)

def even_odd_coord_graph(N):
    C = []
    Z = []
    for n in N:
        C.append(num_to_even_odd_coord(n))
        Z.append(collatz(n)['iter'])
    X = [c[0] for c in C]
    Y = [c[1] for c in C]
    m = max(Z)
    Z = [z/m for z in Z]
    plt.scatter(X, Y, c=[(z,1 - z,0) for z in Z])
    plt.show()

def even_odd_coord_graph_grid(x_max, y_max):
    X = []
    Y = []
    Z = []
    i = 0
    while i < x_max:
        j = 0
        while j < y_max:
            X.append(i)
            Y.append(j)
            Z.append(collatz(even_odd_coord_to_num(i, j))['iter'])
            j += 1
        i += 1
    m = max(Z)
    Z = [z/m for z in Z]
    plt.scatter(X, Y, c=[(z,1 - z,0) for z in Z])
    plt.show()

def mod_inv(x, n):
    for i in range(0, n):
        if i*x % n == 1:
            return i
    return -1

def check_finite_collatz(p, l=1):
    N = np.arange(2, p**l, 1)
    B = []
    for n in N:
        if collatz_mod(n, mod=p**l)['failed']:
            B.append(n)

    if len(B) > 0:
        return False, B
    return True, B

def has_odd_fixed_point(n):
    for i in range(2, n):
        if (3*i + 1) % n == i:
            if i % 2 == 1:
                return True, i
    return False, -1

def has_odd_congruent_to_n(n):
    for i in range(2, n):
        if (3*i + 1) % n == 0:
            if i % 2 == 1:
                return True, i
    return False, -1

def check_finite_conjectures(n, l):
    t, v = has_odd_fixed_point(n**l)
    if gcd(n**l, 2) == 1 and (n**l - 1) % 4 != 0 and n != 3:
    #if t:
        #print("{} failed odd fixed points".format(n**l))
        return False, v
    t, v = has_odd_congruent_to_n(n**l)
    if gcd(n**l, 3) == 1 and (n**l - 1) % 6 != 0 and n != 4:
    #if gcd(n**l, 3) == 1:
    #if t:
        #print("{} failed odd congruence to n".format(n**l))
        return False, v
    return True, -1

def gcd(x, y):
   while(y):
       x, y = y, x % y
   return x

def eea(a, b):
    """return (g, x, y) such that a*x + b*y = g = gcd(a, b)"""
    x0, x1, y0, y1 = 0, 1, 1, 0
    while a != 0:
        (q, a), b = divmod(b, a), a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return b, x0, y0

def generate_finite_graph(m, save=False):
    N = np.arange(1, m)
    D = {}
    for n in N:
        if n % 2 == 1:
            D[n] = [(3 * n + 1) % m]
        else:
            D[n] = [n // 2]

    nx.draw(nx.DiGraph(D).to_directed(), node_size=10, arrows=True)
    if not save:
        plt.show()
    else:
        plt.savefig('./images/finite-graph/{}.png'.format(m))
    plt.close()

def polar_to_xy(theta, r):
    return r * np.cos(theta), r * np.sin(theta)

def collatz_mod_step(x, n):
    if x % 2 == 0:
        return x // 2
    return (3*x + 1) % n

def generate_circle_diagram(n, save=False):
    angles = np.linspace(0, 2*np.pi, n)
    coords = [polar_to_xy(a, 1) for a in angles]

    dest = [collatz_mod_step(i, n) for i, a in enumerate(angles)]

    X = [c[0] for c in coords]
    Y = [c[1] for c in coords]

    for i, (x, y, d) in enumerate(zip(X, Y, dest)):
        plt.plot([x, X[d]], [y, Y[d]])
        #plt.arrow(x, y, X[d] - x, Y[d] - y, width=0.01)

    plt.scatter(X, Y)
    if save:
        plt.savefig('./images/zn-circle/{}.png'.format(n))
    else:
        plt.show()
    plt.close()

#Checks highest n such thtt 2^n | 3^(n-1)(x + 1)
def divergence_degree(x):
    n = 1
    while True:
        if not (3**(n - 1) * (x + 1)) % 2**n == 0:
            return n - 1
        n += 1
