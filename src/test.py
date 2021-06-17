import time

from math import floor
from multiprocessing.pool import Pool
from functools import partial
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

from utils import *
from scaling_factors import get_wab
from test_r import *

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
            x = (3*x + 1) // 2**vp(3*x+1, 2)
            code += '1'
        else:
            x = x // 2**vp(x, 2)
            code += '0'
        l.append(x)

    res = {
    'iter': n,
    'code': code,
    'values': l
    }

    return res

def collatz(x):
    code = ''
    n = 0

    l = [x]
    while x != 1:
        n += 1
        if x % 2 == 1:
            x = (3*x + 1)
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

#gcd by euclidean algorithm
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

#modular inverse
def mod_inv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        return None
    else:
        return x % m

def first_n_primes(n):
    L = [2]
    m = 3
    while len(L) < n:
        if is_prime(m):
            L.append(m)
        m += 1
    return L

def pairwise_gcd(M):
    L = []
    i = 0
    while i < len(M):
        j = i + 1
        while j < len(M):
            L.append(egcd(M[i], M[j])[0])
            j += 1
        i += 1
    return L

def is_power(p, n=2):
    while p % n == 0:
        p //= n
    if p == 1:
        return True
    return False


def get_as_congruence_class(n, w):
    r = n % w
    q = (n - r) // w
    return q, r

def at_least_n_equal(L, n):
    V = list(set(L))
    C = [L.count(v) for v in V]
    for c in C:
        if c >= n:
            return True
    return False

def generate_circle_graph(w):
    if w % 3 == 0:
        return -1
    d = w % 3
    r1 = (d*w - 1) // 3
    if d == 1:
        r2 = (2*w - 1) // 3
    else:
        r2 = (w - 1) // 3
    O = collatz(w)['values']
    Ow = [o % w for o in O]
    P = {}
    angle_delta = 2*np.pi / w
    i = 0
    while i < w:
        angle = (np.pi / 2)  - i * angle_delta
        P[i] = (1 * np.cos(angle), 1 * np.sin(angle))
        i += 1

    #xy axis
    plt.plot([0, 0], [-1.2, 1.2], 'k')
    plt.plot([-1.2, 1.2], [0, 0], 'k')

    for p in P:
        if p == r1:
            plt.plot([P[p][0]], [P[p][1]], 'r.', markersize=6)
        elif p == r2:
            plt.plot([P[p][0]], [P[p][1]], 'g.', markersize=6)
        else:
            plt.plot([P[p][0]], [P[p][1]], 'b.', markersize=5)

    q_max = max([get_as_congruence_class(o, w)[0] for o in O])
    print(q_max)

    i = 1
    while i < len(Ow):
        o = O[i]
        q,r = get_as_congruence_class(o, w)
        #plt.plot([P[Ow[i-1]][0], P[Ow[i]][0]], [P[Ow[i-1]][1], P[Ow[i]][1]], 'r')
        x1 = P[Ow[i - 1]][0]
        y1 = P[Ow[i - 1]][1]
        x2 = P[Ow[i]][0]
        y2 = P[Ow[i]][1]
        plt.arrow(x1, y1, x2 - x1, y2 - y1, head_width=0.02, head_length=0.05, length_includes_head=True)
        i += 1
    #plt.show()

def generate_cylinder_graph(w):
    if w % 3 == 0:
        return -1
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    O = collatz(w)['values']
    V = [get_as_congruence_class(o, w) for o in O]
    lower_q = min([v[0] for v in V])
    upper_q = max([v[0] for v in V])
    P = {}
    angle_delta = 2*np.pi / w
    i = 0
    while i < w:
        angle = (np.pi / 2)  - i * angle_delta
        P[i] = (1 * np.cos(angle), 1 * np.sin(angle))
        i += 1

    z = lower_q
    while z <= upper_q:
        ax.scatter3D([P[i][0] for i in P], [P[i][1] for i in P], [z for _ in range(len(P))])
        z += 1

    i = 1
    while i < len(O):
        #plt.plot([P[Ow[i-1]][0], P[Ow[i]][0]], [P[Ow[i-1]][1], P[Ow[i]][1]], 'r')
        x1 = P[V[i - 1][1]][0]
        y1 = P[V[i - 1][1]][1]
        z1 = V[i - 1][0]
        x2 = P[V[i][1]][0]
        y2 = P[V[i][1]][1]
        z2 = V[i][0]
        ax.plot3D([x1, x2], [y1, y2], [z1, z2])
        i += 1
    #plt.show()

def tuple_diff(a, b):
    return tuple([i - j for i, j in zip(a, b)])

#ensure sequential makes sure that this is also valid for a, a + 1, ..., a + (ensure_sequential)
#if ensure_sequential == 0, we don't check.
def get_wa(a, w_bound=False, ensure_sequential=0, include_even=False, verbose=False):
    wab = []
    b = 2*a
    if verbose:
        print((a, b))

    r = 2**b - 3**a
    skip = 2
    if include_even:
        skip = 1
    W = range(1, 3**a, skip)
    if w_bound:
        W = range(1, min(3**a, w_bound), skip)
    for w in W:
        if in_Rn(w*r, a):
            if verbose:
                print('{} | {}: {}'.format(w, w * r, inv_r(w*r, a)))
            in_all = True
            for k in range(ensure_sequential):
                bp = 2*(a + k + 1)
                rp = 2**bp - 3**(a + k + 1)
                if not in_Rn(w*rp, a + k + 1):
                    in_all = False
                    break
            if in_all:
                wab.append(w)

    return wab

def seq_diff(a):
    b = [a[0]] + [a[i] - a[i - 1] for i in range(1, len(a))]
    return tuple(b)

def most_common(lst):
    return max(set(lst), key=lst.count)

def get_abnormal_w(upper_a, upper_w, verbose=False):
    Ws = []
    for a in range(1, upper_a):
        W = []
        b = 2 * a
        for w in range(1, min(upper_w, 3**a), 2):
            r = w * (2**b - 3**a)
            if in_Rn(r, a):
                W.append(w)
        if verbose:
            print('{}: {}'.format((a, b), W))
        Ws.append(W)

    M = []

    for i, W in enumerate(Ws[:-1]):
        for w in W:
            j = i + 1
            found = False
            while j < len(Ws):
                if w in Ws[j]:
                    found = True
                    break
                j += 1
            if not found:
                M.append(w)


    M = sorted(M)
    if verbose:
        print('Found once: {}'.format(M))

    S = []

    for i, W in enumerate(Ws[:-1]):
        for w in W:
            j = i + 1
            if w not in Ws[j]:
                j += 1
                found_again = False
                while j < len(Ws):
                    if w in Ws[j]:
                        found_again = True
                        break
                    j += 1
                if found_again:
                    S.append(w)

    S = sorted(list(set(S)))
    if verbose:
        print('Skipped: {}'.format(S))

    D = {}
    for w in S:
        A = []
        for a in range(1, upper_a):
            b = 2*a
            r1 = 1*(2**b - 3**a)
            I1 = inv_r(r1, a)

            rw = w*(2**b - 3**a)
            if not in_Rn(rw, a):
                continue
            A.append(a)
            Iw = inv_r(rw, a)
        A = seq_diff(A)[1:]
        A = tuple([a for a in A if a != 1])
        if A in D:
            D[A].append(w)
        else:
            D[A] = [w]

    if verbose:
        for d in D:
            print('{}: {}'.format(d, D[d]))

    W = []
    for w in Ws:
        W += w
    W = sorted(list(set(W).difference(set(M)).difference(set(S))))

    return W, M, S, D

def index_of_longest_seq_of(L, s):
    longest_seq = -1
    longest_index = -1

    for i, l in enumerate(L):
        if l == s:
            c = 1
            for k in range(1, len(L) - i):
                if L[i + k] == s:
                    c += 1
                else:
                    break
            if c > longest_seq:
                longest_seq = c
                longest_index = i

    return longest_index, longest_seq


if __name__ == '__main__':

    #Print table for LaTeX
    '''w = 85
    for a in range(2, 11):
        b = 2*a
        r1 = 1*(2**b - 3**a)
        I1 = inv_r(r1, a)

        rw = w*(2**b - 3**a)
        if not in_Rn(rw, a):
            continue
        Iw = inv_r(rw, a)
        #print('{}{}: {}'.format(a, ' ' * (2 - len(str(a))), tuple_diff(Iw, I1)))
        print('{} & {} & {} & {} \\\\'.format(a, tuple([i - 1 for i in I1]), tuple([i - 1 for i in Iw]), tuple_diff(Iw, I1)))

    exit()'''

    '''W = get_wa(10, 500)
    for w in W:
        print('w = {}'.format(w))
        for l in range(10):
            for t in range(10):
                for s in range(10):
                    n = s + 2*l
                    if not in_Rn(w, t+1):
                        continue
                    if not in_Rn(2**n - 3**l*w, l):
                        continue
                    print((l, s, t))
        print('')
    exit()'''

    #Check if fits pattern
    '''w = 635
    for l in range(20):
        for t in range(20):
            for s in range(20):
                n = s + 2*l
                if not in_Rn(w, t+1):
                    continue
                if not in_Rn(2**n - 3**l*w, l):
                    continue
                print((l, s, t))
    exit()'''

    #Check if each w with w(2^2a - 3^a) a-special has (t+1)-special representation of w
    '''a = 40
    b = 2*a
    W = get_wa(a, 50000, ensure_sequential=0)[1:]
    #print(W)
    I1 = [i - 1 for i in inv_r(1*(2**b - 3**a), a)]
    I1_minus = [i - 1 for i in inv_r(1*(2**(b - 2) - 3**(a - 1)), a - 1)]
    for w in W:
        Iw = [i - 1 for i in inv_r(w*(2**b - 3**a), a)]
        Idiff = tuple_diff(Iw, I1)
        s_guess_1 = most_common(Idiff)

        Iw_minus = [i - 1 for i in inv_r(w*(2**(b - 2) - 3**(a - 1)), a - 1)]
        Idiff_minus = tuple_diff(Iw_minus, I1_minus)
        for i, (p, q) in enumerate(zip(Idiff, Idiff_minus)):
            if p != q:
                break
        s_guess_2 = Idiff[i]
        if s_guess_1 != s_guess_2:
            print('s guess doesnt match for {}: {} | {}, {}'.format(w, Idiff, s_guess_1, s_guess_2))
            continue
        s = s_guess_2
        index, seq_len = index_of_longest_seq_of(Idiff, s)
        t = a - (index + seq_len)
        l  = index
        #print('{}: {} | ({}, {}, {})'.format(w, Idiff, l, s, t))
        Iacc = []
        for i,v in enumerate(Iw[a - t:][::-1]):
            Iacc.append(v - 2*a)
        Iacc.append(s-2*t)
        Iacc = Iacc[::-1]
        Iacc = tuple(Iacc)
        if not in_Rn(w, t + 1):
            print('{}: -- | ({}, {}, {}) **'.format(w, l, s, t))
            print('  We get {}, but w does not have an {}-special rep.'.format(Iacc, t+1))
            print()
        else:
            Ir1 = [i - 1 for i in inv_r(w, t + 1)]
            Ir1 = tuple(Ir1)
            if Ir1[-1] + 2 != Idiff[-1]:
                print('{}: -- | ({}, {}, {}) **'.format(w, l, s, t))
                print('  Should have {}, but we get {}'.format(Ir1, Iacc))
                print('  These give us {} and {}'.format(r_func([i + 1 for i in Ir1]), r_func([i + 1 for i in Iacc])))
                print()


    exit()'''



    #See outputs of pattern
    '''w = 4835
    a = 12
    l, s, t = (3,11,5)
    n = s + 2*l
    I1 = tuple([i - 1 for i in inv_r(w, t+1)])
    I2 = tuple([i - 1 for i in inv_r(2**n - 3**l*w, l)])
    Vsub = [i - 2*j for j, i in enumerate(I2)]
    for _ in range(a - l - t):
        Vsub.append(s)
    for j, i in enumerate(I1[1:]):
        Vsub.append(2*(t - j) + i)
    Vsub = tuple(Vsub)
    Vw = tuple([2*j + i for j,i in enumerate(Vsub)])
    print('w = {}'.format(w))
    print('(l, s, t) = {}'.format((l, s, t)))
    print('n = {}'.format(n))
    print('I1 = {}'.format(I1))
    print('I2 = {}'.format(I2))
    print('V{} - V1 = {}'.format(w, Vsub))
    print('V{} = {}'.format(w, Vw))
    V = tuple([i + 1 for i in Vw])
    if r_func(V) == w*(4**a - 3**a):
        print("MATCH")
    exit()'''



    #Find s, t to satisfy condition 1
    #w = 725
    #w = 821
    '''w = 397
    print('w = {}\n'.format(w))
    for t in range(int(np.log(w)/np.log(3)) + 1):
        for s in range(2*t, 2*t + 20):
            r = w - 3**t * 2**(s - 2*t)
            if in_Rn(r, t):
                print('s = {} | t = {}'.format(s, t))
                print(tuple([i - 1 for i in inv_r(r, t)]))
                print()
    exit()'''

    #Find l, s to satisfy condition 2
    Sol = []
    w = 821
    for l in range(20):
        for s in range(20):
            r = 2**(s + 2*l) - 3**l * w
            if in_Rn(r, l):
                Sol.append((l, s))
                print('l = {} | s = {}'.format(l, s))
                print(tuple([i - 1 for i in inv_r(r, l)]))
                print()
    S = sorted(list(set([sol[1] for sol in Sol])))
    print('{}: {}'.format(w, S))
    exit()

    #Just a single w
    '''#for w in [1, 2, 5, 11, 19, 53, 85, 149, 151, 325, 331, 341, 397]:
    w = 19
    print('w = {}'.format(w))
    for a in range(1, 50):
        b = 2*a
        r1 = 1*(2**b - 3**a)
        I1 = inv_r(r1, a)
        I1 = tuple([i - 1 for i in I1])

        rw = w*(2**b - 3**a)
        if not in_Rn(rw, a):
            continue
        Iw = inv_r(rw, a)
        Iw = tuple([i - 1 for i in Iw])
        #print('{}{}: {}'.format(a, ' ' * (2 - len(str(a))), tuple_diff(Iw, I1)))
        #print('{} & {} & {} & {} \\\\'.format(a, tuple([i - 1 for i in I1]), tuple([i - 1 for i in Iw]), tuple_diff(Iw, I1)))
        print('{}{}: {}'.format(a, ' ' * (2 - len(str(a))), tuple_diff(Iw, I1)))
        #print('{}{}: {}'.format(a, ' ' * (2 - len(str(a))), Iw))'''

    #For a bunch of numbers
    '''#2
    W = [325, 397, 3095, 3941, 5165, 5357, 6853, 6965, 7375, 8405, 8747, 8945, 11149, 11285, 14677, 15245, 16405, 16621, 16813, 17197, 17323, 18829, 20365, 20725, 21061, 21191, 21427, 21571, 21677, 21877, 22807, 26197, 26965, 27491, 29095, 29581, 30085, 34061, 36551, 36659, 36743, 37445, 38509, 40811, 41131, 47045, 48899, 49799, 51797, 52117, 52885, 57389, 57677, 61525, 63253, 66043, 66577, 66947, 68159, 69061, 69269, 69409, 69493, 69517, 69859, 69889, 72103, 72397, 82261, 84037, 84053, 84341, 84367, 84877, 85061, 85141, 85205, 86257, 86627, 87797, 88499, 89335, 91201, 91237, 92101, 92981, 93917, 105557, 110515, 111637, 112405, 113557, 114269, 115597, 118135, 118919, 119215, 120377, 122645, 122773, 123157, 123445, 123509, 127597, 131087, 137549, 138125, 138325, 141155, 145781, 146317, 146393, 146701, 147269, 148661, 149141, 149347, 149465, 149909, 150641, 150925, 152371, 152821, 156359, 156815, 157331, 157519, 157571, 158627, 164117, 164405, 164551, 164693, 165613, 166963, 179597, 188419, 195491, 196949, 197635, 199349, 199829, 201005, 201233, 202133, 204509, 206741, 214019, 220301, 223157, 229655, 229799, 230087, 234983, 236657, 238477, 240853, 241037, 241625, 250157, 255545, 255707, 255923, 263125, 264037, 264293, 264325, 264373, 264481, 265201, 265345, 265493, 265805, 266125, 267077, 267637, 268181, 268357, 268501, 269845, 270061, 270805, 270853, 270877, 272497, 272621, 272677, 273031, 276205, 278431, 279475, 279799, 279823, 280469, 282055, 283799, 285457, 286765, 288493, 289159, 289345, 290005, 291805, 292165, 292693, 292885, 294541]
    #5
    #W = [545, 547, 2161, 2801, 5219, 8717, 20921, 40183, 44597, 46739, 46913, 47381, 47611, 60401, 65965, 67513, 79825, 83341, 83495, 84277, 90553]
    #6, 2
    #W = [1261, 11975]
    #2, 5
    #W = [2243]
    #7
    #W = [2987, 4153, 5291, 6721, 12065, 16651, 36131, 47971, 48259, 67009]
    #4
    #W = [4621, 45397, 46901, 70807, 71749, 73453, 79399, 79669, 82397]
    #3, 2
    #W = [11851, 21067, 29761, 35915, 36821, 37451, 84821]
    #4, 2
    #W = [37037]
    #3
    #W = [59765, 67717, 70597, 72001, 72541, 80245]
    #5, 2
    #W = [69773]
    #4, 5
    #W = [79681]
    #2, 2, 3
    #W = [81613]
    #6
    #W = [96493]
    for w in W:
        #print('w = {}'.format(w))
        A = []
        T = []
        c = 0
        for a in range(2, 40):
            if c > 8:
                break
            b = 2*a
            r1 = 1*(2**b - 3**a)
            I1 = inv_r(r1, a)

            rw = w*(2**b - 3**a)
            if not in_Rn(rw, a):
                continue
            A.append(a)
            c += 1
            Iw = inv_r(rw, a)
            #print('{}{}: {}'.format(a, ' ' * (2 - len(str(a))), tuple_diff(Iw, I1)))
            T.append(tuple_diff(Iw, I1))
            #print('{}{}: {}'.format(a, ' ' * (2 - len(str(a))), tuple_diff(seq_diff(Iw), seq_diff(I1))))
        A = seq_diff(A)[1:]
        A = [a for a in A if a != 1]

        type = -1
        if A == [2]:
            i = 0
            while T[0][i] == T[1][i]:
                i += 1
            if len(T) > 2:
                if T[2][i] != T[1][i]:
                    i -= 1
            leading = tuple(T[1][:i])
            skip = tuple(T[1][i:i + 2])
            if len(T) > 2:
                repeat = T[2][i + 1]
            else:
                repeat = 'n/a'
            i += 2
            trailing = tuple(T[1][i:])
            orig_trailing = tuple(T[0][i-2:])
            if skip[0] + 3 == skip[1]:
                type = 1
            if skip[0] + 6 == skip[1]:
                type = 2
            if skip[0] == 9 and skip[1] == 10:
                type = 3

        if type == -1:
            print('w = {}'.format(w))
            for t in T:
                print('{}{}: {}'.format(len(t), ' ' * (2 - len(str(len(t)))), t))
            if orig_trailing == trailing:
                print('({}, {}, {}, {})'.format(leading, skip, repeat, trailing))
            else:
                print('({}, {}, {}, {})**'.format(leading, skip, repeat, trailing))
                print(orig_trailing)
            print('')

    exit()'''

    '''W, M, S, D = get_abnormal_w(40, 3*10**5, verbose=False)

    for d in D:
        print('{}: {}'.format(d, D[d]))'''
