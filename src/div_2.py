import time

from math import floor
from multiprocessing.pool import Pool
from functools import partial
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from itertools import combinations, permutations, product

from utils import *
from scaling_factors import get_wab
from garb import nu
from scaling_factors import *
from test_r import get_all_r

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

def is_power_of_p(n, p):
    a = 0
    while n % p**a == 0:
        a += 1
    a -= 1
    if n == p**a:
        return True
    return False

def tuple_diff(a, b):
    return tuple([i - j for i, j in zip(a, b)])


def seq_diff(a):
    b = [a[0]] + [a[i] - a[i - 1] for i in range(1, len(a))]
    return tuple(b)

#Uses 0 indexing
def inv_r(r, n):
    if n == 0:
        if r == 0:
            return tuple()
        else:
            return (-1,)
    if n == 1:
        v = nu(r, 2)
        if 2**v == r:
            return (v,)
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
        L.append(nu(r - s) + 1)
        i += 1
    return tuple([i - 1 for i in L])

#Uses 0 indexing
def r_func(A, right_offset=0):
    r = 0
    for i, v in enumerate(A[::-1]):
        r += 2**(v) * 3**(i + right_offset)
        #r += 3**(len(A) - 1 - i - offset) * 2**(A[i])
    return r

#Considers which parts of original rep are given by comb
def r_func_comb(A, comb, total_len=-1):
    if total_len == -1:
        total_len = len(A)
    r = 0
    for c, v in zip(comb[::-1], A[::-1]):
        r += 2**(v) * 3**(total_len - c - 1)
        #r += 3**(len(A) - 1 - i - offset) * 2**(A[i])
    return r

def r_func_comb(A, comb):
    total_len = len(A)
    A = [A[i] for i in comb]
    r = 0
    for c, v in zip(comb[::-1], A[::-1]):
        r += 2**(v) * 3**(total_len - c - 1)
        #r += 3**(len(A) - 1 - i - offset) * 2**(A[i])
    return r

def get_all_r(a, b, include_even = False, starting_with=tuple()):
    R = []
    Is = []
    if len(starting_with) > 0:
        comb = [(*starting_with, *c) for c in list(combinations(list(range(starting_with[-1] + 1, b + 1)), a - len(starting_with)))]
    elif include_even:
        comb = list(combinations(list(range(0, b + 1)), a))
    else:
        comb = [(0, *c) for c in list(combinations(list(range(1, b + 1)), a - 1))]
    if len(starting_with) > 0:
        starting_length = len(starting_with)
        comb = [c for c in comb if c[:starting_length] == starting_with]
    for I in comb:
        R.append(r_func(I))
        Is.append(I)
    return R, Is

def get_all_comb_sum(v1, total, max_len):
    C = []
    for l in range(1, max_len):
        nums = list(range(1, total - v1 - l + 2))
        result = [(v1, *seq) for seq in product(nums, repeat=l) if sum(seq) == total - v1]
        C += result
    C = sorted(list(set(C)))
    D = {}
    for c in C:
        if not len(c) in D:
            D[len(c)] = [c]
        else:
            D[len(c)].append(c)
    return C, D

def get_matching_patterns(I, Maps):
    Patterns = []
    for pattern in Maps:
        match = True
        sub_index = pattern.index('_')
        v = I[sub_index]
        if v not in Maps[pattern][1]:
            match = False
        i = sub_index + 1
        while i < len(pattern):
            if pattern[i] == '*':
                i += 1
                continue
            v = v + pattern[i]
            if I[i] != v:
                match = False
                break
            i += 1
        if match:
            Patterns.append(pattern)
    return Patterns

def apply_pattern(I, pattern, Maps):
    #Assume pattern applies, no checking
    Is = []
    sub_index = pattern.index('_')
    v = I[sub_index]
    if not v in Maps[pattern][1]:
        return []
    Res_patterns = Maps[pattern][0]
    for res_pattern in Res_patterns:
        J = [v]
        for j in res_pattern[1:]:
            J.append(J[-1] + j)
        J = tuple(J)
        Is.append(J)
    return Is

def pattern_to_comb(pattern):
    c = []
    for i, j in enumerate(pattern):
        if j != '*':
            c.append(i)
    return c

def get_Idw(check_d, w):
    a_offset = 0
    while True:
        a = max(5, check_d - 2) + a_offset
        R, Is = get_all_r(a, 2*a + 1, include_even=False)
        L = []
        for r, I in zip(R, Is):
            r_ = r + 3**a * w
            d = nu(r, 2)
            d_ = nu(r_, 2)
            #d = nu(r, 4)
            #d_ = nu(r_, 4)
            if d_ == check_d:
                L.append(I)
        if len(L) in [0, 1]:
            a_offset += 1
            continue
        a_offset = 0
        matches_up_to = -1
        for i in range(a):
            match = True
            v = L[0][i]
            for I in L:
                if I[i] != v:
                    match = False
                    break
            if match:
                matches_up_to += 1
            else:
                break
        break
    return L[0][:matches_up_to + 1]


if __name__ == '__main__':

    #DELTE ME LATER
    R = {(1,1): 0, (1,2): -1, (2,1): 2, (2,2): 1}
    C = []
    S = []
    a = 8
    comb = list(product([1,2], repeat=2*a))
    for c in comb:
        if c[-1] == 1:
            continue
        s = []
        if c[0] == 1:
            s.append(2)
        else:
            s.append(1)
        for i in range(len(c) - 1):
            s.append(s[-1] + R[tuple(c[i:i + 2])])
            if s[-1] > a + 1:
                break
        if s[-1] != a:
            continue
        C.append(c)
        S.append(s)

    with_sub_count = 0
    for c, s in zip(C, S):
        K = []
        for k in range(a - 1):
            if s[2 * k + 1] == k + 1:
                K.append(k)
        if len(K) == 0:
            continue
        with_sub_count += 1
        print('{}: {}'.format(c, s))
        for k in K:
            print('  a = {} gives {}: {} | {}'.format(k + 1, c[:2*k + 2], s[:2*k + 2], c[2*k + 2:]))
        print()

    print(len(C))
    print('{} comb in total, {} of which can be broken into subseqs'.format(len(C), with_sub_count))

    exit()


    #Checking for repreating 1,2,1,2
    '''N = 30
    F = []
    Dx = []
    Dy = []
    #Check for oscillating 1,2,1,2
    for w in range(1, 10**7, 1):
        w0 = w
        v1 = 1 if w % 4 in [0, 1] else 2
        if w % 2 == 0:
            w //= 2
        else:
            w = (3*w + 1) // 2
        v2 = 1 if w % 4 in [0, 1] else 2
        follows = True
        i = 1
        while i < N:
            if (v1, v2) not in [(1,2), (2,1)]:
                follows = False
                break
            v1 = v2
            if w % 2 == 0:
                w //= 2
            else:
                w = (3*w + 1) // 2
            v2 = 1 if w % 4 in [0, 1] else 2
            i += 1
        if i - 1 > 16:
            Dx.append(w0)
            Dy.append(i - 1)
        if follows:
            F.append(w0)
    print(F)

    for x in Dx:
        O = collatz_T(x)
        code = O['code']
        code = [int(x) for x in list(code)]
        O = O['values']
        types = [(1 if o%4 in [0,1] else 2, o%4) for o in O]
        print('{}: {}'.format(x, types))
        print('  {}'.format(code))
        print()

    plt.plot(Dx,Dy)
    plt.show()
    exit()'''

    #Checking for repreating 1,2,1,2 - but bigger
    '''N = 30
    F = []
    Dx = []
    Dy = []
    #Check for oscillating 1,2,1,2
    for w in range(1, 10**6, 1):
        w0 = w
        v1 = (1 if w % 4 in [0, 1] else 2, w%4)
        if w % 2 == 0:
            w //= 2
        else:
            w = (3*w + 1) // 2
        v2 = (1 if w % 4 in [0, 1] else 2, w%4)
        follows = True
        i = 0
        while i < N:
            if (v1, v2) not in [((1,1),(2,2)), ((2,2),(1,1))] and i > 0:
                follows = False
                break
            v1 = v2
            if w % 2 == 0:
                w //= 2
            else:
                w = (3*w + 1) // 2
            v2 = (1 if w % 4 in [0, 1] else 2, w%4)
            i += 1
        if i - 1 > 12:
            Dx.append(w0)
            Dy.append(i - 1)
        if follows:
            F.append(w0)
    print(F)

    for x in Dx:
        O = collatz_T(x)
        code = O['code']
        code = [int(x) for x in list(code)]
        O = O['values']
        types = [(1 if o%4 in [0,1] else 2, o%4) for o in O]
        print('{}: {}'.format(x, [(t, c) for  t,c in zip(types, code)]))
        print()

    exit()'''


    '''R = {(1,1): 0, (1,2): -1, (2,1): 2, (2,2): 1}
    a = 8
    w = 565
    O = collatz_T(w)['values']
    T = [1 if o % 4 in [0, 1] else 2 for o in O]
    while len(T) < 2*a:
        T += [2, 1]
    S = []
    if T[0] == 1:
        S.append(2)
    else:
        S.append(1)
    for i in range(0, len(T) - 1):
        S.append(S[-1] + R[tuple(T[i: i + 2])])
    print('a = {}, w = {}'.format(a, w))
    print(T)
    print(S)
    print(S[2*a - 1])
    if S[2*a - 1] < a:
        print('Has solution for m < a (m = a - {})'.format(a - S[2*a - 1]))
    elif S[2*a - 1] == a and T[2*a - 1] == 1:
        print('Has solution for m = a, but either stabilizes or decreases')
    elif S[2*a - 1] == a + 1 and T[2*a - 1] == 1:
        print('Has solution for m = a, but then r/4^a is even')
    exit()'''


    '''Residues = []
    Residues.append([1])
    Residues.append([1, 19, 29])
    Residues.append([1, 11, 25, 29, 37, 67, 81, 87, 109, 115])
    Residues.append([1, 15, 33, 37, 49, 51, 59, 89, 115, 135, 145, 147, 153, 163, 185, 203, 209, 221, 259, 267, 279, 285, 299, 317, 323, 333, 343, 375, 389, 429, 449, 457, 485, 493, 501])
    Residues.append([1, 19, 43, 49, 51, 65, 67, 79, 119, 153, 171, 193, 195, 203, 217, 247, 271, 279, 291, 293, 345, 371, 381, 399, 419, 421, 429, 437, 457, 475, 499, 501, 517, 523, 551, 573, 599, 609, 635, 645, 651, 657, 659, 661, 701, 711, 727, 733, 749, 761, 801, 811, 827, 835, 863, 877, 887, 899, 903, 929, 953, 955, 961, 973, 977, 1027, 1031, 1039, 1053, 1061, 1081, 1101, 1105, 1113, 1123, 1125, 1139, 1181, 1187, 1201, 1203, 1253, 1281, 1291, 1329, 1331, 1341, 1349, 1355, 1367, 1385, 1409, 1413, 1431, 1433, 1483, 1517, 1545, 1559, 1561, 1569, 1581, 1595, 1653, 1671, 1709, 1721, 1737, 1745, 1747, 1781, 1795, 1805, 1809, 1821, 1865, 1869, 1885, 1899, 1937, 1939, 1967, 1975, 2013, 2021, 2033])
    Residues.append([1, 25, 57, 59, 65, 67, 87, 89, 105, 157, 171, 203, 257, 259, 271, 289, 329, 361, 363, 371, 387, 389, 391, 431, 459, 475, 493, 501, 531, 557, 561, 563, 573, 581, 587, 589, 609, 633, 661, 665, 689, 691, 697, 711, 735, 757, 797, 811, 847, 861, 877, 881, 885, 933, 969, 973, 977, 997, 1015, 1067, 1081, 1103, 1113, 1149, 1169, 1171, 1181, 1197, 1239, 1271, 1273, 1281, 1293, 1297, 1303, 1319, 1369, 1375, 1385, 1405, 1413, 1419, 1441, 1461, 1473, 1483, 1497, 1501, 1517, 1543, 1573, 1575, 1581, 1601, 1603, 1669, 1675, 1721, 1723, 1773, 1781, 1797, 1807, 1821, 1847, 1879, 1885, 1907, 1911, 1977, 1979, 2021, 2077, 2081, 2091, 2109, 2127, 2183, 2197, 2211, 2277, 2295, 2315, 2327, 2329, 2373, 2379, 2381, 2393, 2403, 2405, 2429, 2485, 2487, 2509, 2513, 2555, 2583, 2585, 2599, 2619, 2623, 2633, 2667, 2677, 2693, 2699, 2701, 2711, 2735, 2755, 2797, 2817, 2819, 2835, 2889, 2923, 2935, 2959, 2971, 2989, 2991, 3001, 3003, 3019, 3059, 3091, 3101, 3117, 3121, 3123, 3191, 3195, 3225, 3235, 3237, 3261, 3289, 3293, 3301, 3313, 3317, 3339, 3395, 3397, 3403, 3421, 3465, 3491, 3493, 3529, 3543, 3577, 3589, 3595, 3599, 3605, 3607, 3609, 3619, 3627, 3643, 3665, 3667, 3679, 3699, 3701, 3729, 3731, 3745, 3799, 3833, 3843, 3847, 3881, 3899, 3901, 3913, 3929, 3935, 3947, 3969, 4001, 4015, 4021, 4033, 4099, 4103, 4105, 4115, 4133, 4145, 4147, 4195, 4197, 4215, 4227, 4229, 4249, 4305, 4307, 4313, 4333, 4401, 4403, 4439, 4443, 4503, 4505, 4515, 4517, 4529, 4531, 4537, 4553, 4577, 4609, 4613, 4619, 4637, 4641, 4679, 4743, 4753, 4755, 4779, 4791, 4809, 4811, 4823, 4837, 4857, 4923, 4933, 4939, 4941, 4959, 5009, 5011, 5025, 5047, 5057, 5059, 5083, 5095, 5105, 5109, 5123, 5127, 5133, 5137, 5143, 5157, 5217, 5219, 5221, 5237, 5263, 5313, 5315, 5319, 5353, 5363, 5411, 5413, 5419, 5425, 5427, 5429, 5441, 5463, 5485, 5499, 5519, 5527, 5529, 5549, 5565, 5619, 5665, 5689, 5719, 5721, 5751, 5755, 5767, 5789, 5821, 5833, 5849, 5853, 5921, 5955, 5965, 5969, 5993, 6019, 6021, 6023, 6033, 6035, 6037, 6051, 6071, 6095, 6125, 6129, 6149, 6155, 6159, 6221, 6225, 6259, 6273, 6321, 6323, 6329, 6337, 6339, 6341, 6349, 6397, 6409, 6429, 6435, 6437, 6459, 6461, 6511, 6529, 6543, 6567, 6573, 6611, 6629, 6643, 6659, 6665, 6735, 6743, 6755, 6757, 6779, 6829, 6845, 6865, 6867, 6877, 6903, 6925, 6929, 6935, 6945, 6957, 6961, 6963, 6979, 7037, 7043, 7047, 7063, 7065, 7083, 7133, 7169, 7183, 7207, 7233, 7235, 7245, 7249, 7261, 7283, 7307, 7339, 7345, 7347, 7369, 7371, 7387, 7439, 7455, 7485, 7521, 7539, 7543, 7553, 7569, 7571, 7611, 7623, 7643, 7665, 7669, 7689, 7741, 7777, 7789, 7829, 7843, 7853, 7861, 7873, 7879, 7889, 7891, 7915, 7947, 7949, 7953, 7971, 7973, 7979, 7983, 7993, 8045, 8093, 8141, 8145, 8149, 8157, 8163])

    R = []
    for r in Residues:
        R += r
    R = set(R)
    D = {r:[] for r in R}
    for j, r in enumerate(Residues):
        for i in r:
            D[i].append(j + 1)

    for r in R:
        print('{}: {}'.format(r, D[r]))
        plt.plot([r] * len(D[r]), D[r], 'b.-')
    plt.show()


    for i in range(len(Residues) - 1):
        for j in range(i + 1, len(Residues)):
            r1 = Residues[i]
            r2 = Residues[j]
            print('a = {} and a = {} share {}'.format(i + 1, j + 1, sorted(list(set(r1).intersection(set(r2))))))

    exit()'''

    '''R = {(1,1): 0, (1,2): -1, (2,1): 2, (2,2): 1}
    S = {}
    residues = []
    a = 3
    comb = list(product([1,2], repeat=2*a))
    for c in comb:
        s = 0
        if c[0] == 1:
            s += 2
        else:
            s += 1
        for i in range(len(c) - 1):
            s += R[tuple(c[i:i + 2])]
        if not s == a:
            continue
        if c[-1] == 1:
            continue
        S[c] = s
    for c in S:
        print('{}: {}'.format(c, S[c]))
    exit()'''

    #find all 1,2 sequences giving us l_(2a-2) >= a
    '''R = {(1,1): 0, (1,2): -1, (2,1): 2, (2,2): 1}
    S = {}
    residues = []
    a = 6
    comb = list(product([1,2], repeat=2*a))
    for c in comb:
        s = 0
        if c[0] == 1:
            s += 2
        else:
            s += 1
        for i in range(len(c) - 1):
            s += R[tuple(c[i:i + 2])]
        if not s == a:
            continue
        if c[-1] == 1:
            continue
        S[c] = s
    V = []
    for c in S:
        #if c[:6] != (2,2,1,1,1,2):
        #    continue
        #print('{}'.format(c))
        #find matching w
        v = None
        #for w in range(1, 10**3, 2):
        for w in range(1, 2**(2*a + 1), 2):
        #for w in range(5, 6, 2):
            O = collatz_T(w)
            code = O['code']
            O = O['values']
            O = [1 if o % 4 in [0, 1] else 2 for o in O]
            if w == 1:
                O = [1 + (i % 2) for i in range(2*a+1)]
            if tuple(O[:len(c)]) == c:
                #print('{}'.format(c))
                if w < 2**(2*a + 1):
                    residues.append(w)
                if v is None:
                    v = get_Idw(2*a, w)
                rv = r_func(v)
                r = rv + 3**a * w
                r //= 4**a
                V.append((w, r))
                if in_Rn(w * (4**a - 3**a), a):
                    print('  w = {} | {} | rem: {} ***'.format(w, v, r))
                    pass
                elif in_Rn(w, a):
                    print('  w = {} | {} | rem: {} *'.format(w, v, r))
                    pass
                else:
                    print('  w = {} | {} | rem: {}'.format(w, v, r))
                    pass
                #print('  w = {}'.format(w))
    m = 3**a / 4**a
    V = sorted(V, key=lambda x:x[0])
    for v in V:
        #if v[0] % 3 == 1:
        #print('{}: (~= {} mod 3) | {} ({})'.format(v, v[0] % 3, int(m * v[0]), v[1] - int(m * v[0])))
        pass
    print(4**a/3**a)
    print('{}/{} comb are valid'.format(len(S), 2**(2*a)))
    residues = list(sorted(residues))
    print('residues: {} mod {}'.format(residues, 2**(2*a + 1)))

    #plt.plot([v[0] for v in V if v[0] % 3 == 0], [v[1] for v in V if v[0] % 3 == 0], 'r', label='0')
    #plt.plot([v[0] for v in V if v[0] % 3 == 1], [v[1] for v in V if v[0] % 3 == 1], 'g', label='1')
    #plt.plot([v[0] for v in V if v[0] % 3 == 2], [v[1] for v in V if v[0] % 3 == 2], 'b', label='2')
    #plt.plot([0, V[-1][0]], [3**a/4**a * 0, 3**a/4**a *V[-1][0]], 'k', label='ref')
    #plt.legend()
    #plt.show()
    exit()'''


    #find all 1,2 sequences giving us l_(2a-2) >= a
    R = {(1,1): 0, (1,2): -1, (2,1): 2, (2,2): 1}
    C = []
    a = 2
    comb = list(product([1,2], repeat=2*a))
    for c in comb:
        if c[-1] == 1:
            continue
        s = 0
        if c[0] == 1:
            s += 2
        else:
            s += 1
        for i in range(len(c) - 1):
            s += R[tuple(c[i:i + 2])]
        if not s == a:
            continue
        C.append(c)
    W = []
    for w in range(1, 2**(2*a + 1), 2):
    #for w in range(5, 6, 2):
        O = collatz_T(w)
        code = O['code']
        O = O['values']
        O = [1 if o % 4 in [0, 1] else 2 for o in O]
        if w == 1:
            O = [1 + (i % 2) for i in range(2*a+1)]
        if len(O) < 2*a:
            O += [1 + ((i + 1) % 2) for i in range(2*a - len(O))]
        O = tuple(O[:2*a])
        if O in C:
            W.append(w)
            print('{}: {}'.format(w, O))
    print(W)
    print(len(W))

    exit()

    #check multiple d
    w = 43
    starting_with = tuple()
    O = collatz_T(w)['values']
    types = [1 if (o % 4) in [0, 1] else 2 for o in O]
    print('w = {} | {}'.format(w, types))
    print('  Should stabalize at {}'.format(len(types) - 1))
    a_offset = 0
    check_d = 1
    V = []
    while check_d < 20:
        a = max(5, check_d - 2) + a_offset
        R, Is = get_all_r(a, 2*a + 1, include_even=False, starting_with=starting_with)
        L = []
        for r, I in zip(R, Is):
            r_ = r + 3**a * w
            d = nu(r, 2)
            d_ = nu(r_, 2)
            #d = nu(r, 4)
            #d_ = nu(r_, 4)
            if d_ == check_d:
                L.append(I)
        if len(L) in [0, 1]:
            a_offset += 1
            continue
        a_offset = 0
        matches_up_to = -1
        for i in range(a):
            match = True
            v = L[0][i]
            for I in L:
                if I[i] != v:
                    match = False
                    break
            if match:
                matches_up_to += 1
            else:
                break
        type = -1

        if L[0][matches_up_to + 1] == L[0][matches_up_to] + 1:
            type = 1
        else:
            type = 2
            k = L[0][matches_up_to + 1]
        if type == 1:
            starting_with = L[0][:matches_up_to]
            print('{} ({}): {}, ...)'.format(check_d, matches_up_to + 1, str(L[0][:matches_up_to + 1])[:-1]))
        else:
            starting_with = L[0][:matches_up_to + 1]
            print('{} ({}): {}, k, ...) : k >= {}'.format(check_d, matches_up_to + 1, str(L[0][:matches_up_to + 1])[:-1], k))
        V.append((check_d, matches_up_to + 1))
        check_d += 1
    print(V)
    exit()
