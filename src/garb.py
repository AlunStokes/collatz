from math import sqrt, log, log2, ceil
from test_r import *
from scaling_factors import *
from utils import *
import itertools
from self_ref import get_self_cont_multi, get_pseudo_self_cont_multi, get_neg_inv_multi
import numpy as np

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
        return None
        #raise Exception('modular inverse does not exist')
    else:
        return x % m

def all_equal(A, n):
    for a in A:
        if a != n:
            return False
    return True

def find_lines(R):
    R = R.copy()
    L = {}
    while len(R) > 0:
        r = R[0]
        del R[0]
        p = r[1] - 2*r[0]
        L[p] = [r]
        R_ = R.copy()
        i = 0
        while i < len(R):
            r_ = R[i]
            if 2*r_[0] + p == r_[1]:
                L[p].append(r_)
                del R[i]
                continue
            i += 1
    return L

def get_invs(w, U):
    I = []
    M = range(0, U, 1)
    for m in M:
        r = 2**m % w
        inv = mod_inv(3, w)
        r *= inv
        r %= w
        r = w - r
        r %= w
        if r not in I:
            I.append(r)
    I = sorted(I)
    return I

def get_as_congruence_class(n, w):
    r = n % w
    q = (n - r) // w
    return q, r

#min n so that 3^a divides 2^b*n - r
def get_min_div(a, b, r):
    n = 0
    while (2**b*n - r) % 3**a != 0:
        n += 1
    return n

def nu(x, b=2):
    k = 1
    while x % b**k == 0:
        k += 1
    return k - 1

def largest_pow_less(x, b=2):
    k = 1
    while x >= b**k:
        k += 1
    return k - 1

def find_ab(w, U):
    L = []
    B = range(1, U)
    for b in B:
        A = range(1, b+1)
        for a in A:
            if in_Rn(w*(2**b - 3**a), a):
                L.append((a, b))
    return L

if __name__ == '__main__':

    #Finds what happens when p processed applied to p in Pw (Important, shows never leads to odd number for appropriate p)
    '''d = 6
    P = [(4**k - 4)//d for k in range(1, 10)]
    print(P)
    for p in P:
        r = p
        L1 = [r]
        i = 0
        while r % 2 == 0 and r > 0:
            if i % 2 == 0:
                r //= 2
                r -= 1
            else:
                r //= 2
            i += 1
            L1.append(r)
        r = p
        L2 = [r]
        i = 0
        while r % 2 == 0 and r > 0:
            if i % 2 == 1:
                r //= 2
                r -= 1
            else:
                r //= 2
            i += 1
            L2.append(r)
        print('{}'.format(p))
        print('{}'.format(L1))
        print('{}'.format(L2))
        print('')

    exit()'''

    #Check largest power of 2 for w(2^2a-3^a)
    '''Ws = {}
    for a in range(50, 51):
        print('a = {}'.format(a))
        b = 2*a
        r1 = 2**b - 3**a
        I1 = inv_r(r1, a)
        W = []
        D = []
        for w in range(1, min(1000, 3**a), 2):
            rw = w * (2**b - 3**a)
            if in_Rn(rw, a):
                W.append(w)
                Iw = inv_r(rw, a)
                d = Iw[-1] - I1[-1]
                D.append(d)
                print('{}: {}'.format(w, d))
        Ws[a] = W
        print('')'''

    '''w = 565
    a = 8
    m = a - 4
    p = tuple([i - 1 for i in inv_r(w * (4**a - 3**a), a)])
    p1 = tuple([i - 2*a for i in p[m:]])
    p2 = p[:m]

    print('w = {}; a = {}; m = {}'.format(w, a, m))
    print('p1 = {}'.format(p1))
    print('p2 = {}'.format(p2))

    r1 = r_func([i + 1 for i in p1])
    r2 = r_func([i + 1 for i in p2])

    if not 3**(a - m) * (r2 + 3**m * w) % 4**a == 0:
        print('4^{} does not divide R(p2) + 3^{}*w'.format(a, m))

    r = 3**(a - m) * (r2 + 3**m * w) // 4**a
    q = w - r1
    print('{} = {}'.format(r, q))
    exit()'''

    w = 565
    a = 7
    print('w = {}'.format(w))
    R = get_all_r(a, 2*a + 9)
    for r in R:
        if (r + 3**a * w) / (4**a * w) > 1:
            print('{}: {}'.format(r, tuple([i - 1 for i in inv_r(r, a)])))
    exit()


    #Check largest power of 2 for w(2^b-3^a)
    '''for b in range(1, 20):
        for a in range(1, b + 1):
            if 2**b - 3**a < 0:
                break
            W = []
            D = []
            for w in range(1, min(100000, 3**a), 2):
                rw = w * (2**b - 3**a)
                if in_Rn(rw, a):
                    W.append(w)
                    Iw = inv_r(rw, a)
                    d = Iw[-1] - (b)
                    D.append(d)
                    if d < 3:
                        print('{} | {}: {}'.format((a, b), w, d))
                        print('\t{}'.format(Iw))
                        print()
            if b == 2*a:
                plt.plot(np.log(W), D, 'r.', markersize=4)
                #plt.plot(W, D, 'r.', markersize=4)
            else:
                plt.plot(np.log(W), D, 'b.', markersize=4)
                #plt.plot(W, D, 'b.', markersize=4)
    #plt.plot(W, D, '.')
    #plt.ylim(0)
    plt.ylabel('v_a - (b - 1)')
    plt.xlabel('log(w)')
    plt.show()'''

    #Largest powers of 2 less than r
    '''a = 6
    b = 2*a
    W = get_wab(a, b, 1000)

    for w in W:
        r = w * (2**b - 3**a)
        I = inv_r(r, a)
        P = []
        k = a - 1
        while k >= 0:
            P.append(largest_pow_less(r//3**k + 1, 2))
            k -= 1
        if P[-1] > I[-1]:
            print('{} : {} | {} ***'.format(r, P, I))
        else:
            print('{} : {} | {}'.format(r, P, I))'''

    #Check largest power of 2 for multiples in R
    '''a = 8
    b = 20
    R = get_all_r(a, b)
    R = [r for r in R if r % 2 == 1]
    R = sorted(R)

    M = []
    i = 0
    while i < len(R) - 1:
        j = i + 1
        while j < len(R):
            if R[i] % R[j] == 0:
                M.append((R[j], R[i], R[i] // R[j]))
            elif R[j] % R[i] == 0:
                M.append((R[i], R[j], R[j] // R[i]))
            j += 1
        i += 1

    for m in M:
        r1 = m[0]
        r2 = m[1]
        q = m[2]
        I1 = inv_r(r1, a)
        I2 = inv_r(r2, a)
        d = I2[-1] - I1[-1]
        plt.plot(q, d, 'b.')
        if d < 0:
            print('{}: {}'.format(m, d))
            print('\t{}'.format(I1))
            print('\t{}'.format(I2))
    plt.show()'''

    '''W = get_pseudo_self_cont_multi(10**9, workers=8)
    print(W)
    exit()'''

    #Checking the form of pseudo-self-contained numbers
    '''W = {}
    w = 5
    while w < 10**4:
        if w % 3 == 0:
            w += 2
            continue
        r = ((w % 3) * w - 1) // 3
        O = collatz(w)['values'][1:]
        for o in O:
            if o % w == r:
                #print('{} ({}): {} = {} * w + {}'.format(w, r, o, *get_as_congruence_class(o, w)))
                if w in W:
                    W[w].append(o)
                else:
                    W[w] = [o]
        w += 2

    for w in W:
        r = ((w % 3) * w - 1) // 3
        r2 = ((2*w % 3) * 2*w - 1) // 3
        print('w = {} ({} | {})'.format(w, r, r2))
        for v in W[w]:
            print('\t{} = {}*w + {} | '.format(v, *get_as_congruence_class(v, w)))

    print(W)

    exit()'''

    #Find orbit length distribution
    '''W = []
    L = []
    w = 5
    while w < 10**4:
        O = collatz(w)['values'][1:]
        W.append(w)
        L.append(len(O))
        w += 1
    R = [l / w for l, w in zip(L, W)]

    k = 50
    R_avg = [sum(R[i: i + k])/k for i in range(len(R) - k)]
    print(R_avg[-5:-1])
    plt.plot(R_avg)
    #plt.ylim(0, 0.5)
    plt.show()


    exit()'''


    #Find distribution of -k^-1 hit by orbits
    '''#How many times each is hit
    D = {}
    #How many eligible to be hit
    E = {}
    #Keeping track of w and orbit length
    w = 5
    a = 6
    #U = 500 + 1
    #U = 10**a
    #U = 10**a // 4
    while w < 10**a:
        K = []
        k = 2
        while k < w:
            #print('({},{}) = {}'.format(w, k, gcd))
            if egcd(w, k)[0] == 1:
                if k in E:
                    E[k] += 1
                else:
                    E[k] = 1
                K.append(k)
            k += 1
        V = [w - mod_inv(k, w) for k in K]
        O = collatz(w)['values'][1:]
        Ow = [o % w for o in O]
        for k, v in zip(K,V):
            if v in Ow:
                if k in D:
                    #D[k] += Ow.count(v)
                    D[k] += 1
                else:
                    #D[k] = Ow.count(v)
                    D[k] = 1
        w += 1

        #Graph every so often
        if w % 500 == 0:
            keysD = sorted(list(D.keys()))
            R_inverses = []
            V = []
            for k in keysD[:w // 10]:
                #R_inverses.append(D[k] / w)
                R_inverses.append(D[k] / E[k])
                if k == 3:
                    plt.plot([k], [R_inverses[-1]], 'r.')
                else:
                    plt.plot([k], [R_inverses[-1]], 'b.')
                #V.append((k, R_inverses[-1]))

            #Plot in order
            #V = sorted(V, key=lambda x:x[1])
            #plt.plot([v[1] for v in V], 'b.')
            #ind = [v[0] for v in V].index(3)
            #plt.plot([ind], [V[ind]], 'r.')


            #R_orblen = [l / w for l, w in zip(L, W)]
            #smoothing = 50
            #R_orblen_avg = [sum(R_orblen[i: i + smoothing])/smoothing for i in range(len(R_orblen) - smoothing)]
            #avg_len = R_orblen_avg[-1]
            #avg_inv = sum(R_inverses) / len(R_inverses)
            M = max(R_inverses)
            m = min(R_inverses)
            m2 = min(set(R_inverses).difference({m}))
            print('Estimated density of -3^(-1) is {:.4f}'.format(3 / 2 * len(get_pseudo_self_cont_multi(w)) / w))
            print('Calculated density of -3^(-1) is {:.4f}'.format(D[3] / E[3]))
            print('Maximum density is {:.4f} ({:.2f}) for -{}^(-1)'.format(M, M / m, keysD[R_inverses.index(M)]))
            print('Minimum density is {:.4f} for -{}^(-1)'.format(m, keysD[R_inverses.index(m)]))
            print('2nd to minimum density is {:.4f} ({:.2f}) for -{}^(-1)'.format(m2, m2 / m, keysD[R_inverses.index(m2)]))
            #print('Average density by orbit is {:.4f}'.format(avg_len))
            #print('Average density by inverses is {:.4f}'.format(avg_inv))
            print()

            plt.ylim(0 - M/10, M + M/10)
            plt.grid()
            plt.xlabel('k')
            plt.ylabel('prop. of w whose orbit contains -k^(-1)')
            #plt.plot([keysD[0], keysD[-1]], [avg_len, avg_len], 'r', label='{}'.format('avg by orbit len'))
            #plt.plot([keysD[0], keysD[-1]], [avg_inv, avg_inv], 'm', label='{}'.format('avg by inv data'))
            #plt.legend()
            #plt.savefig('./images/neg_inv_dist/no_multiplicity_all/{}.png'.format(w))
            #plt.savefig('./images/neg_inv_dist/multiplicity_all/{}.png'.format(w))
            plt.savefig('./images/neg_inv_dist/no_multiplicity/{}.png'.format(w))
            #plt.savefig('./images/neg_inv_dist/multiplicity/{}.png'.format(w))
            plt.clf()

    exit()'''


    #How many times each is hit
    '''P = primes_less_than(10**3)
    Vs = []
    D = {}
    #How many eligible to be hit
    E = {}
    #Keeping track of w and orbit length
    w = 5
    a = 8
    #U = 500 + 1
    #U = 10**a
    #U = 10**a // 4
    while w < 10**a:
        K = []
        for k in P:
            #print('({},{}) = {}'.format(w, k, gcd))
            if w % k != 0:
                if k in E:
                    E[k] += 1
                else:
                    E[k] = 1
                K.append(k)
            k += 1
        V = [w - mod_inv(k, w) for k in K]
        O = collatz_accel(w)['values'][1:]
        Ow = [o % w for o in O]
        for k, v in zip(K,V):
            if v in Ow:
                if k in D:
                    #D[k] += Ow.count(v)
                    D[k] += 1
                else:
                    #D[k] = Ow.count(v)
                    D[k] = 1
        w += 1

        #Graph every so often
        if w % 1000 == 0:
            keysD = sorted(list(D.keys()))
            R_inverses = []
            V = []
            for k in keysD[:w // 10]:
                #R_inverses.append(D[k] / w)
                R_inverses.append(D[k] / E[k])
                #R_inverses.append(D[k] / ((w - k) * (1 - 1/k)))
                #if k == 3:
                #    plt.plot([k], [R_inverses[-1]], 'r.')
                #else:
                #    plt.plot([k], [R_inverses[-1]], 'b.')
                V.append((k, R_inverses[-1]))
            Vs.append(V)

            Y = []
            plt.rcParams["figure.figsize"] = (12,8)
            for v in V:
                k = v[0]
                y = v[1]
                if k == 3:
                    plt.plot([k], [y], 'r.')
                else:
                    plt.plot([k], [y], 'b.')

            print(w)

            M = max([v[1] for v in V])
            #m = min([v[1] for v in Vs[i - 1]])

            #plt.ylim(0 - M/10, max(M + M/10, 0.025))
            #plt.ylim(0 - M/10, M + M/10)
            plt.ylim(0, M + M/10)
            plt.grid()
            plt.xlabel('k')
            plt.ylabel('prop. of w whose orbit contains -k^(-1)')

            plt.savefig('./images/neg_inv_dist/no_multiplicity_p/{}.png'.format(w))
            plt.clf()

    exit()'''

    '''L = 10**5
    w = L
    while w < L + 20:
        if w % 3 == 0:
            w += 1
            continue
        O = collatz(w)['values']
        Ow = [o % w for o in O]
        Ow = sorted(Ow)
        r = ((w % 3)*w - 1)//3
        print('{} ({}): {}'.format(w, r, Ow))
        print()
        w += 1'''


    #Find distribution of residue classes
    '''#How many times each is hit
    D = {}
    #How many eligible to be hit
    E = {}
    #Keeping track of w and orbit length
    w = 5
    a = 4
    #U = 100 + 1
    #Keep track of -3^(-1) proportion
    T = []
    U = 10**a
    while w < 10**a:
        V = list(range(0, U))
        for v in V:
            if v in E:
                E[v] += 1
            else:
                E[v] = 1
        O = collatz(w)['values'][1:]
        Ow = [o % w for o in O]
        for v in V:
            if v in Ow:
                if v in D:
                    D[v] += Ow.count(v)
                    #D[k] += 1
                else:
                    D[v] = Ow.count(v)
                    #D[k] = 1
        if w % 3 != 0:
            r = ((w % 3) * w - 1) // 3
            if r in Ow:
                T.append(Ow.count(r))
                #T.append(1)
            else:
                T.append(0)

        w += 1

        #Graph every so often
        if w % 10 == 0:
            keysD = sorted(list(D.keys()))
            R = []
            for k in keysD[:U]:
                R.append(D[k] / E[k])
                plt.plot([k], [R[-1]], 'b.')
            M = max(R)
            m = min(R)
            m2 = min(set(R).difference({m}))
            print('Estimated density of -3^(-1) is {:.4f}'.format(3 / 2 * len(get_pseudo_self_cont_multi(w)) / w))
            print('Calculated density of -3^(-1) is {:.4f}'.format(D[3] / E[3]))
            print('Maximum density is {:.4f} ({:.2f}) for {}'.format(M, M / m, keysD[R.index(M)]))
            print('Minimum density is {:.4f} for {}'.format(m, keysD[R.index(m)]))
            print('2nd to minimum density is {:.4f} ({:.2f}) for {}'.format(m2, m2 / m, keysD[R.index(m2)]))
            print()

            plt.ylim(0 - M/10, M + M/10)
            t_avg = sum(T) / len(T)
            plt.plot([keysD[0], keysD[-1]], [t_avg, t_avg], 'r')
            #plt.savefig('./images/residue_dist/no_multiplicity/{}.png'.format(w))
            plt.savefig('./images/residue_dist/multiplicity/{}.png'.format(w))
            plt.clf()'''

    #Look at which -k^{-1} in orbit as w grows
    '''upper_bound = 10**3

    w = 5
    while w < upper_bound:
        E = []
        i = 2
        while i < w:
            if egcd(i, w)[0] == 1:
                E.append(i)
            i += 1
        V = [w - mod_inv(e, w) for e in E]
        for e,v in zip(E, V):
            print('{} -> {}'.format(e, v))
        O = collatz(w)['values'][1:]
        Ow = [o % w for o in O]
        T = list(sorted(list(set(Ow).intersection(set(E)))))
        if 3 in T:
            plt.plot(Ow)
            plt.show()
            #print('w = {} ({})'.format(w, E))
            print('w = {}'.format(w))
            print('  {}'.format(T))

        w += 1'''

    #The next couple have to do with -k^-1 stopping showing up
    '''k = 347
    U = 10**3
    w = 5
    while w < 10**5:
        E = []
        for k in range(2, min(U, w)):
            if egcd(k, w)[0] == 1:
                E.append(k)
        V = [w - mod_inv(k, w) for k in E]
        VE = {}
        for v,e in zip(V, E):
            VE[v] = e
        O = collatz(w)['values']
        Ow = [o % w for o in O]
        T = set(Ow).intersection(set(V))
        T = [VE[v] for v in T]
        T = list(sorted(list(T)))
        if k in T:
            print('{}: {}'.format(w, T[:20]))
        w += 1
    exit()'''


    '''K = [293]
    a = 7
    U = 10**a
    print('Up to 10^{}'.format(a))
    for k in K:
        W = get_neg_inv_multi(U, k)
        print('k = {}'.format(k))
        #print('{}: {}'.format(k, W))
        print(W[:5])
        print(W[-5:])
        print(len(W))
        print(len(W) / (U * (1 - 1/k)))
        print()
    exit()'''

    '''k = 347
    U = 10**8
    W = get_neg_inv_multi(U, k)
    print(W)
    print(len(W))
    print(len(W) / (U / k))
    plt.plot(np.log(W), label='{}'.format(k))
    plt.legend()
    plt.show()'''


    #Finding pseudo-self-contained numbers
    '''#W = get_pseudo_self_cont_multi(10**7, workers=8)
    W = [5, 7, 11, 13, 19, 25, 31, 43, 47, 49, 55, 61, 62, 67, 71, 73, 79, 83, 95, 109, 121, 125, 133, 137, 145, 157, 166, 175, 259, 265, 275, 293, 319, 335, 347, 377, 553, 586, 619, 671, 694, 733, 763, 781, 823, 925, 967, 1001, 1342, 1579, 1951, 2005, 2111, 2371, 2599, 2684, 2813, 2929, 3463, 3901, 4315, 5197, 5341, 5467, 6925, 6943, 7157, 8203, 8617, 9715, 10379, 12343, 19151, 21865, 22369, 26899, 35327, 36479, 38302, 39071, 44903, 58025, 106495, 189097, 205135, 470503, 1423699, 2025797, 4051594]
    Wsc = [31, 62, 83, 166, 293, 347, 586, 671, 694, 1342, 2684, 19151, 38302, 2025797, 4051594]
    #W = [w for w in W if w not in Wsc]
    print(W)
    print(len(W))

    W_mult = []
    for w in W:
        K = [0]
        k = 1
        while k < 100:
            wk = 2**k * w
            r = ((wk % 3) * wk - 1) // 3
            O = collatz(wk)['values']
            Ow = [o % wk for o in O]
            if r in Ow:
                K.append(k)
            k += 1
        W_mult.append((w, K))

    for w in W_mult:
        if len(w[1]) > 1:
            print(w)

    m = {}
    m[1] = 0
    m[2] = 0
    for w in W:
        m[w %  3] += 1
    print('{}/{} ~= {:.4f} are equivalent to 1 mod 3'.format(m[1], len(W), m[1] / len(W)))
    print('{}/{} ~= {:.4f} are equivalent to 1 mod 3'.format(m[2], len(W), m[2] / len(W)))

    for w in W:
        O = collatz(w)['values'][1:]
        r = ((w % 3) * w - 1) // 3
        Ow = [o % w for o in O]
        print('{} ~={}'.format(w, w % 3))
        for i, n in enumerate(Ow):
            if n == r:
                #o0 = O[i - 2]
                #q0, r0 = get_as_congruence_class(o0, w)
                #o1 = O[i - 1]
                #q1, r1 = get_as_congruence_class(o1, w)
                #o2 = O[i - 0]
                #q2, r2 = get_as_congruence_class(o2, w)
                #print('  {} = {}*w + {} ({}) -> {} = {}*w + {} ({}) - > {} = {}*w + {} ({})'.format(o0, q0, r0, mod_inv(w - r0, w), o1, q1, r1, mod_inv(w - r1, w), o2, q2, r2, mod_inv(w - r2, w)))
                Q = [get_as_congruence_class(O[i - k], w) for k in range(min(4, i), -1, -1)]
                Q = list(reversed(Q))
                print('  {}'.format(Q))
                plt.plot([q[0] for q in Q])
        print()
    plt.show()'''


    '''W_rand = np.random.randint(5, 10**6, 80)
    for w in W_rand:
        O = collatz(w)['values'][1:]
        ind = len(O) // 2
        Q = [get_as_congruence_class(O[ind - k], w) for k in range(min(4, ind), -1, -1)]
        print('  {}'.format(Q))
        plt.plot([q[0] for q in Q])
    plt.show()'''



    #Checking against inscreasing modulo
    '''S = {}
    S[1] = [(1 * (4**k - 1)) // 3 for k in range(0, 12)]
    S[2] = [(2 * (4**k - 1)) // 3 for k in range(0, 12)]
    W = []
    M = [5**0 * 3**1 * 2**k for k in range(0, 6)]
    T = {}
    T[1] = []
    T[2] = []
    w = 5
    while w < 10**6:
        if w % 3 == 0:
            w += 2
            continue
        d = w % 3
        O = collatz(w)['values'][1:]
        for n in O:
            q, r = get_as_congruence_class(n, w)
            #Checks 1 and 2
            #if r == (d*w - 1)//3 and n % 2 == 1:
            #Checks 1 and 3
            if r == (d*w - 1)//3 and q in S[d]:
            #if r == (d*w - 1)//3:
            #Checks 2 and 3
            #if n % 2 == 1 and q in S[d]:
                #print('{} (~={}, {} mod 3, 6): {} * w + {}'.format(w, w % 3, w % 6, q, r))
                #print('{} (~={}, {}, {} mod 3, 6, {}): {} * w + {}'.format(w, w % 3, w % 6, w % m, m, q, r))
                V = [w % m for m in M]
                follows_pattern = True
                i = 1
                while i < len(V):
                    diff = V[i] - V[i - 1]
                    if diff not in [0, M[i - 1]]:
                        follows_pattern = False
                        break
                    i += 1
                T[w % 3].append(q)
                if follows_pattern:
                    print('{} (~= {} mod {}): {} * w + {} ***'.format(w, V, M, q, r))
                else:
                    print('{} (~= {} mod {}): {} * w + {}'.format(w, V, M, q, r))
                W.append(w)

                #print('{}: C({}) = {}'.format(w, r, collatz(r)['values']))
        #print('')
        w += 2

    T[1] = sorted(list(set(T[1])))
    T[2] = sorted(list(set(T[2])))
    print(T[1])
    print(T[2])
    plt.plot(np.log(W))
    plt.show()
    exit()'''

    '''W = [31, 62, 83, 166, 293, 586, 347, 694, 671, 1342, 2684, 19151, 38302, 2025797, 4051594]
    W = [w for w in W if w % 2 == 1]
    for w in W:
        d = w % 3
        r1 = (d*w - 1) // 3
        if d == 1:
            r2 = 2 * r1
        else:
            r2 = (r1 - 1) / 2
        O = collatz(w)['values'][1:]
        Ow = [o % w for o in O]
        index = Ow.index(0)
        print('{}: {} ({}) -> {} ({}) -> {} ({})'.format(w, get_as_congruence_class(O[index - 2], w), r2, get_as_congruence_class(O[index - 1], w), r1, get_as_congruence_class(O[index], w), 0))
'''
    #Interesting distribution vs random
    '''fig = plt.figure(figsize=(16,4))
    #W = [31, 62, 83, 166, 293, 347, 586, 671, 694, 1342, 2684, 19151, 38302, 2025797, 4051594]
    W = [31,83]
    #W = [np.random.randint(1, 10**6) for _ in range(25)]
    W = [w for w in W if w % 2 == 1]
    for w in W:
        printed_line = False
        print('{}'.format(w))
        Q = []
        R = []
        O = collatz(w)['values'][1:]
        longest_num = len(str(max(O)))
        Ow = [o % w for o in O]
        d = w % 3
        S = [(d * (4**k - 1)) // 3 for k in range(0, 6)]
        for n in O:
            q, r = get_as_congruence_class(n, w)
            if n % w == 0:
                if not printed_line:
                    print('---------------------')
                    printed_line = True
            if n in S:
                #print('True  {} = {} * {} + {} ***'.format(n, q, w, r))
                print('* {} = {} * w + {}'.format(n, q, r))
                pass
            else:
                #print('False {} = {} * {} + {} ***'.format(n, q, w, r))
                print('  {} = {} * w + {}'.format(n, q, r))
                pass
            #print('{}{} = {} * w + {}'.format(n, ' ' * (longest_num - len(str(n))), q, r))
            Q.append(q)
            R.append(r)


        index = len(Ow) - list(reversed(Ow)).index(0)
        Q1 = Q[:index]
        R1 = R[:index]
        Q2 = Q[index:]
        R2 = R[index:]
        #print('{}'.format(R))
        print('--------------------------------------------')
        S = []
        for (q1,r1), (q2,r2) in zip(zip(Q1, R1), zip(Q2, R2)):
            n1 = q1 * w + r1
            n2 = q2 * w + r2
            s1 = '{} (~= {}) = {} * w + {}'.format(n1, n1 % 3, q1, r1)
            s2 = '{} (~= {}) = {} * w + {}'.format(n2, n2 % 3, q2, r2)
            S.append((s1, s2))
        max_len = max([len(s[0]) for s in S])
        for s in S:
            print(s[0] +  (max_len-len(s[0]) + 1)*" " + '|' + " " + s[1])
        #print(Q)
        Ri = [mod_inv(r, w) for r in O]
        Ri = [w - r if r is not None else r for r in Ri]
        #print(Ri)
        #print('2 * {} - 3 = {} = {}'.format(Ri[Ri.index(3) - 1], Ri[Ri.index(3) - 1] * 2 - 3, w))
        print([Q[i] - Q[i-1] for i in range(1, len(Q))])
        print('')
        #max_q = max(Q)
        #plt.plot([q / max_q for q in Q], label='q{}'.format(w))
        plt.plot(Q, label='q{}'.format(w))
        index = Ow.index(0)
        plt.plot(index, Q[index], 'k.')
        #plt.title('dividends for w = {}'.format(w))
        #plt.show()
        #max_r = max(R)
        #plt.plot([r / max_r for r in R], label='r{}'.format(w))
        #plt.plot(R, label='r{}'.format(w))
        #plt.title('remainders for w = {}'.format(w))
        #plt.show()
    plt.legend()
    plt.show()'''





    '''w = 1
    while w < 10**7:
        if w % 3 == 0:
            w += 2
            continue
        d = 3
        if w % 3 == 2:
            d = 6
        O = collatz_accel(w)['values']
        K = range(1, 10)
        for k in K:
            r = (4**k * w - d//3)//d
            if r in O:
                print((w, k))
                print('')
        w += 2'''

    '''w = 31

    B = range(1, 40)
    for b in B:
        A = range(1, b + 1)
        for a in A:
            K = range(1, b)
            for k in K:
                r = w * (2**(b + 2*k) - 3**(a + 1)) - 2**b
                if r % 3 != 0:
                    continue
                r //= 3
                if in_Rn(r, a):
                    print((a, b, k))
                    print(r)
                    I = inv_r(r, a)
                    print(I)
                    print('')'''




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
