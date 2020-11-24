from math import log2, log

import matplotlib.pyplot as plt

#2-adic valuation function on the natural numbers
def v2(x):
    if x == 0:
        return 0
    i = 1
    while True:
        if x % 2**i != 0:
            return i - 1
        i += 1

#determines if value r is in the image of R_a
def in_Rn(r, n):
    if n == 0 and r == 0:
        return True
    for i in range(n-1, 0, -1):
        r = r // 2**v2(r)
        r -= 3**i

    if 2**v2(r) != r:
        return False
    return True

#Gives the sequence A that generates r with length of A = n
def inv_r(r, n):
    L = []
    i = 0
    while i < n:
        s = 0
        j = 0
        while j < i:
            s += 3**(n-(j+1)) * 2**(L[j] - 1)
            j += 1
        L.append(v2(r - s) + 1)
        i += 1
    return tuple(L)

#Gets the set (given as a list) of scaling factors w given a for b = 2a
#includes even w if include_even = True
#prints if verbose = True
def get_wa(a, w_bound=False, include_even=False, verbose=False):
    return get_wab(a, 2*a, w_bound, include_even, verbose)

#Gets the set (given as a list) of scaling factors w given a, b
#includes even w if include_even = True
#prints if verbose = True
#prints out the following:
#w | w * r : A
#where A is the sequence that generates w * r
def get_wab(a, b, w_bound=False, include_even=False, verbose=False):
    wab = []

    if verbose:
        print((a, b))

    r = 2**b - 3**a
    W = range(1, 3**a)
    if w_bound:
        W = range(1, min(3**a, w_bound))
    for w in W:
        if in_Rn(w*r, a):
            if not include_even:
                if w % 2 == 1:
                    if verbose:
                        print('{} | {}: {}'.format(w, w * r, inv_r(w*r, a)))
                    wab.append(w)
            else:
                if verbose:
                    print('{} | {}: {}'.format(w, w * r, inv_r(w*r, a)))
                wab.append(w)

    return wab

def seq_to_diff(s):
    N = []
    i = 1
    while i < len(s):
        N.append(s[i] - s[i - 1])
        i += 1
    return tuple(N)

def subtract_tuple(a, b):
    return tuple([i - j for i,j in zip(a, b)])

if __name__ == '__main__':

    #A is the set of lengths of |A| we are looking at scaling factors of 4^a - 3^a for
    A = range(3, 30)
    #W holds the set of w such that w(4^a-3^a) in Image(Ra) for each a, up to w = 2000
    W = []
    for a in A:
        W.append(get_wa(a, 100))

    #P maps a scaling factor, w, to a list of tuples, where each tuple is the sequence that generates w(4^a - 3^a)
    P = {}
    for w, a in zip(W, A):
        #This is the generating sequence for 4^a - 3^a
        ref_tup = tuple(range(1, 2*a + 1, 2))
        for i in w:
            if i not in P:
                P[i] = []
            #P[i].append(subtract_tuple(seq_to_diff(inv_r(i * (4**a - 3**a), a)), seq_to_diff(ref_tup)))
            P[i].append(seq_to_diff(inv_r(i * (4**a - 3**a), a)))

    #This prints out the elements of P in a formatted way
    for p in P:
        print('{}'.format(p))
        if len(P[p]) == 1:
            pass
        elif len(P[p][0]) + 1 == len(P[p][1]):
            pass
            #Continue should be uncommneted to show only the anomolies where there is a skip between the first and
            #second occurance of w being a scaling factor
            #Note that other anomolies will show up (2 skips. for example)
            #continue
        for i in P[p]:
            print('\t{}: {}'.format(len(i) + 1, i))
