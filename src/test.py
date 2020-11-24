from test_r import *
from utils import *

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

def collatz_b(x, b):
    for i in b:
        if i == 0:
            x /= 2
        else:
            x *= 3
            x += 1
            x /= 2
    return x


if __name__ == '__main__':

    working = []

    B = range(1, 11)
    for b in B:
        A = range(0, b + 1)
        for a in A:
            W = range(0, 3**a)
            for w in W:
                R = range(3**a - 2**a, int(2**b * ((3/2)**a - 1)) + 1)
                for r in R:
                    kt = (w * (3**a - 2**b) + r)
                    kb = (3**a * (2**b - 3**a))
                    if kt % kb == 0:

                        if not in_Rn(r, a):
                            #print('{} fails by r not in Im(R_{})'.format((a,b,w,r), a))
                            pass
                            continue
                        inv = inv_r(r, a)
                        if len(inv) > 0:
                            m = max(inv)
                        else:
                            m = 0
                        if m > b:
                            #print('{} fails by max(A) > |b| (max(A) is {})'.format((a,b,w,r), m))
                            pass
                            continue
                        k = kt // kb
                        if k < 0:
                            #print('{} fails by k < 0 (k = {})'.format((a,b,w,r), k))
                            pass
                            continue
                        n = smallest_n(3**a, 2**b, -r)
                        if not w == n:
                            #print('{} fails by w not right (correct is {})'.format((a,b,w,r), n))
                            pass
                            continue

                        working.append((a,b,w,r))

                        t = (2**b*w - r) // 3**a
                        print('{}'.format((a,b,w,r)))
                        #print('|A|: {}'.format(a))
                        #print('|b|: {}'.format(b))
                        #print('w: {}'.format(w))
                        #print('r: {}'.format(r))
                        #print('k: {}'.format(k))
                        print('x = {}'.format(2**b * k + t))
                        print('T(x) = {}'.format(3**a * k + w))
                        #print('r = {}'.format(r))
                        #print('w = {}'.format(w))
                        print('max(A) = {}'.format(m))
                        print('A = {}'.format(inv_r(r, a)))
                        #print('r / (2^b - 3^a) = {}'.format(r // (2**b - 3**a)))
                        print('')

    print('')
    for w in working:
        print(w)
 
