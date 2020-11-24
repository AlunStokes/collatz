from test_r import *

def collatz_accel(x):
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

if __name__ == '__main__':
    '''A = range(1, 16)
    for a in A:
        W = range(1, 1000)
        for w in W:
            V = range(1, 1000)
            for v in V:

                if w * 2**a - v * 3 == 1:
                    print((a, w, v), v / w)'''

    '''K = range(1, 16)
    for k in K:
        W = range(1, 1000)
        for w in W:
            #V = collatz_accel(w)['values']
            V = range(1, 1000)
            for v in V:
                if w * 2**k == 1 + 3*v:
                    print((k, w, v), prime_fac(w), prime_fac(3*v + 1))'''

    '''K = range(1, 7)
    for k in K:
        print(k)
        W = range(1, 30)
        for w in W:
            if (2**k * w - 1) % 3 == 0:
                print((k, w))
        print('')'''

    N = range(1, 10)
    for n in N:
        K = range(0, 10)
        for k in K:
            M = range(1, n + 1)
            for m in M:
                W = range(1, 20)
                for w in W:
                    r = 2**k * w - 1
                    if r % 3 == 0:
                        r //= 3
                        r *= 2**n
                        r -= 3**m * w
                        if in_Rn(r, m):
                        #if 3 not in prime_fac(r):
                            print((n,k,m,w), r, inv_r(r, m))


    '''X = range(0, 50000)
    for x in X:
        w = 3*x + 2
        v = 2*x + 1
        if v in collatz_accel(w)['values']:
            print(x)
            print((w, v))
            print(collatz_accel(w)['values'])'''

    '''P = range(3, 150)
    for p in P:
        C = collatz_accel(p)['values']
        i = 1
        while i < len(C):
            if int(abs(C[i] - C[0])) in {2,4,8,16,32,64,128,256,512,1024,2048}:
                print(p)
                break
            i += 1'''

    '''L = []
    A = range(2, 30)
    for a in A:
        K = range(1, 30)
        for k in K:
            S = [2*a - k - (a - 2)]
            i = 1
            while i < a:
                S.append(1)
                i += 1
            L.append(S)

    R = []
    for a in L:
        x0 = 1
        x = x0
        for i in a:
            x = (3 * x + 1) / 2**i
        R.append((a, x / x0))

    for r in R:
        print(r)'''
