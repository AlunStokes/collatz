from utils import *
from check_in_r import *

def is_power_of_2(x):
    if 2**v2(x) == x:
        return True
    return False

def print_eq_default(a, b, w, wp, iter):
    s = '{}*2^{} - {}*3^{} = '.format(w, b, wp, a)
    if a == 1:
        s += '2^m{}'.format(a + iter)
    else:
        i = iter + 1
        while i < a + iter:
            s += '2^m{} * ('.format(i + 1)
            i += 1
        s = s[:-4]

        i = 1
        while i < a:
            s += ' + 3^{})'.format(i)
            i += 1
        s = s[:-1]

    print(s)

def print_eq_r2(a, w, wp, r2, iter):
    s = '{} - {}*3^{}*2^{} = '.format(w, wp, a, r2)
    if a == 1:
        s += '2^m{}'.format(a + iter)
    else:
        i = iter + 1
        while i < a + iter:
            s += '2^m{} * ('.format(i + 1)
            i += 1
        s = s[:-4]

        i = 1
        while i < a:
            s += ' + 3^{})'.format(i)
            i += 1
        s = s[:-1]

    print(s)

def solve(a, b, w, wp, r2 = 0, iter = 0):
    if iter == 0:
        print_eq_default(a, b, w, wp, iter)
    if r2 == 0:
        d = v2(3*wp + 1)
        if d < b:
            b -= d
            wp = (3*wp + 1) // 2**d
            a -= 1
            iter += 1
            print('m{} = {}'.format(1 + iter, d))
            if is_power_of_2(w*2**b - wp*3**a):
                #print_eq_lhs(a, w - 3**a * wp, iter)
                return
            else:
                print_eq_default(a, b, w, wp, iter)
        elif d == b:
            b = 0
            a -= 1
            m = d
            wp = (3*wp + 1) // 2**d
            r = w - 3**a * wp
            m += v2(r)
            r //= 2**v2(r)
            iter += 1
            print('m{} = {}'.format(1 + iter, m))
            
        else:
            m = b
            b = 0
            a -= 1
            #d not m so that we elimate the remaining pwoer of 2 and store it separately
            wp = (3*wp + 1) // 2**d
            r2 = d - m
            iter += 1
            print('m{} = {}'.format(1 + iter, m))
            print_eq_r2(a, w, wp, r2, iter)

        solve(a, b, w, wp, r2=r2, iter=iter)
    else:
        r = w - wp * 3**a * 2**r2 - 3**(a - 1)
        m = v2(r)
        iter += 1
        if r != 2**m:
            print(r)
            print(2**m)
            raise Exception("PORBLEMS")
        print('m{} = {}'.format(1 + iter, m))
        return

if __name__ == '__main__':
    a = 4
    b = 8
    w = 85
    wp = w

    solve(a, b, w, wp)
