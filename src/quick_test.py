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

#Biggest n such that (r - 2^n) / 3 in R_(a-1)
def biggest_n_st(r, a):
    n = int(log2(r)) + 1
    while True:
        if in_Rn((r - 2**n) // 3, 2):
            return n
        n -= 1

if __name__ == '__main__':
    f = lambda x : 2*x + 4
    a = 3
    b = f(a)

    R = get_all_r(a, b)
    R = [r for r in R if r % 2 == 0]
    for r in R:
        print(r)
        print(inv_r(r, a - 1))
        print(biggest_n_st(r, a))
        print('')
