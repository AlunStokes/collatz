from utils import *

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def func(x, a, b):
    return a * np.exp(-b * x)

if __name__ == '__main__':
    W = [31, 62, 83, 166, 293, 347, 586, 671, 694, 1342, 2684, 19151, 38302, 2025797, 4051594]
    Wo = [w for w in W if w % 2 == 1]
    n = 3

    U = 10**4
    E = []
    while n < U:
        O = collatz(n)['values']
        e = max(O) / n
        E.append((n, e))
        n += 1
    plt.plot([e[0] for e in E], [e[1] for e in E], 'b')
    Ew = [e for e in E if e[0] in W]
    L = []
    for w in Wo:
        if w not in [e[0] for e in Ew]:
            continue
        i = 0
        while i < len(Ew):
            if Ew[i][0] == w:
                e1 = Ew[i]
            if Ew[i][0] == 2*w:
                e2 = Ew[i]
            i += 1
        b = np.log2(e2[1]/e1[1]) / (e2[0] - e1[0])
        a = e1[1] / 2**(b * e1[0])
        L.append((a, b))
        print((a,b))
    for l, w in zip(L, Wo):
        X = np.linspace(3, U, 1000)
        plt.plot(X, l[0] * 2**(l[1] * X), label='{}'.format(w))
    plt.plot([e[0] for e in Ew if e[0] % 2 == 1], [e[1] for e in Ew if e[0] % 2 == 1], 'r.')
    plt.plot([e[0] for e in Ew if e[0] % 2 == 0], [e[1] for e in Ew if e[0] % 2 == 0], 'g.')
    plt.ylim(bottom=-50)
    plt.legend()
    plt.show()
