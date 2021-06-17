import os

import numpy as np
import matplotlib.pyplot as plt

def get_as_congruence_class(n, w):
    r = n % w
    q = (n - r) // w
    return q, r

if __name__ == '__main__':
    files = os.listdir('./neg_inv')

    D = {}

    for file in files:
        if '.sh' in file:
            continue
        ind_a = file.index('a')
        ind_dot = file.index('.')
        k = int(file[1:ind_a])
        a = int(file[ind_a + 1:ind_dot])

        with open('./neg_inv/{}'.format(file), 'r') as f:
            r = f.read()
        lines = r.split('\n')[1:-2]

        W = []
        for line in lines:
            w = int(line.split('\t')[-1])
            W.append(w)
        W = sorted(W)
        #print('For k = {} up to 10^{}:'.format(k, a))
        #print(W)
        #print(len(W))
        #print(W[-1] / 10**a)

        if k in D:
            D[k][a] = W
        else:
            D[k] = {}
            D[k][a] = W

    longest_w = -1

    K = sorted(list(D.keys()))

    for k in K:
        A = sorted(list(D[k].keys()))

        '''print('k = {}'.format(k))
        for a in A:
            W = D[k][a]
            if len(W) > longest_w:
                longest_w = len(W)
            print('  Up to 10^{}, we have {} such w'.format(a, len(W)))
            print('  Largest element {} ({})'.format(W[-1], W[-1] / 10**a))
            print('  Proportion of candidate w is {}'.format(len(W) / (10**a * (1 - 1/k))))

            #plt.title('k={} 10^{}'.format(k, a))
            #plt.plot(np.log(W))
            #plt.plot([len(W)], [np.log(10**a)], 'r.')
            #plt.show()

            print()'''

        #For LaTeX
        a = A[-1]
        W = D[k][a]
        while W[0] < 5:
            del W[0]
        s = '${}$ & ${}$ & ${} \\approx 10^{{ {} }}$ & $\\leq 10^{{ {} }}$ \\\\'
        max_order = int(np.log10(W[-1]))
        prop = len(W) / ((10**a - k) * (1 - 1/k))
        prop_order = int(np.log10(prop))
        s = s.format(k, len(W), W[-1], max_order, prop_order)
        print(s)
