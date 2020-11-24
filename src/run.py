import utils

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    '''D = {}

    N = np.arange(2, 5000, 1)
    for n in N:
        coord = utils.get_coord(n)
        if coord in D:
            D[coord].append(n)
        else:
            D[coord] = [n]

    for d in D:
        print('{}: {}'.format(d,D[d]))

    print(len(list(D.keys())))'''

    '''seq = utils.collatz(201)['values']
    m = max(seq)
    l = max([len(str(s)) for s in seq])
    print(seq)

    V = []
    for n in seq:
        V.append(utils.dict_to_vec(utils.get_factor_vector(n, pad_to = m)))

    Vn = np.array(V)

    #print('{}{}: {}'.format(seq[0], ' ' * (l - len(str(seq[0]))), V[0]))
    i = 1
    while i < len(seq):
        if seq[i] > seq[i - 1]:
            #print('{}{}: {}'.format(seq[i], ' ' * (l - len(str(seq[i]))), V[i]))
            print(Vn[i].nonzero())
        i += 1'''

    '''N = np.arange(3000, 3200, 1)
    print(N)
    utils.plot_paths(N)

    for n in N:
        print(utils.collatz(n))'''

    '''for i in range(100, 200):
        #N = np.arange(i, i + 100, 2)
        N = np.random.randint(100, 500, 100)
        utils.plot_ordered_iterations(N, 'random-100', save=True)'''

    #Check number of consecutive (3x + 1) / 2 setps
    '''N = np.arange(200000001, 201000001, 2)
    L = []
    for n in N:
        j = 0
        while n % 2 == 1:
            n = (3 * n + 1) // 2
            j += 1
        L.append(j)

    R = [(n, l) for n, l in zip (N, L)]
    R = sorted(R, key = lambda x: x[1])[::-1]
    print(R[:10])'''

    #Check consecutive prime factorizations
    '''N = np.arange(2, 100, 1)
    for n in N:
        v = utils.dict_to_vec(utils.get_factor_vector(n))
        s = np.sum(v)
        print(s, v)'''

    '''N = np.arange(2, 100)
    for n in N:
        r = utils.mod_inv(2, n)
        if r > 0:
            m = (r * (n - 1)) % n
            if m % 2 == 0:
                print(n)

    exit()'''

    '''N = np.arange(3, 1000000, 2)
    #num, degree
    m = (0, 0)
    for n in N:
        r = utils.divergence_degree(n)
        if r > m[1]:
            m = (n, r)
    print(m)

    exit()'''


    #Check working fininte fields for collatz
    #P = utils.get_primes_up_to(500)
    P = np.arange(3, 160)
    L = [1]
    W = []
    #print(N)
    for p in P:
        for l in L:
            t, b = utils.check_finite_collatz(p, l)
            if t:
                W.append((p, l))
            else:
                pass
                #print("{}: {}".format((p, l), b))
    print("Works, known: ", W)
    V = []
    #print(N)
    for p in P:
        for l in L:
            t, v = utils.check_finite_conjectures(p, l)
            if t:
                V.append((p, l))
            else:
                pass
                #print("{}: {}".format((p, l), b))
    print("Wokrs, conjecture: ", V)
    print("")

    X = [x for x in V if x not in W]
    Y = [y for y in W if y not in V]
    print("Works, missed by conj.: ", Y)
    print("Doesnt work, picked by conj.: ", X)

    for x in X:
        print(x, ': ', utils.check_finite_collatz(x[0], 1))


    #utils.generate_circle_diagram(24)


    '''Z = []
    N = np.arange(0, 63)
    for n in N:
        Z.append(int(-1/3 * (-1 + (-1/2)**(-n))))
    print(Z)
    for n, z in zip(N, Z):
        L = []
        M = np.arange(2**(n + 1) - z, 10 * 2**(n + 1) - (2**(n + 1) - z), 2**(n + 1))
        #print(M)
        for m in M:
            L.append(utils.phi(m))
            #print('{}: {}'.format(n, utils.phi(n)))
        print(L)
    #plt.plot(L)
    #plt.show()'''

    #utils.even_odd_coord_graph_grid(10, 60)
