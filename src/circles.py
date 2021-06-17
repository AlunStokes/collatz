import time

from math import floor
from multiprocessing.pool import Pool
from functools import partial
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

from utils import *
from test_r import *

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

def collatz(x):
    code = ''
    n = 0

    l = [x]
    while x != 1:
        n += 1
        if x % 2 == 1:
            x = (3*x + 1)
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

#gcd by euclidean algorithm
def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

#modular inverse
def mod_inv(a, m):
    g, x, y = egcd(a, m)
    if g != 1:
        return None
    else:
        return x % m

def get_as_congruence_class(n, w):
    r = n % w
    q = (n - r) // w
    return q, r

def generate_circle_graph(w):
    if w % 3 == 0:
        return -1
    d = w % 3
    r1 = (d*w - 1) // 3
    if d == 1:
        r2 = (2*w - 1) // 3
    else:
        r2 = (w - 1) // 3
    O = collatz(w)
    parity = list([int(x) for x in list(O['code'])])
    O = O['values']
    Ow = [o % w for o in O]
    P = {}
    angle_delta = 2*np.pi / w
    i = 0
    while i < w:
        angle = (np.pi / 2)  - i * angle_delta
        P[i] = (1 * np.cos(angle), 1 * np.sin(angle))
        i += 1

    #xy axis
    #plt.plot([0, 0], [-1.2, 1.2], 'k')
    #plt.plot([-1.2, 1.2], [0, 0], 'k')

    for p in P:
        if p == r1:
            plt.plot([P[p][0]], [P[p][1]], 'r.', markersize=6)
        elif p == r2:
            plt.plot([P[p][0]], [P[p][1]], 'g.', markersize=6)
        else:
            plt.plot([P[p][0]], [P[p][1]], 'b.', markersize=5)

    q_max = max([get_as_congruence_class(o, w)[0] for o in O])

    i = 1
    while i < len(Ow):
        o = O[i]
        q,r = get_as_congruence_class(o, w)
        #plt.plot([P[Ow[i-1]][0], P[Ow[i]][0]], [P[Ow[i-1]][1], P[Ow[i]][1]], 'r')
        x1 = P[Ow[i - 1]][0]
        y1 = P[Ow[i - 1]][1]
        x2 = P[Ow[i]][0]
        y2 = P[Ow[i]][1]
        if parity[i - 1] == 0:
            plt.arrow(x1, y1, x2 - x1, y2 - y1, head_width=0.02, linewidth=0.4, head_length=0.05, length_includes_head=True, color='b')
        else:
            plt.arrow(x1, y1, x2 - x1, y2 - y1, head_width=0.02, linewidth=0.4, head_length=0.05, length_includes_head=True, color='r')
        i += 1
    #plt.show()

def generate_cylinder_graph(w):
    if w % 3 == 0:
        return -1
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    O = collatz(w)['values']
    V = [get_as_congruence_class(o, w) for o in O]
    lower_q = min([v[0] for v in V])
    upper_q = max([v[0] for v in V])
    P = {}
    angle_delta = 2*np.pi / w
    i = 0
    while i < w:
        angle = (np.pi / 2)  - i * angle_delta
        P[i] = (1 * np.cos(angle), 1 * np.sin(angle))
        i += 1

    z = lower_q
    while z <= upper_q:
        ax.scatter3D([P[i][0] for i in P], [P[i][1] for i in P], [z for _ in range(len(P))])
        z += 1

    i = 1
    while i < len(O):
        #plt.plot([P[Ow[i-1]][0], P[Ow[i]][0]], [P[Ow[i-1]][1], P[Ow[i]][1]], 'r')
        x1 = P[V[i - 1][1]][0]
        y1 = P[V[i - 1][1]][1]
        z1 = V[i - 1][0]
        x2 = P[V[i][1]][0]
        y2 = P[V[i][1]][1]
        z2 = V[i][0]
        ax.plot3D([x1, x2], [y1, y2], [z1, z2])
        i += 1
    #plt.show()

def generate_multi_circle_graph(w):
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    O = collatz(w)['values']
    V = [get_as_congruence_class(o, w) for o in O]

    max_q = max([v[0] for v in V])
    theta = 2 * np.pi / w
    phi = 2 * np.pi / (max_q + 1)

    P = {}
    for i in range(0, w):
        angle = (-np.pi / 2)  - i * theta
        P[i] = (1 * np.cos(angle), 1 * np.sin(angle) + 1)
        i += 1

    Ps = {}
    for i in range(0, max_q + 1):
        C = {}
        for j in range(0, w):
            C[j] = (P[j][0], np.cos(i * phi) * 0 - np.sin(i * phi) * P[j][1], np.sin(i * phi) * 0 + np.cos(i * phi) * P[j][1])
        Ps[i] = C

    for i in range(0, max_q + 1):
        X = []
        Y = []
        Z = []
        for j in range(0, w):
            X.append(Ps[i][j][0])
            Y.append(Ps[i][j][1])
            Z.append(Ps[i][j][2])
        #ax.scatter3D(X, Y, Z)

    ax.scatter3D([0], [0], [0], 'k.')

    i = 1
    while i < len(O):
        #plt.plot([P[Ow[i-1]][0], P[Ow[i]][0]], [P[Ow[i-1]][1], P[Ow[i]][1]], 'r')
        q1, r1 = V[i - 1]
        q2, r2 = V[i]
        x1 = Ps[q1][r1][0]
        x2 = Ps[q2][r2][0]
        y1 = Ps[q1][r1][1]
        y2 = Ps[q2][r2][1]
        z1 = Ps[q1][r1][2]
        z2 = Ps[q2][r2][2]

        ax.plot3D([x1, x2], [y1, y2], [z1, z2])
        #ax.scatter3D([x1], [y1], [z1], '.')
        i += 1

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim(-1.5, 1.5)
    plt.show()

def generate_even_circle(w):
    P = {}
    angle_delta = 2*np.pi / w
    i = 0
    while i < w:
        angle = (np.pi / 2)  - i * angle_delta
        P[i] = (1 * np.cos(angle), 1 * np.sin(angle))
        i += 1

    #xy axis
    #plt.plot([0, 0], [-1.2, 1.2], 'k')
    #plt.plot([-1.2, 1.2], [0, 0], 'k')

    for p in P:
        if w % 3 != 0:
            if p == w - mod_inv(3, w):
                plt.plot([P[p][0]], [P[p][1]], 'r.', markersize=5)
            else:
                plt.plot([P[p][0]], [P[p][1]], 'b.', markersize=5)
        else:
            plt.plot([P[p][0]], [P[p][1]], 'b.', markersize=5)

    inv2 = mod_inv(2, w)

    for v in range(w):
        x2 = P[v][0]
        y2 = P[v][1]
        v_even = inv2 * v % w
        #v_odd = (3 * v + 1) % w
        x_even = P[v_even][0]
        y_even = P[v_even][1]
        #x_odd = P[v_odd][0]
        #y_odd = P[v_odd][1]
        plt.arrow(x2, y2, x_even - x2, y_even - y2, head_width=0.02, linewidth=0.4, head_length=0.05, length_includes_head=True, color='b')
        #plt.arrow(x2, y2, x_odd - x2, y_odd - y2, head_width=0.02, linewidth=0.4, head_length=0.05, length_includes_head=True, color='r')
        i += 1

    #plt.show()

def generate_odd_circle(w):
    P = {}
    angle_delta = 2*np.pi / w
    i = 0
    while i < w:
        angle = (np.pi / 2)  - i * angle_delta
        P[i] = (1 * np.cos(angle), 1 * np.sin(angle))
        i += 1

    #xy axis
    #plt.plot([0, 0], [-1.2, 1.2], 'k')
    #plt.plot([-1.2, 1.2], [0, 0], 'k')

    for p in P:
        if w % 3 != 0:
            if p == w - mod_inv(3, w):
                plt.plot([P[p][0]], [P[p][1]], 'r.', markersize=5)
            else:
                plt.plot([P[p][0]], [P[p][1]], 'b.', markersize=5)
        else:
            plt.plot([P[p][0]], [P[p][1]], 'b.', markersize=5)

    inv2 = mod_inv(2, w)

    for v in range(w):
        x2 = P[v][0]
        y2 = P[v][1]
        #v_even = inv2 * v % w
        v_odd = (3 * v + 1) % w
        #x_even = P[v_even][0]
        #y_even = P[v_even][1]
        x_odd = P[v_odd][0]
        y_odd = P[v_odd][1]
        #plt.arrow(x2, y2, x_even - x2, y_even - y2, head_width=0.02, linewidth=0.4, head_length=0.05, length_includes_head=True, color='b')
        plt.arrow(x2, y2, x_odd - x2, y_odd - y2, head_width=0.02, linewidth=0.4, head_length=0.05, length_includes_head=True, color='r')
        i += 1

    #plt.show()

def generate_neg_inv_circle(w, include_ops=False):
    E = [1]
    i = 2
    while i < w:
        if egcd(i, w)[0] == 1:
            E.append(i)
        i += 1
    V = [(mod_inv(e, w)) % w for e in E]
    V_neg = [(w - mod_inv(e, w)) % w for e in E]

    P = {}
    angle_delta = 2*np.pi / w
    i = 0
    while i < w:
        angle = (np.pi / 2)  - i * angle_delta
        P[i] = (1 * np.cos(angle), 1 * np.sin(angle))
        i += 1

    #xy axis
    #plt.plot([0, 0], [-1.2, 1.2], 'k')
    #plt.plot([-1.2, 1.2], [0, 0], 'k')

    for p in P:
        if w % 3 != 0:
            if p == w - mod_inv(3, w):
                plt.plot([P[p][0]], [P[p][1]], 'r.', markersize=5)
            else:
                plt.plot([P[p][0]], [P[p][1]], 'b.', markersize=5)
        else:
            plt.plot([P[p][0]], [P[p][1]], 'b.', markersize=5)


    for e, v in zip(E, V_neg):
        x1 = P[e][0]
        y1 = P[e][1]
        x2 = P[v][0]
        y2 = P[v][1]
        plt.arrow(x1, y1, x2 - x1, y2 - y1, head_width=0.02, head_length=0.05, length_includes_head=True, color='k')
        if include_ops:
            inv2 = mod_inv(2, w)
            v_even = inv2 * v % w
            v_odd = (3 * v + 1) % w
            x_even = P[v_even][0]
            y_even = P[v_even][1]
            x_odd = P[v_odd][0]
            y_odd = P[v_odd][1]
            plt.arrow(x2, y2, x_even - x2, y_even - y2, head_width=0.02, linewidth=0.4, head_length=0.05, length_includes_head=True, color='b')
            plt.arrow(x2, y2, x_odd - x2, y_odd - y2, head_width=0.02, linewidth=0.4, head_length=0.05, length_includes_head=True, color='r')
        i += 1

    '''for e, v in zip(E, V):
        x1 = P[e][0]
        y1 = P[e][1]
        x2 = P[v][0]
        y2 = P[v][1]
        #v_even =
        plt.arrow(x1, y1, x2 - x1, y2 - y1, head_width=0.02, head_length=0.05, length_includes_head=True, color='g')
        i += 1'''
    #plt.show()

def plot_connections(C, w, r=-1, highlight_neg3inv = True):
    P = {}
    angle_delta = 2*np.pi / w
    i = 0
    while i < w:
        angle = (np.pi / 2)  - i * angle_delta
        P[i] = (1 * np.cos(angle), 1 * np.sin(angle))
        i += 1

    for p in P:
        markersize = 5
        if p == r:
            markersize = 10
        if w % 3 != 0 and highlight_neg3inv:
            if p == w - mod_inv(3, w):
                plt.plot([P[p][0]], [P[p][1]], 'r.', markersize=markersize)
            else:
                plt.plot([P[p][0]], [P[p][1]], 'b.', markersize=markersize)
        else:
            plt.plot([P[p][0]], [P[p][1]], 'b.', markersize=markersize)

    for c in C:
        x1 = P[c[0][0]][0]
        y1 = P[c[0][0]][1]
        x2 = P[c[0][1]][0]
        y2 = P[c[0][1]][1]
        if c[1] == 0:
            plt.arrow(x1, y1, x2 - x1, y2 - y1, head_width=0.02, linewidth=0.4, head_length=0.05, length_includes_head=True, color='b')
        else:
            plt.arrow(x1, y1, x2 - x1, y2 - y1, head_width=0.02, linewidth=0.4, head_length=0.05, length_includes_head=True, color='r')

def generate_circle_paths(w):
    E = [1]
    i = 2
    while i < w:
        if egcd(i, w)[0] == 1:
            E.append(i)
        i += 1
    V = [(mod_inv(e, w)) % w for e in E]
    V_neg = [(w - mod_inv(e, w)) % w for e in E]

    P = {}
    angle_delta = 2*np.pi / w
    i = 0
    while i < w:
        angle = (np.pi / 2)  - i * angle_delta
        P[i] = (1 * np.cos(angle), 1 * np.sin(angle))
        i += 1

    for r in range(1, w):
        Cs = []
        for q in range(1, 10):
            C = []
            wq = q * w + r
            O = collatz(wq)
            Par = list(O['code'])
            O = O['values']
            Ow = [o % w for o in O]

            i = 1
            while i < len(Ow):
                C.append(((Ow[i - 1], Ow[i]), int(Par[i - 1])))
                i += 1

            plot_connections(C, w, r=r)
            plt.savefig('./images/circles/test/{}/{}_r{}_q{}.png'.format(w, w, r, q))
            plt.clf()
            Cs.append(C)
        Cm = set(Cs[0])
        for i in range(1, len(Cs)):
            Cm = Cm.intersection(set(Cs[i]))

        plot_connections(Cm, w, r=r)
        plt.savefig('./images/circles/test/{}/{}_r{}.png'.format(w, w, r))
        plt.clf()



    '''for e, v in zip(E, V):
        x1 = P[e][0]
        y1 = P[e][1]
        x2 = P[v][0]
        y2 = P[v][1]
        #v_even =
        plt.arrow(x1, y1, x2 - x1, y2 - y1, head_width=0.02, head_length=0.05, length_includes_head=True, color='g')
        i += 1'''
    #plt.show()

def random_walk_on_circle(w, iter, num_generate=10):
    count_with = 0
    avg_prop_odd_with = 0
    count_without = 0
    avg_prop_odd_without = 0
    i = 0
    while i < num_generate:
        steps = [1 * (np.random.rand() < 0.5)]
        while len(steps) < iter:
            if steps[-1] == 1:
                steps.append(0)
            else:
                steps.append(1 * (np.random.rand() < 0.5))

        r = w - mod_inv(3, w)
        inv2 = mod_inv(2, w)
        C = []
        found = False
        v0 = 1
        for p in steps:
            if p == 0:
                v1 = v0 * inv2 % w
                C.append(((v0, v1), 0))
                v0 = v1
            else:
                v1 = (3 * v0 + 1) % w
                C.append(((v0, v1), 1))
                v0 = v1
            if v1 == r:
                count_with += 1
                avg_prop_odd_with += sum([c[1] for c in C]) / len(C)
                found = True
                break
        if not found:
            count_without += 1
            avg_prop_odd_without += sum([c[1] for c in C]) / len(C)


        #plot_connections(C, w)
        #plt.savefig('./images/circles/random_walk/{}_iter{}_({}).png'.format(w, iter, i))
        #plt.clf()
        i += 1
    avg_prop_odd_with /= count_with
    avg_prop_odd_without /= count_without
    print('{}/{} ~= {} contained -3^(-1)'.format(count_with, num_generate, count_with/num_generate))
    print('These had an average proportion of {} odd steps'.format(avg_prop_odd_with))
    print('This compares to {} for those without -3^1'.format(avg_prop_odd_without))

def generate_semicircle(center_x, center_y, radius, stepsize=0.1, orientation='right'):
    """
    generates coordinates for a semicircle, centered at center_x, center_y
    """

    x = np.arange(center_x, center_x+radius+stepsize, stepsize)
    x[x > radius] = radius
    try:
        y = np.sqrt(radius**2 - x**2)
    except:
        print(radius)
        print(x)
        exit()


    # since each x value has two corresponding y-values, duplicate x-axis.
    # [::-1] is required to have the correct order of elements for plt.plot.
    if orientation == 'left':
        x = -np.concatenate([x,x[::-1]])
    else:
        x = np.concatenate([x,x[::-1]])

    # concatenate y and flipped y.
    y = np.concatenate([y,-y[::-1]])

    return x, y + center_y

def linear_diagram(w, orbit=False):
    r = w - mod_inv(3, w)

    if orbit:
        O = collatz(w)
        parity = [int(n) for n in list(O['code'])]
        O = O['values']
        Ow = [o % w for o in O]
        i = 1
        while i < len(Ow):
            n = Ow[i - 1]
            v = Ow[i]
            if v == n:
                i += 1
                continue

            y = (n + v) / 2
            rad = abs(n - v) / 2

            if parity[i - 1] == 0:
                plt.plot(*generate_semicircle(0, y, rad, 0.01), 'b')
            else:
                plt.plot(*generate_semicircle(0, y, rad, 0.01, orientation='left'), 'r')
            i += 1
    else:
        V = {}
        inv2 = mod_inv(2,w)
        for n in range(w):
            v_even = inv2*n % w
            v_odd = (3*n + 1) % w
            V[n] = [v_even, v_odd]

        for n in range(w):
            v = V[n][0]
            if v == n:
                continue
            y = (n + v) / 2
            rad = abs(n - v) / 2
            plt.plot(*generate_semicircle(0, y, rad, 0.01), 'b')

            v = V[n][1]
            y = (n + v) / 2
            rad = abs(n - v) / 2
            plt.plot(*generate_semicircle(0, y, rad, 0.01, orientation='left'), 'r')

    plt.plot(w * [0], list(range(w)), 'k.', markersize=5)

    plt.plot([0], [r], 'r.', markersize=7)

    plt.xlim((-w + 1) // 2, (w - 1) // 2)
    #plt.show()


if __name__ == '__main__':
    '''for w in range(5, 10**3):
        if w % 2 == 0 or w % 3 == 0:
            continue
        linear_diagram(w, orbit=True)
        #plt.savefig('./images/linear_residues/all/{}.png'.format(w))
        plt.savefig('./images/linear_residues/orbit/{}.png'.format(w))
        plt.clf()'''

    '''for w in range(11, 45):
        if w % 2 == 0 or w % 3 == 0:
            continue
        if not os.path.exists('./images/circles/test/{}'.format(w)):
            os.makedirs('./images/circles/test/{}'.format(w))
        generate_even_circle(w)
        generate_odd_circle(w)
        plt.savefig('./images/circles/test/{}/{}.png'.format(w, w))
        plt.clf()
        generate_circle_paths(w)'''


    #8n^2
    '''w = 2
    n = 1
    while n < 12:
        r = 2**n * w
        generate_neg_inv_circle(r)
        plt.savefig('./images/circles/n^2/{}.png'.format(r))
        plt.clf()
        n += 1
    exit()'''

    #Check if pseudo even and odd operations maintain residue across all multiples of w for each coprime residue
    '''w = 15
    inv2 = mod_inv(2,w)
    E = [k for k in range(1, w) if egcd(k, w)[0] == 1]
    for e in E:
        print('e = {}'.format(e))
        for k in range(0, 5):
            v = k * w + e
            v_even = v * inv2
            q_even, r_even = get_as_congruence_class(v_even, w)
            v_odd = 3*v + 1
            q_odd, r_odd = get_as_congruence_class(v_odd, w)
            print('{}*w + {} -> {}*w + {}  |  {}*w + {}'.format(k, e, q_even, r_even, q_odd, r_odd))
        print()
    exit()'''

    '''w = 7
    while w < 10**3:
        if w % 2 == 0:
            w += 3
            continue
        generate_neg_inv_circle(w, include_ops=True)
        plt.savefig('./images/circles/mixed/{}_all.png'.format(w))
        #plt.savefig('./images/circles/all/{}.png'.format(w))
        plt.clf()

        generate_neg_inv_circle(w, include_ops=False)
        plt.savefig('./images/circles/mixed/{}_neg_inv.png'.format(w))
        #plt.savefig('./images/circles/neg_inv/{}.png'.format(w))
        plt.clf()

        generate_even_circle(w)
        plt.savefig('./images/circles/mixed/{}_even.png'.format(w))
        #plt.savefig('./images/circles/even/{}.png'.format(w))
        plt.clf()

        generate_odd_circle(w)
        plt.savefig('./images/circles/mixed/{}_odd.png'.format(w))
        #plt.savefig('./images/circles/odd/{}.png'.format(w))
        plt.clf()
        w += 3'''

    w = 7
    while w < 10**3:
        if w % 2 == 0:
            w += 3
            continue
        generate_circle_graph(w)
        plt.savefig('./images/circles/orbit/{}.png'.format(w))
        #plt.savefig('./images/circles/all/{}.png'.format(w))
        plt.clf()
        w += 3
