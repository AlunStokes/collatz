import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    u = 15
    W = [31, 62, 83, 166, 293, 347, 586, 671, 694, 1342, 2684, 19151, 38302, 2025797, 4051594]
    W = [w for w in W if w % 2 == 1]

    plt.plot(np.log10(W), 'b.-', label='known w')
    plt.plot([len(W)], [u], 'r.-', label='upper bound')
    plt.xlabel('index')
    plt.ylabel('log10(w)')
    plt.grid()
    plt.legend()
    plt.title('Growth rate of w < 10^{}'.format(u))
    plt.show()
