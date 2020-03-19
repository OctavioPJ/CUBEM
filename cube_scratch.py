import re
import numpy as np
from glob import glob
from Cube import Charm
import matplotlib.pyplot as plt

if __name__ == '__main__':
    alf = [(np.loadtxt(file), file) for file in glob('evolution/*.txt')]
    exact = []
    numer = []
    N = []
    charming = Charm(n0=100E6)
    for found, file in alf:
        charming.ExactSolution(found[-1, 0])
        exact.append(charming.Solution)
        numer.append(found[-1, 1])
        N.append(int(re.findall('evol_([0-9]+)(?:_True)?\.txt', file)[0]))

    plt.loglog(N, abs(np.array(numer)-np.array(exact))/np.array(exact)*100, 'o')
