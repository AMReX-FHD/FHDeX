#!/usr/bin/env python

import numpy as np
import pymbar

from pymbar import testsystems

#[x_n, u_kn, N_k, s_n] = testsystems.HarmonicOscillatorsTestCase().sample()

#mbar = MBAR(u_kn, N_k)


N = 2 #number of umbrellas
files = []
phi_initial = 0.05
phi_increment = 0.0025
for y in range (0,N):
    print('data'+str((y)*phi_increment+phi_initial))
    files.append('data'+str((y)*phi_increment+phi_initial))

x_kn = []
N_k = []

for filename in files:
    with open('data/{}'.format(filename)) as f:
        data = [float(x) for x in f.read().strip().split('\n')]
    x_kn.append(data[::500]) #skipping every 500 frames
    N_k.append(len(data[::500]))
#    x_kn.append(data[::100])
#    N_k.append(len(data[::100]))


N_k = np.array(N_k)
y_kn = np.zeros((len(N_k), N_k.max()))

for i in range(len(N_k)):
    for j in range(N_k[i]):
        y_kn[i, j] = x_kn[i][j]

x_n = np.array(sum(x_kn, []))

def U(k, x):
    #  lam = [0.003, 0.002, 0.001, 0.0, -0.001]
    lam = [] #input bias
    return lam[k] * x

u_kn = np.array([U(k, x_n) for k in range(len(N_k))])

mbar = pymbar.MBAR(u_kn, N_k, verbose = True)

x_n_min = x_n.min()

resolution = 1

x_n_adj = np.array([int((x - x_n_min) / resolution) for x in x_n])

f_k, df_k = mbar.computePMF(u_kn, x_n_adj, nbins = x_n_adj.max() + 1)

with open('FE.txt', 'w') as f:
    for k, (fe, err) in enumerate(zip(f_k, df_k)):
        print('{:4.2f} {:9.6f} {:8.6f}'.format(x_n_min + (k * resolution), fe, err), file = f)