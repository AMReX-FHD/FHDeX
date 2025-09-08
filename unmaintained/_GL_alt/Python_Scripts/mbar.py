#!/usr/bin/env python

import numpy as np
import pymbar

from pymbar import testsystems




N =999  #number of umbrellas
files = [] # inmitialize array of name of data files to be read
phi_initial = 0.001
phi_increment = 0.001
Umbrella=100


phi_0_arr=[0]*N
for index1 in range (0 ,N):
        phi_0_arr[index1]=phi_initial + index1*phi_increment
for y in range (0,N):
       place_holder =str(format(((y)*phi_increment+phi_initial), '.3f'))
       files.append('data'+place_holder) #add names of files to be read to array. Assuing form is datax.xx
       print('data'+place_holder) #print out for verifying correct data file names have been chosen
x_kn = [] #1D harmonic oscillator positions
N_k = [] # Number of samples drawn from each state


for filename in files: #for loop in entries of files array
    with open('data/{}'.format(filename)) as f: # open a particular data file in loop
         data = [float(x) for x in f.read().strip().split('\n')]
    x_kn.append(data[::500]) #skipping every 500 frames. Set to get one position per umbrella state
    N_k.append(len(data[::500])) #append number of position from umbrella. ( set to 1 6/13/19)
#    x_kn.append(data[::100])
#    N_k.append(len(data[::100]))


N_k = np.array(N_k)# Array of where each entry is the number of positions in each umbrella
y_kn = np.zeros((len(N_k), N_k.max()))


for i in range(len(N_k)):
    for j in range(N_k[i]):
        y_kn[i, j] = x_kn[i][j]

x_n = np.array(sum(x_kn, []))


def U(k, x):
     return 0.5 *Umbrella*((x - phi_0_arr[k]) ** 2)

u_kn = np.array([U(k, x_n) for k in range(len(N_k))])


mbar = pymbar.MBAR(u_kn, N_k, verbose = True)


#x_n_min = x_n.min()
#resolution = 1
#x_n_adj = np.array([int((x - x_n_min) / resolution) for x in x_n])
#f_k, df_k = mbar.computePMF(u_kn[4,:], x_n_adj, nbins = x_n_adj.max() + 1)
#with open('FE.txt', 'w') as f:
#    for k, (fe, err) in enumerate(zip(f_k, df_k)):
#        print('{:4.2f} {:9.6f} {:8.6f}'.format(x_n_min + (k * resolution), fe, err), file = f)

results  = mbar.getFreeEnergyDifferences(compute_uncertainty=True, uncertainty_method=None, warning_cutoff=1e-10, return_theta=True)
Delta_f=results[0] # Exract free energy differences from results tuple
Final_Dimless_Energy=Delta_f[998,:] #This extracts the free energy difference between the the first state and the rest of the states only

f= open("FE.txt","w+")
for i in range(0,N-1):
     f.write("%+3.16E %+3.16E   \n"% (Final_Dimless_Energy[i],x_n[i] ))
f.close()
