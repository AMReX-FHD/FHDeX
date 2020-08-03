import numpy as np

# number density of each species

n1 = 1.817055e+19
n2 = 3.604090e+18
n3 = 1.820603e+18
n4 = 8.678933e+17

# adsorption/desorption parameters

cads1 = 2.015564e-09
rdes1 = 3.296152e+11

cads2 = 8.976568e-10
rdes2 = 2.911713e+10

rdes3 = 2.206933e+10
cads3 = 6.379993e-10

rdes4 = 7.263829e+09
cads4 = 4.404998e-10

# equilibrium coverage

rads1 = cads1*n1
rads2 = cads2*n2
rads3 = cads3*n3
rads4 = cads4*n4

A = np.array([[rads1+rdes1,rads1,rads1,rads1],[rads2,rads2+rdes2,rads2,rads2],[rads3,rads3,rads3+rdes3,rads3],[rads4,rads4,rads4,rads4+rdes4]])
b = np.array([rads1,rads2,rads3,rads4])
x = np.linalg.solve(A,b)

# output

nsite = 40*40
vac = 1-sum(x)
p = np.append(vac,x)

mean = nsite*p
var = nsite*np.multiply(p,1-p)

dens = np.array([n1,n2,n3,n4])
cads = np.array([cads1,cads2,cads3,cads4])
rdes = np.array([rdes1,rdes2,rdes3,rdes4])
print "number density: ", dens
print "cads: ", cads
print "rdes: ", rdes
print

print "number of sites: %d" % nsite
print "p    (vac,spec1,spec2,spec3,spec4): ", p
print "mean (vac,spec1,spec2,spec3,spec4): ", mean
print "var  (vac,spec1,spec2,spec3,spec4): ", var
print

