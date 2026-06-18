import numpy as np

# number density of each species

n1 = 1.817055e+19
n2 = 3.604090e+18
n3 = 1.820603e+18
n4 = 8.678933e+17

# collision rates onto a site (this is proportional to site area)

rcol1 = 3.662392e+11
rcol2 = 3.235236e+10
rcol3 = 1.161544e+10
rcol4 = 3.823068e+09

# adsorption/desorption parameters

cads1 = 4.031128e-09
rdes1 = 6.592305e+11

cads2 = 2.692970e-09
rdes2 = 1.844085e+11

cads3 = 1.594998e-09
rdes3 = 1.161544e+10

cads4 = 2.642999e-09
rdes4 = 2.293841e+09

# conversion probabilities

rads1 = cads1*n1
rads2 = cads2*n2
rads3 = cads3*n3
rads4 = cads4*n4

print "conversion probabilities: ", [rads1/rcol1,rads1/rcol1,rads1/rcol1,rads1/rcol1]

# single-component equilibirium coverage

print "single-component equilibrium coverage: ", [rads1/(rads1+rdes1),rads2/(rads2+rdes2),rads3/(rads3+rdes3),rads4/(rads4+rdes4)]
print

# equilibrium coverage

A = np.array([[rads1+rdes1,rads1,rads1,rads1],[rads2,rads2+rdes2,rads2,rads2],[rads3,rads3,rads3+rdes3,rads3],[rads4,rads4,rads4,rads4+rdes4]])
b = np.array([rads1,rads2,rads3,rads4])
x = np.linalg.solve(A,b)

# output

nsite = 80*80
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

