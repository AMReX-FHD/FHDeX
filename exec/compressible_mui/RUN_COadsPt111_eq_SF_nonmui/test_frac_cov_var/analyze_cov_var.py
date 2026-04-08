import sys
import numpy as np

assert len(sys.argv) == 3, "len(sys.argv) == 3"

vv = np.loadtxt(sys.argv[1],usecols=(2),unpack=True)
Nsample = len(vv)

Ntot = 300*300
ra = 7.562496e+06
rd = 2.926784e+07
meancov = ra/(ra+rd)

exact_var = meancov*(1-meancov)/Ntot
mean_var = np.mean(vv)
std_var = np.std(vv,ddof=1)

print "exact_var = %e" % exact_var
print "mean_var  = %e" % mean_var
print "std_var   = %e" % std_var

T = 1e-12*1e4
a = (ra+rd)*T
ratio = 1-2*(a-1+np.exp(-a))/a**2

print "var/exact_var  = %f +/- %f" % (mean_var/exact_var,std_var/exact_var/np.sqrt(Nsample))
print "expected ratio = %f\n" % ratio

#####

vv2 = np.loadtxt(sys.argv[2],usecols=(2),unpack=True)

mean_var2 = np.mean(vv2)
std_var2 = np.std(vv2,ddof=1)

print "exact_var = %e" % exact_var
print "mean_var2 = %e" % mean_var2
print "std_var2  = %e" % std_var2

print "var2/exact_var = %f +/- %f" % (mean_var2/exact_var,std_var2/exact_var/np.sqrt(Nsample))
print "expected ratio = %f" % 1
