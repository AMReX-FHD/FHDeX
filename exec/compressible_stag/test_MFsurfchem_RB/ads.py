import sys
import numpy as np
import math

file = sys.argv[1]
Ncell = 64

data = np.loadtxt(file)

x1 = data[0][1]
s1 = data[0][2]
x2 = data[0][3]
s2 = data[0][4]
x3 = data[0][5]
s3 = data[0][6]
x4 = data[0][7]
s4 = data[0][8]

meantheta1 = 1.350497e-01
vartheta1 = 1.297903e-06

print "E[theta1] = %e\tstd= %e" % (x1,s1/math.sqrt(Ncell))
print "E[theta2] = %e\tstd= %e" % (x2,s2/math.sqrt(Ncell))
print "expected  = %e" % meantheta1
print "Var[theta1] = %e\tstd= %e" % (x3,s3/math.sqrt(Ncell))
print "Var[theta2] = %e\tstd= %e" % (x4,s4/math.sqrt(Ncell))
print "expected    = %e" % vartheta1
