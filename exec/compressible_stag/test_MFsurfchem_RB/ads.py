import sys
import numpy as np
import math

file = sys.argv[1] 
Ncell = 64

data = np.loadtxt(file)

x1 = data[0][1]
v1 = data[0][2]
x2 = data[0][3]
v2 = data[0][4]
x3 = data[0][5]
v3 = data[0][6]
x4 = data[0][7]
v4 = data[0][8]

meantheta1 = 1.350497e-01
vartheta1 = 3.100865e-05

print "E[theta1] = %e (%e)" % (x1,v1/math.sqrt(Ncell))
print "E[theta2] = %e (%e)" % (x2,v2/math.sqrt(Ncell))
print "expected  = %e" % meantheta1
print "Var[theta1] = %e (%e)" % (x3,v3/math.sqrt(Ncell))
print "Var[theta2] = %e (%e)" % (x4,v4/math.sqrt(Ncell))
print "expected  = %e" % vartheta1
