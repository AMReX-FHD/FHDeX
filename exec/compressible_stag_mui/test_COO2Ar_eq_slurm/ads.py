import sys
import numpy as np
import math

file = sys.argv[1]
Ncellx = 4
Ncelly = 4

data = np.loadtxt(file)

x1 = data[0][1]
s1 = data[0][2]
x2 = data[0][3]
s2 = data[0][4]

x3 = data[0][5]
s3 = data[0][6]
x4 = data[0][7]
s4 = data[0][8]

meantheta1 = 4.036202e-01/2
vartheta1 = 5.349132e-06
meantheta2 = 4.349240e-01/2
vartheta2 = 5.461447e-06

print("E1[theta]   = %e\tstd= %e" % (x1,s1/math.sqrt(Ncellx*Ncelly)))
print("expected    = %e" % meantheta1)
print("Var1[theta] = %e\tstd= %e" % (x2,s2/math.sqrt(Ncellx*Ncelly)))
print("expected    = %e" % vartheta1)

print("E2[theta]   = %e\tstd= %e" % (x3,s3/math.sqrt(Ncellx*Ncelly)))
print("expected    = %e" % meantheta2)
print("Var2[theta] = %e\tstd= %e" % (x4,s4/math.sqrt(Ncellx*Ncelly)))
print("expected    = %e" % vartheta2)
