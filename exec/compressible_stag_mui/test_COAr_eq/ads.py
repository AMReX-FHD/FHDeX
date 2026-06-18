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

meantheta1 = 2.053333e-01
vartheta1 = 1.813017e-06

print("E[theta]   = %e\tstd= %e" % (x1,s1/math.sqrt(Ncellx*Ncelly)))
print("expected   = %e" % meantheta1)
print("Var[theta] = %e\tstd= %e" % (x2,s2/math.sqrt(Ncellx*Ncelly)))
print("expected   = %e" % vartheta1)
