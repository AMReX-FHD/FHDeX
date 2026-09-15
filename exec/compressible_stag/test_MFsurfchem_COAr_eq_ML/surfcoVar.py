import sys
import numpy as np
import math

file = sys.argv[1]
Ncellx = 16
Ncelly = 16

data = np.loadtxt(file)

x1 = data[0][1]
s1 = data[0][2]
x2 = data[0][3]
s2 = data[0][4]
x3 = data[0][5]
s3 = data[0][6]
x4 = data[0][7]
s4 = data[0][8]
x5 = data[0][9]
s5 = data[0][10]
x6 = data[0][11]
s6 = data[0][12]

varrhoYkT = 0
varthetarhoYk = 0
varthetavx = 0
varthetavy = 0
varthetavz = 0
varthetaT = 0

print("Var[rhoYk-T]   = %e\tstd= %e" % (x1,s1/math.sqrt(Ncellx*Ncelly)))
print("expected       = %e" % varrhoYkT)
print("Var[theta-rhoYk] = %e\tstd= %e" % (x2,s2/math.sqrt(Ncellx*Ncelly)))
print("expected       = %e" % varthetarhoYk)
print("Var[theta-vx]  = %e\tstd= %e" % (x3,s3/math.sqrt(Ncellx*Ncelly)))
print("expected       = %e" % varthetavx)
print("Var[theta-vy]  = %e\tstd= %e" % (x4,s4/math.sqrt(Ncellx*Ncelly)))
print("expected       = %e" % varthetavy)
print("Var[theta-vz]  = %e\tstd= %e" % (x5,s5/math.sqrt(Ncellx*Ncelly)))
print("expected       = %e" % varthetavz)
print("Var[theta-T]   = %e\tstd= %e" % (x6,s6/math.sqrt(Ncellx*Ncelly)))
print("expected       = %e" % varthetaT)
