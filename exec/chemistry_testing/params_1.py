import math

k1 = 0.3
k2 = 0.5

n1 = 10.
n2 = 10.

A = n1+2*n2

# k1*(A-2*n2)^2 = k2*n2
# --> 4*k1*n2^2 - (4*k1*A+k2)*n2 + k1*A^2 = 0

a = 4*k1
b = -(4*k1*A+k2)
c = k1*A**2

dV = 1.e3

sol1 = (-b+math.sqrt(b**2-4.*a*c))/2./a
sol2 = (-b-math.sqrt(b**2-4.*a*c))/2./a


n1_bar = (A-2.*sol2)

Var_n1 = (4.*k1*(n1_bar**2))/(dV*((4.*k1*n1_bar)+k2))
Var_n2 = (1./4.)*Var_n1

print("n1 = %e\tn2 = %e" % (A-2.*sol1,sol1))
print("n1 = %e\tn2 = %e" % (A-2.*sol2,sol2))

# print variance
print("Var1 = %e\tVar2 = %e" % (Var_n1,Var_n2))
