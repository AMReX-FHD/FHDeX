from itertools import combinations as comb
import numpy as np

# works up to 11x11 (as far as I can tell)

def solve_system(ads,des,dens):
    R = np.zeros(len(ads),dtype=np.longdouble)
    for i in range(0,len(dens)):
        R[i] = (des[i] / (ads[i]*dens[i])) + 1
    n = len(R)
    x = np.zeros(n)
    res = 0
    if n == 1:
        return np.array([1/float(R[0])])
    elif n == 2:
        denom = (R[0] * R[1] - 1)
        x1 = float(R[1] - 1) / denom
        x2 = float(R[0] - 1) / denom
        return np.array([x1,x2])
    else:
        denom = (n-2)*np.sum(R) + (-1)**n * reduce(lambda z, y: z*y, R) - (n-1)
        for k in range(2,n-1):
            c = list(comb(R,k))
            res = 0
            for tup in c:
                prod = reduce(lambda z, y: z*y, tup)
                res += prod
            denom += (-1)**(k+1) * (n - (k+1))*res
        for i in range(0,n):
            r = np.delete(np.array(R),i)
            numer = np.sum(r) - 1
            m = n-1
            for k in range(2,m+1):
                c = list(comb(r,k))
                res = 0
                for tup in c:
                    prod = reduce(lambda z, y: z*y, tup)
                    res += prod
                numer += (-1)**(k+1) * res
            x[i] = float(numer) / denom
    return x
