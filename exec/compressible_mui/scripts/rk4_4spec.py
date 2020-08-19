import numpy as np
import matplotlib.pyplot as plt

def solve(ads,des,dens,t0,T,spec_init,N_0,N):
    
    # calculate h (step-size)
    h = (T - t0)/N

    # define sol'n vectors
    t = np.zeros(N+1)
    spec1 = np.zeros(N+1,dtype=np.longdouble)
    spec2 = np.zeros(N+1,dtype=np.longdouble)
    spec3 = np.zeros(N+1,dtype=np.longdouble)
    spec4 = np.zeros(N+1,dtype=np.longdouble)

    t[0] = t0
    spec1[0] = spec_init[0]
    spec2[0] = spec_init[1]
    spec3[0] = spec_init[2]
    spec4[0] = spec_init[3]

    f = lambda x,y,z,u: (1 - x - y - z - u)*ads[0]*dens[0] - des[0]*x
    g = lambda x,y,z,u: (1 - x - y - z - u)*ads[1]*dens[1] - des[1]*y
    p = lambda x,y,z,u: (1 - x - y - z - u)*ads[2]*dens[2] - des[2]*z
    q = lambda x,y,z,u: (1 - x - y - z - u)*ads[3]*dens[3] - des[3]*u

    F = [f,g,p,q]
    num_F = len(F)

    k = np.zeros([num_F,4],dtype=np.longdouble)

    for i in range(0,N):
        for j in range(0,num_F):
            k[j,0] = h*F[j](spec1[i],spec2[i],spec3[i],spec4[i])
        for j in range(0,num_F):
            k[j,1] = h*F[j](spec1[i] + k[j,0]/2,spec2[i] + k[j,0]/2,spec3[i] + k[j,0]/2,spec4[i] + k[j,0]/2)
        for j in range(0,num_F):   
            k[j,2] = h*F[j](spec1[i] + k[j,1]/2,spec2[i] + k[j,1]/2,spec3[i] + k[j,1]/2,spec4[i] + k[j,1]/2)
        for j in range(0,num_F):
            k[j,3] = h*F[j](spec1[i] + k[j,2],spec2[i] + k[j,2],spec3[i] + k[j,2],spec4[i] + k[j,2])
        
        spec1[i+1] = spec1[i] + (1./6)*(k[0,0] + 2*k[0,1] + 2*k[0,2] + k[0,3])
        spec2[i+1] = spec2[i] + (1./6)*(k[1,0] + 2*k[1,1] + 2*k[1,2] + k[1,3])
        spec3[i+1] = spec3[i] + (1./6)*(k[2,0] + 2*k[2,1] + 2*k[2,2] + k[2,3])
        spec4[i+1] = spec4[i] + (1./6)*(k[3,0] + 2*k[3,1] + 2*k[3,2] + k[3,3])
        t[i+1] = t[i] + h
    
    spec1 = N_0 * spec1
    spec2 = N_0 * spec2
    spec3 = N_0 * spec3
    spec4 = N_0 * spec4
    
    return t,spec1,spec2,spec3,spec4
