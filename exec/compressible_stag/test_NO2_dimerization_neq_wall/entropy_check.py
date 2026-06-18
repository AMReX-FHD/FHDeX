import numpy as np
import math
import entropy_check_funcs



ns = 2

# read data

infile = "res.zavg3"

[z_in,rho_in,rhoE_in,rhoYk0_in,rhoYk1_in,T_in,p_in,Yk0_in,Yk1_in,Xk0_in,Xk1_in] = np.loadtxt(infile,usecols=(0,1,3,5,7,9,11,13,15,17,19),unpack=True)

nz = len(z_in)
dz = z_in[1]-z_in[0]

cu0 = rho_in
cu4 = rhoE_in
cu5 = rhoYk0_in
cu6 = rhoYk1_in

prim0 = rho_in
prim4 = T_in
prim5 = p_in
prim6 = Yk0_in
prim7 = Yk1_in
prim8 = Xk0_in
prim9 = Xk1_in



# compute transport coefficients
[eta,kappa,zeta,Dij,chi] = entropy_check_funcs.ComputeTransportCoefficients(prim0,prim4,prim5,prim6,prim7)

# compute chemical potentials and chemical production rates
[mu1,mu2,Omega1,Omega2] = entropy_check_funcs.ComputeChemPotChemProdRate(cu5,cu6,prim4)

# output for additional layer info
outfile1 = "res.layer_info"
np.savetxt(outfile1,np.column_stack((z_in,eta,kappa,zeta,Dij[:,0],Dij[:,1],Dij[:,2],Dij[:,3],chi[:,0],chi[:,1],mu1,mu2,Omega1,Omega2)))
print("** %s generated" % outfile1)



# compute fluxes Q and Fk
[Qz,F1z,F2z,QQz] = entropy_check_funcs.ComputeQzFkz(dz,kappa,Dij,chi,prim4,prim5,prim6,prim7,prim8,prim9)

# output for fluxes
outfile2 = "res.fluxes"
np.savetxt(outfile2,np.column_stack((z_in,Qz,F1z,F2z,QQz)))
print("** %s generated" % outfile2)



# div(Q/T)
Tface = np.zeros(nz)
for k in range(1,nz):
    Tface[k] = 0.5*(prim4[k-1]+prim4[k])

term1 = np.zeros(nz)
for k in range(1,nz-1):
    tmp2 = Qz[k+1]/Tface[k+1]
    tmp1 = Qz[k]/Tface[k]
    term1[k] = (tmp2-tmp1)/dz

# sum div(muk*Fk/T)
mu1face = np.zeros(nz)
mu2face = np.zeros(nz)
for k in range(1,nz):
    mu1face[k] = 0.5*(mu1[k-1]+mu1[k])
    mu2face[k] = 0.5*(mu2[k-1]+mu2[k])

term2 = np.zeros(nz)
for k in range(1,nz-1):
    tmp2 = mu1face[k+1]*F1z[k+1]/Tface[k+1]
    tmp1 = mu1face[k]*F1z[k]/Tface[k]
    term2[k] += (tmp2-tmp1)/dz
    tmp2 = mu2face[k+1]*F2z[k+1]/Tface[k+1]
    tmp1 = mu2face[k]*F2z[k]/Tface[k]
    term2[k] += (tmp2-tmp1)/dz

# (Q.gradT)/T^2
gradTz = np.zeros(nz)
for k in range(1,nz):
    gradTz[k] = (prim4[k]-prim4[k-1])/dz

term3 = np.zeros(nz)
for k in range(1,nz-1):
    tmp1 = Qz[k]*gradTz[k]/Tface[k]/Tface[k]
    tmp2 = Qz[k+1]*gradTz[k+1]/Tface[k+1]/Tface[k+1]
    term3[k] = 0.5*(tmp1+tmp2)

# sum Fk.div(muk/T)
gradmu1Tz = np.zeros(nz)
gradmu2Tz = np.zeros(nz)
for k in range(1,nz):
    gradmu1Tz[k] = (mu1[k]/prim4[k]-mu1[k-1]/prim4[k-1])/dz
    gradmu2Tz[k] = (mu2[k]/prim4[k]-mu2[k-1]/prim4[k-1])/dz

term4 = np.zeros(nz)
for k in range(1,nz-1):
    tmp1 = F1z[k]*gradmu1Tz[k]
    tmp2 = F1z[k+1]*gradmu1Tz[k+1]
    term4[k] = 0.5*(tmp1+tmp2)
    tmp1 = F2z[k]*gradmu2Tz[k]
    tmp2 = F2z[k+1]*gradmu2Tz[k+1]
    term4[k] += 0.5*(tmp1+tmp2)

# (sum muk*Mk*Omegak)/T
M1 = 46.0055
M2 = 92.0110
term5 = np.zeros(nz)
for k in range(1,nz-1):
    term5[k] = (mu1[k]*M1*Omega1[k]+mu2[k]*M2*Omega2[k])/prim4[k]



# div(F1)
term11 = np.zeros(nz)
for k in range(1,nz-1):
    term11[k] = (F1z[k+1]-F1z[k])/dz

# M1*Omega1
term12 = np.zeros(nz)
for k in range(1,nz-1):
    term12[k] = M1*Omega1[k]

# div(Qz)
term13 = np.zeros(nz)
for k in range(1,nz-1):
    term13[k] = (Qz[k+1]-Qz[k])/dz

# (div.Q)/T
term14 = np.zeros(nz)
for k in range(1,nz-1):
    term14[k] = (Qz[k+1]-Qz[k])/dz/prim4[k]

# (sum muk*(div(Fk)-Mk*Omegak))/T
term15 = np.zeros(nz)
for k in range(1,nz-1):
    tmp = (F1z[k+1]-F1z[k])/dz
    term15[k] = mu1[k]*(tmp-M1*Omega1[k])/prim4[k]
    tmp = (F2z[k+1]-F2z[k])/dz
    term15[k] += mu2[k]*(tmp-M2*Omega2[k])/prim4[k]



# output for final terms
outfile3 = "res.entropy_check"
np.savetxt(outfile3,np.column_stack((z_in,term1,term2,term3,term4,term5,term11,term12,term13,term14,term15,mu1,mu2,Omega1,Omega2)))
print("** %s generated" % outfile3)
