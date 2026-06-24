import numpy as np
import math



# constants

k_B = 1.38065e-16
Runiv = 8.31446e7
pi = 3.1415926535897932



# species properties

ns = 2

molmass = np.array([46.0055,92.0110])
diameter = np.array([3.765e-08,4.621e-08])

e0 = np.array([4.685624e+09,-1.730782e+09])
dof = np.array([7.279586,18.049072])

hcp = np.zeros(ns)
for i in range(ns):
    hcp[i] = 0.5*(2.+dof[i])*Runiv/molmass[i]

Tst = 3.500000e+02
Pst = 1.000000e+06
sst = np.array([5.349074e+07,3.447942e+07])

# reaction parameters

k1 = 4.065753e+11
k2 = 6.236899e+07
a1 = -5.337990e+10
a2 = 5.370000e+11
b1 = 0.644950
b2 = -1.100000



# functions

def GetEnthalpies(T):

    global e0
    global hcp

    return e0+hcp*T



def GetChemicalPotentials(T,ck):

    global Runiv
    global molmass
    global Tst
    global Pst
    global sst

    hk = GetEnthalpies(T)

    return hk - T*sst - T*math.log(T/Tst)*hcp + np.divide(Runiv*T*np.log(Runiv*T*ck/Pst),molmass)



def GetChemicalProductionRates(T,ck):

    global Runiv
    global Tst
    global k1
    global k2
    global a1
    global a2
    global b1
    global b2

    rate1 = k1*math.exp(-a1/Runiv*(1/T-1/Tst))*math.pow(T/Tst,b1)*ck[0]*ck[0]
    rate2 = k2*math.exp(-a2/Runiv*(1/T-1/Tst))*math.pow(T/Tst,b2)*ck[1]

    Omega1 = -2*rate1+2*rate2
    Omega2 = rate1-rate2

    return np.array([Omega1,Omega2])



def D_GIO1(Dbin,Ykp,Xkp,D_tilde,kloc):

    global ns

    Diff_ij = np.zeros(ns*ns)
    Zmat = np.zeros(ns*ns)
    Pmat = np.zeros(ns*ns)
    Jmat = np.zeros(ns*ns)
    PJ = np.zeros(ns*ns)
    matrix1 = np.zeros(ns*ns)
    matrix2 = np.zeros(ns*ns)

    Di = np.zeros(ns)
    Minv = np.zeros(ns)
    Mmat = np.zeros(ns)

    jmax = 3

    # Find Di matrix
    for i in range(ns):
        term2 = 0.
        for j in range(ns):
            if (j != i):
                term2 += Xkp[j]/Dbin[i,j]
        Di[i] = (1.-Ykp[i])/term2

    # Compute Mmat and Minv
    for i in range(ns):
        Mmat[i] = Xkp[i]/Di[i]
        Minv[i] = Di[i]/Xkp[i]

    # Compute P matrix
    for i in range(ns):
        for j in range(ns):
            Pmat[i*ns+j] = -Ykp[j]
            if (i == j):
                Pmat[i*ns+j] += 1.

    # Compute Deltamat
    for i in range(ns):
        for j in range(ns):
            if (i == j):
                term1 = 0.
                for k in range(ns):
                    if (k != i):
                        term1 += Xkp[i]*Xkp[k]/Dbin[i,k]
                Deltamat = term1
            else:
                Deltamat = -Xkp[i]*Xkp[j]/Dbin[i,j]
            Zmat[i*ns+j] = -Deltamat

    # Compute Zmat
    for i in range(ns):
        Zmat[i*ns+i] += Mmat[i]

    # Compute Jmat
    for i in range(ns):
        for j in range(ns):
            Jmat[i*ns+j] = Minv[i]*Zmat[i*ns+j]

    # Compute PJ
    for i in range(ns):
        for j in range(ns):
            PJ[i*ns+j] = 0.
            for k in range(ns):
                PJ[i*ns+j] += Pmat[i*ns+k]*Jmat[k*ns+j]

    # Compute P M^-1 Pt; store it in matrix2
    for i in range(ns):
        for j in range(ns):
            scr = 0.
            for k in range(ns):
                scr += Pmat[i*ns+k]*Minv[k]*Pmat[j*ns+k]
                # notice the change in indices for Pmat to represent Pmat^t
            matrix2[i*ns+j] = scr
            Diff_ij[i*ns+j] = scr

    if (jmax > 0):
        for ll in range(jmax):
            # matrix1=0
            for i in range(ns):
                for j in range(ns):
                    scr = 0.
                    for k in range(ns):
                        scr += PJ[i*ns+k]*Diff_ij[k*ns+j]
                    matrix1[i*ns+j] = scr+matrix2[i*ns+j]
            for i in range(ns):
                for j in range(ns):
                    Diff_ij[i*ns+j] = matrix1[i*ns+j]

    # Compute D_tilde
    for i in range(ns):
        for j in range(ns):
            # indicies on LHS are switched due to C++/F90 memory layout
            # this is the only location this needs to happen I believe
            D_tilde[kloc,j*ns+i] = Diff_ij[i*ns+j]*Ykp[i]



# Transport Coefficients from Giovangigli

def IdealMixtureTransportGIO(kloc,density,temperature,pressure,Yk,diff_ij,chitil):

    global pi
    global k_B
    global Runiv
    global ns
    global molmass
    global diameter

    mu = np.zeros((ns,ns))
    diam = np.zeros((ns,ns))
    Dbin = np.zeros((ns,ns))
    amat = np.zeros((ns,ns))

    molecular_mass = np.zeros(ns)
    Xk = np.zeros(ns)
    Xkp = np.zeros(ns)
    Ykp = np.zeros(ns)
    eta1 = np.zeros(ns)
    rhs = np.zeros(ns)
    Bi = np.zeros(ns)

    # compute molecular_mass by dividing molmass by Avogadro's
    for n in range(ns):
        molecular_mass[n] = molmass[n]*(k_B/Runiv)

    # mixture molecular weight
    mbar = 0.0
    for n in range(ns):
        mbar += Yk[n]/molecular_mass[n]
    mbar = 1.0/mbar

    # mole fracions
    Xksum = 0.0
    for n in range(ns):
        Xk[n] = Yk[n]*mbar/molecular_mass[n]
        Xksum += Xk[n]

    # Perturb Xk for zero or negative values -- store in Xkp
    for n in range(ns):
        Xkp[n] = Xk[n] + 1.0e-16*(Xksum/float(ns) - Xk[n])

    # Perturbed mass fractions Ykp
    for n in range(ns):
        Ykp[n] = Xkp[n]*molecular_mass[n]/mbar

    for i in range(ns):
        for j in range(ns):
            diam[i,j] = 0.5*(diameter[i] + diameter[j])
            mu[i,j] = molecular_mass[i]*molecular_mass[j]/(molecular_mass[i] + molecular_mass[j])

    # Binary Diffusion Coefficients
    kT = k_B*temperature
    for i in range(ns):
        for j in range(ns):
            Dbin[i,j] = (3.0/16.0)*math.sqrt(2.0*pi*kT*kT*kT/mu[i,j])/(pressure*pi*diam[i,j]*diam[i,j])

    AKL = 1.0
    BKL = 1.0
    CKL = 1.0

    # Viscosity
    for i in range(ns):
        eta1[i] = 5.0/(16.0*diam[i,i]*diam[i,i])*math.sqrt(molecular_mass[i]*k_B*temperature/pi)

    for i in range(ns):
        rhs[i] = Xkp[i]
        Bi[i] = rhs[i]
        for j in range(ns):
            if (i==j):
                amat[i,i] = Xkp[i]*Xkp[i]/eta1[i]
                for k in range(ns):
                    if (k!=i):
                        amat[i,i] += 2.0*kT*Xkp[i]*Xkp[k]/(pressure*Dbin[i,k]*(molecular_mass[i]+molecular_mass[k]))*(1.0+0.6*molecular_mass[k]/molecular_mass[i]*AKL)
            else:
                amat[i,j] = 2.0*kT*Xkp[i]*Xkp[j]/(pressure*Dbin[i,j]*(molecular_mass[i]+molecular_mass[j]))*(0.6*AKL-1.0)

    for j in range(ns-1):
        for i in range(j+1,ns):
            fact1 = amat[i,j]/amat[j,j]
            for k in range(ns):
                amat[i,k] -= fact1*amat[j,k]
            Bi[i] -= fact1*Bi[j]
    Bi[ns-1] /= amat[ns-1,ns-1]

    for i in range(ns-2,-1,-1):
        sum1 = 0.0
        for j in range(ns-1,i,-1):
            sum1 += amat[i,j]*Bi[j]
        Bi[i] = (Bi[i] - sum1)/amat[i,i]

    eta = 0.0
    for i in range(ns):
        eta += Xkp[i]*Bi[i]

    # zero bulk viscosity for monoatomic HS
    zeta = 0.

    # Thermal Conductivity
    for i in range(ns):
        rhs[i] = 2.5*Xkp[i]
        Bi[i] = rhs[i]
        for j in range(ns):
            if (i==j):
                amat[i,i] = Xkp[i]*Xkp[i]/Dbin[i,j]*2.0*AKL
                for k in range(ns):
                    if (k!=i):
                        amat[i,i] += Xkp[i]*Xkp[k]/Dbin[i,k]*mu[i,k]/(molecular_mass[i]+molecular_mass[k])*(7.5*molecular_mass[i]/molecular_mass[k]+6.25*molecular_mass[k]/molecular_mass[i]-3.0*molecular_mass[k]/molecular_mass[i]*BKL+4.0*AKL)
            else:
                amat[i,j] = -1.0*Xkp[i]*Xkp[j]/Dbin[i,j]*mu[i,j]/(molecular_mass[i]+molecular_mass[j])*(55.0/4.0-3.0*BKL-4.0*AKL)

    for j in range(ns-1):
        for i in range(j+1,ns):
            fact1 = amat[i,j]/amat[j,j]
            for k in range(ns):
                amat[i,k] -= fact1*amat[j,k]
            Bi[i] -= fact1*Bi[j]
    Bi[ns-1] /= amat[ns-1,ns-1]

    for i in range(ns-2,-1,-1):
        sum1 = 0.0
        for j in range(ns-1,i,-1):
            sum1 += amat[i,j]*Bi[j]
        Bi[i] = (Bi[i] - sum1)/amat[i,i]

    kappa = 0.0
    for i in range(ns):
        kappa += rhs[i]*Bi[i]
    kappa *= pressure/temperature

    # Thermal Diffusion Ratios
    for i in range(ns):
        chitil[kloc,i] = 0.0
        for j in range(ns):
            if (i==j):
                amat[i,i] = 0.0
                for k in range(ns):
                    if (k!=i):
                        amat[i,i] -= 0.5*Xkp[k]/Dbin[i,k]*molecular_mass[k]/(molecular_mass[i]+molecular_mass[k])*(6.0*CKL-5.0)
            else:
                amat[i,j] = 0.5*Xkp[j]/Dbin[i,j]*molecular_mass[i]/(molecular_mass[i]+molecular_mass[j])*(6.0*CKL-5.0)

    for i in range(ns):
        for j in range(ns):
            chitil[kloc,i] += amat[i,j]*Bi[j]

    D_GIO1(Dbin,Ykp,Xkp,diff_ij,kloc)

    return [eta,kappa,zeta]



def ComputeTransportCoefficients(prim0,prim4,prim5,prim6,prim7):

    global ns

    nz = len(prim0)

    eta = np.zeros(nz)
    kappa = np.zeros(nz)
    zeta = np.zeros(nz)

    Dij = np.zeros((nz,ns*ns))
    chi = np.zeros((nz,ns))

    for k in range(nz):
        Yk = np.array([prim6[k],prim7[k]])
        [eta[k],kappa[k],zeta[k]] = IdealMixtureTransportGIO(k,prim0[k],prim4[k],prim5[k],Yk,Dij,chi)

        # want this multiplied by rho for all times
        for kk in range(ns):
            for ll in range(ns):
                n = kk*ns + ll
                Dij[k,n] *= prim0[k]

    return [eta,kappa,zeta,Dij,chi]



def ComputeChemPotChemProdRate(cu5,cu6,prim4):

    global molmass

    nz = len(cu5)

    mu1 = np.zeros(nz)
    mu2 = np.zeros(nz)
    Omega1 = np.zeros(nz)
    Omega2 = np.zeros(nz)

    for k in range(nz):

        ck = np.array([cu5[k]/molmass[0],cu6[k]/molmass[1]])

        muk = GetChemicalPotentials(prim4[k],ck)
        mu1[k] = muk[0]
        mu2[k] = muk[1]

        Omegak = GetChemicalProductionRates(prim4[k],ck)
        Omega1[k] = Omegak[0]
        Omega2[k] = Omegak[1]

    return [mu1,mu2,Omega1,Omega2]



def ComputeQzFkz(dz,kappa,Dij,chi,prim4,prim5,prim6,prim7,prim8,prim9):

    global Runiv
    global ns
    global molmass

    nz = len(kappa)

    Qz = np.zeros(nz)
    QQz = np.zeros(nz)      # only heat conduction
    F1z = np.zeros(nz)
    F2z = np.zeros(nz)

    for k in range(1,nz):

        # compute dk

        dk = np.zeros(ns)
        meanXk = np.zeros(ns)
        meanYk = np.zeros(ns)
        soret = np.zeros(ns)

        meanT = 0.5*(prim4[k-1]+prim4[k])
        meanP = 0.5*(prim5[k-1]+prim5[k])

        term1 = (prim8[k]-prim8[k-1])/dz
        meanXk[0] = 0.5*(prim8[k-1]+prim8[k])
        meanYk[0] = 0.5*(prim6[k-1]+prim6[k])
        term2 = (meanXk[0]-meanYk[0])*(prim5[k]-prim5[k-1])/dz/meanP
        dk[0] = term1 + term2
        ChiX = 0.5*(chi[k-1,0]*prim8[k-1]+chi[k,0]*prim8[k])
        soret[0] = ChiX*(prim4[k]-prim4[k-1])/dz/meanT

        term1 = (prim9[k]-prim9[k-1])/dz
        meanXk[1] = 0.5*(prim9[k-1]+prim9[k])
        meanYk[1] = 0.5*(prim7[k-1]+prim7[k])
        term2 = (meanXk[1]-meanYk[1])*(prim5[k]-prim5[k-1])/dz/meanP
        dk[1] = term1 + term2
        ChiX = 0.5*(chi[k-1,1]*prim9[k-1]+chi[k,1]*prim9[k])
        soret[1] = ChiX*(prim4[k]-prim4[k-1])/dz/meanT

        # compute Fk (based on Eqn. 2.5.24, Giovangigli's book)

        Fk = np.zeros(ns)

        for kk in range(ns):
            Fk[kk] = 0.
            for ll in range(ns):
                Fks = 0.5*(Dij[k-1,ll*ns+kk]+Dij[k,ll*ns+kk])*(dk[ll]+soret[ll])
                Fk[kk] -= Fks

        F1z[k] = Fk[0]
        F2z[k] = Fk[1]

        # compute Q (based on Eqn. 2.5.25, Giovangigli's book)

        hk = GetEnthalpies(meanT)

        Q5 = 0.0
        for n in range(ns):
            Q5s = (hk[n] + 0.5*Runiv*meanT*(chi[k,n]+chi[k,n])/molmass[n])*Fk[n]
            Q5 += Q5s
        Qz[k] = Q5

        kzp = 0.5*(kappa[k-1]+kappa[k])
        Qz[k] -= kzp*(prim4[k]-prim4[k-1])/dz

        QQz[k] -= kzp*(prim4[k]-prim4[k-1])/dz

    return [Qz,F1z,F2z,QQz]
