import numpy as np

# parameters
dx = 1.4976e-4 / 16
k_B = 1.38064852e-16
T_0 = 700.0
V = dx**3
dt = 1e-12

cv1 = 8.160776e6   # CO
cv2 = 3.121786e6   # Ar
cvmix = 0.5 * cv1 + 0.5 * cv2

pres = 1.013250e6
p1 = 0.587846 * pres
p2 = 0.412154 * pres

theta1 = 7.466307e-1
surf_site_num_dens = 1.027285e15

# if you want the number of sites on one wall-cell face:
Ntot = surf_site_num_dens * dx * dx

rho0 = 5.733118e-4
rho1 = rho0 * 0.5
rho2 = rho0 * 0.5

# standard deviations
theta_std = np.sqrt(theta1 * (1.0 - theta1) / Ntot)
temp_std = np.sqrt(k_B * T_0**2 / (rho0 * V * cvmix))
pres_std = np.sqrt(k_B * T_0 * p1 / V)

rng = np.random.default_rng(12345)

def theta_dist(n):
    return rng.normal(loc=theta1, scale=theta_std, size=n)

def temp_dist(n):
    return rng.normal(loc=T_0, scale=temp_std, size=n)

def pres_dist(n):
    return rng.normal(loc=p1, scale=pres_std, size=n)

n = 100000
theta_samples = theta_dist(n)
tempratio_samples = temp_dist(n)/T_0
pres_samples = pres_dist(n)

ads_rate = 1.831671e2 * pres_samples * (1 - theta_samples) * (tempratio_samples)**(-0.5)
des_rate = 3.702336e7 * theta_samples

data = np.column_stack([pres_samples, tempratio_samples, theta_samples, ads_rate, des_rate])

np.savetxt("mfsurfchem_samples.txt", data, header="pressure temp theta rads rdes")
