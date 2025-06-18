import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

all_theta = []

for i in range(1, 401):
    filename = f"data{i}.txt"
    theta = np.loadtxt(filename,usecols = (-4))
    all_theta.append(theta)  

all_theta = np.array(all_theta)

### plot mean of sample trajectories of theta
plt.figure(figsize=(10, 5))
plt.plot(all_theta.mean(axis=0))
plt.title(r"Mean of sample trajectories of $\theta$")
plt.savefig("Sample trajectories of theta (mean)")
plt.close()

### plot variance of sample trajectories of theta
plt.figure(figsize=(10, 5))
plt.plot(all_theta.var(axis=0))
plt.title(r"Variance of sample trajectories of $\theta$")
plt.savefig("Sample trajectories of theta (var)")
plt.close()

### plot 20 sample trajectories of theta
plt.figure(figsize=(10, 5))
for i in range(20):
    plt.plot(all_theta[i, :])
plt.xlabel('time step')
plt.ylabel(r'$\theta$')
plt.savefig("Sample trajectories of theta")
plt.close()

######### define ACF based on bruteforce method #########
def acf_numpy(x, lags=None):
    """
    Compute the unnormalized autocovariance function for an array of shape (n_trajectories, n_timepoints)

    Parameters:
        x: ndarray of shape (N, T) -- N trajectories, T timepoints
        lags: int or array-like or None

    Returns:
        corr: ndarray of shape (len(lags),)
    """
    x = np.asarray(x)
    x = x - x.mean(axis=1, keepdims=True)  # center each trajectory

    N, T = x.shape
    if lags is None:
        lags = np.arange(T)
    elif isinstance(lags, int):
        lags = np.arange(lags)
    else:
        lags = np.array(lags)

    corr = np.zeros(len(lags))
    for i, lag in enumerate(lags):
        if lag >= T:
            continue
        u = x[:, :T - lag]
        v = x[:, lag:]
        corr[i] = np.sum(u * v) / (N * (T - lag))

    return corr



### plot acf of theta
print("generating ACF of theta")
acf_vals = acf_numpy(all_theta[:,::100])
plt.figure(figsize=(10, 5))
plt.plot(acf_vals)
plt.xlabel('lags')
plt.ylabel('ACF')
plt.title(r"ACF of $\theta$")
plt.grid(True)
plt.tight_layout()
plt.savefig('acf_theta', dpi = 300, bbox_inches = 'tight')
plt.close()

