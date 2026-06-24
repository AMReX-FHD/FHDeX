import numpy as np
import matplotlib.pyplot as plt
import sys

def getData(data):
  #extract data from a flux file

  #declare lists
  lFlux = []
  rFlux = []
  netFlux = []

  #counter for first few arguments
  counter = 0

  #open file with data
  file = open(data,"r")

  #Read in each line
  for line in file:

    entries = line.split()

    if counter == 0:
      dt = float(entries[0])
    elif counter == 1:
      lN = int(entries[0])
      rN = int(entries[1])
    elif counter == 2:
      lT = float(entries[0])
      rT = float(entries[1])
    else:
      lFlux.append(int(entries[0]))
      rFlux.append(int(entries[1]))

    counter = counter + 1

  #close file
  file.close()

  #get the net flux - particles going right - left
  N = len(lFlux)
  for i in xrange(0,N):
    netFlux.append((rFlux[i]) - (lFlux[i]))

  return dt, lN, rN, lT, rT, lFlux, rFlux, netFlux

def makeHist(flux, x, D):
  #make a histogram with the data in flux

  #make some histograms
  hist, bin_edges = np.histogram(flux)

  #plot some histograms
  b = 50*np.max(flux)-np.min(flux)*50
  # An "interface" to matplotlib.axes.Axes.hist() method
  n, bins, patches = plt.hist(x=flux, bins=int(np.ceil(b+1)), color='#0504aa',
                              alpha=0.7, rwidth=0.85, density=1)
  plt.grid(axis='y', alpha=0.75)
  plt.xlabel('Value')
  plt.ylabel('Frequency')
  plt.title('Flux Histogram')
  maxfreq = n.max()
  print bins
  # Set a clean upper y-axis limit.
  #plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

  #plot the theoretical distribution
  plt.plot(x,D)

  #show plots on one axes
  plt.show()

  return

def getAverageFlux(N, ts):
  #get the average flux as a function of timestep

  #construct storage for data - first element is 0
  cumFlux = [0.0 for _ in range(ts+1)]

  #base for data location
  base = "samples/fluxes"

  #loop over all the data
  for i in xrange(0,N):
    #get the data
    loc = base + str(i) + ".txt"
    dt, lN, rN, lT, rT, lFlux, rFlux, netFlux = getData(loc)

    #compute a cumulative sum of netFlux
    cs = np.cumsum(netFlux)
    cs = np.insert(cs,0,0) #append a zero as the first element

    #add the cumsum to cumFlux
    cumFlux = np.add(cumFlux, cs)

  #divide cumFlux by number of samples
  Ninv = 1.0 / N
  #cumFlux.astype(float)
  cumFlux[:] = [x * Ninv for x in cumFlux]

  return cumFlux

def plotAverageFlux(avgFlux, N, ts, dt, eBarFlag, theory):
  #plot the average flux

  #get the theory curve
  t = xrange(0,ts+1)
  Th = np.ones(ts+1)*theory

  if (eBarFlag == 1):
    eBars = []

    #get error bars every 500th time step
    for i in xrange(0,ts+1):
      if (i % 500 == 0):
        eBars.append(getErrorBar(N, i))
        print "Got error bar at ", i
      else:
        eBars.append(0)

    #plot the data
    plt.errorbar(t, avgFlux, eBars, errorevery = 500)
  else:
    #plot the data
    plt.plot(t, avgFlux)

  plt.plot(t,Th)
  plt.title('Average Flux over time')
  plt.xlabel('Time Step')
  plt.ylabel('Flux')
  plt.legend(['DSMC', 'Theory Average'])
  plt.show()

  return

def getErrorBar(N, time):
  #get error bar on average flux at time step = time

  #construct storage for data - first element is 0
  cumFlux = [0.0 for _ in range(ts+1)]

  #base for data location
  base = "samples/fluxes"

  #hold all cumFLuxes at given time step
  measurements = []

  #loop over all the data
  for i in xrange(0,N):
    #get the data
    loc = base + str(i) + ".txt"
    dt, lN, rN, lT, rT, lFlux, rFlux, netFlux = getData(loc)

    #compute cumsum
    cs = np.cumsum(netFlux)
    cs = np.insert(cs,0,0) #append a zero as the first element

    measurements.append(cs[time]);

  eBar = np.std(measurements)

  return eBar

def estimateRates(lFlux, rFlux, netFlux, lN, rN, dt):
  #estimate a rate of crossing using data
  #A->B is rFlux, B->A is lFLux. N is original amount

  #create number of particles in each reservoir as fn of time
  lcount = lN
  rcount = rN
  lP = []
  rP = []
  lP.append(lcount)
  rP.append(rcount)
  ts = len(netFlux)
  for i in xrange(0,ts):
    lcount = lcount - netFlux[i]
    rcount = rcount + netFlux[i]
    lP.append(lcount)
    rP.append(rcount)

  #get rates
  lrate = []
  rrate = []
  ltime = 0
  rtime = 0
  for i in xrange(0,ts):
    ltime = ltime + 1
    rtime = rtime + 1
    if rFlux[i] > 0:
      #print rP[i-1], rP[i], rP[i+1]
      rrate.append((float(rFlux[i])/lP[i])/(rtime*dt))
      rtime = 0
    if lFlux[i] > 0:
      lrate.append((float(lFlux[i])/rP[i])/(ltime*dt))
      ltime = 0

  return lrate, rrate

def estimateRatesAll(N):
  #estimate the rates using every trial

  lrates = []
  rrates = []

  #base for data location
  base = "samples/fluxes"

  #loop over all the data
  for i in xrange(0,N):
    #get the data
    loc = base + str(i) + ".txt"
    dt, lN, rN, lT, rT, lFlux, rFlux, netFlux = getData(loc)

    lrate, rrate = estimateRates(lFlux, rFlux, netFlux, lN, rN, dt)
    lsample = np.mean(lrate)
    rsample = np.mean(rrate)

    lrates.append(lsample)
    rrates.append(rsample)

  return lrates, rrates

def fluxSamples(N,t):
  #get all flux samples at time t

  samples = []

  #base for data location
  base = "samples/fluxes"

  #loop over all the data
  for i in xrange(0,N):
    #get the data
    loc = base + str(i) + ".txt"
    dt, lN, rN, lT, rT, lFlux, rFlux, netFlux = getData(loc)

    cs = np.cumsum(netFlux)
    cs = np.insert(cs,0,0) #append a zero as the first element

    samples.append(cs[t]/50.0)

  return samples

def theoryDistribution(xpts, time, pL, pR, tL, tR):
  #evaluate the theoretical distribution for given time and parameters

  P = pL + pR                    #total number of particles
  c = np.sqrt(float(tR)/tL)             #ratio of birth rate to death rate
  x = np.linspace(-pR+0.01, P-pR-0.01, xpts)  #x domain
  pDist = []
  dx = x[1]-x[0]

  #compute the distribution pointwise
  Z = 0.0
  for i in xrange(0,xpts):
    xp = x[i]
    phi = (xp+pR)*np.log(1/c*((xp+pR)/(P-(xp+pR))))+P*np.log(P-(xp+pR))+P*np.log((1+c)/P)
    Zi = np.exp(-time*phi)
    pDist.append(Zi)
    Z = Z + Zi*dx

  #normalize the distribution
  for i in xrange(0,xpts):
    pDist[i] = pDist[i]/Z

  return x, pDist








if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: python fluxStats.py <numSamples> <numTimeSteps> <plotEbar>")

    #store inputs
    N = int(sys.argv[1])
    ts = int(sys.argv[2])
    eBarFlag = int(sys.argv[3])

    #get an example run for parameters
    dt, lN, rN, lT, rT, lFlux, rFlux, netFlux = getData("samples/fluxes5.txt")
    print "Temperatures:     ", lT, rT
    print "Starting numbers: ", lN, rN

    #estimate rates
    #lrate, rrate = estimateRates(lFlux, rFlux, netFlux, lN, rN, dt)

    #lrates, rrates = estimateRatesAll(N)
    #makeHist(lrates)
    #makeHist(rrates)
    #lm = np.mean(lrates)
    #rm = np.mean(rrates)
    #print lm, rm

    #compute expected flux from theory using temperatures
    theory = np.sqrt(lT)/(np.sqrt(lT)+np.sqrt(rT)) * (lN+rN) - rN
    print "The theoretical mean is: ", theory

    #get the average flux as a function of timestep
    print "Getting Average Flux..."
    avgFlux = getAverageFlux(N, ts)
    plotAverageFlux(avgFlux, N, ts, dt, eBarFlag, theory)

    #get a histogram of net flux at a chosen time
    time = 10000
    P = lN+rN
    print "Plotting histogram at time step ", time
    samples = fluxSamples(N, time)
    x, D = theoryDistribution(100, 50, lN/50.0, rN/50.0, rT, lT)
    makeHist(samples, x, D)








