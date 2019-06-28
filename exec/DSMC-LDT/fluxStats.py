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
      dt = entries[0]
    elif counter == 1:
      lN = entries[0]
      rN = entries[1]
    elif counter == 2:
      lT = entries[0]
      rT = entries[1]
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

def makeHist(flux):
  #make a histogram with the data in flux

  #make some histograms
  hist, bin_edges = np.histogram(flux)

  #plot some histograms
  # An "interface" to matplotlib.axes.Axes.hist() method
  n, bins, patches = plt.hist(x=flux, bins='auto', color='#0504aa',
                              alpha=0.7, rwidth=0.85)
  plt.grid(axis='y', alpha=0.75)
  plt.xlabel('Value')
  plt.ylabel('Frequency')
  plt.title('Flux Histogram')
  maxfreq = n.max()
  # Set a clean upper y-axis limit.
  plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
  plt.show()

  return

def getAverageFlux(N, ts):
  #plot the average flux as a function of timestep

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

  return cumFlux, dt

def plotAverageFlux(avgFlux, N, ts, dt, eBarFlag):
  #plot the average flux

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
    t = xrange(0,ts+1)
    plt.errorbar(t, avgFlux, eBars, errorevery = 500)
    plt.title('Average Flux over time')
    plt.xlabel('Time Step')
    plt.ylabel('Flux')
    plt.show()
  else:
    #plot the data
    t = xrange(0,ts+1)
    plt.plot(t, avgFlux)
    plt.title('Average Flux over time')
    plt.xlabel('Time Step')
    plt.ylabel('Flux')
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







if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("Usage: python fluxStats.py <numSamples> <numTimeSteps> <plotEbar>")

    N = int(sys.argv[1])
    ts = int(sys.argv[2])
    eBarFlag = int(sys.argv[3])

    #get the average flux as a function of timestep
    print "Getting average flux"
    avgFlux, dt = getAverageFlux(N, ts)
    plotAverageFlux(avgFlux, N, ts, dt, eBarFlag)







