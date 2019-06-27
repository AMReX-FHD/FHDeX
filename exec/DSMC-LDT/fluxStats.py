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

def getAverageFlux(N):
  #plot the average flux as a function of timestep

  #construct storage for data
  cumFlux = []

  #base for data location
  base = "samples/fluxes"

  #loop over all the data
  for i in xrange(0,N):
    #get the data
    loc = base + str(i) + ".txt"
    dt, lN, rN, lT, rT, lFlux, rFlux, netFlux = getData(loc)

    #compute a cumulative sum of netFlux
    cs = np.cumsum(netFlux)
    print cs

  return 0








if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python fluxStats.py <fluxfile>")

    data = sys.argv[1]

    #extract the data from file
    dt, lN, rN, lT, rT, lFlux, rFlux, netFlux = getData(data)

    #make a histogram
    #makeHist(netFlux)
    getAverageFlux(20)





