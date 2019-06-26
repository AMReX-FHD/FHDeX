import numpy as np
import matplotlib.pyplot as plt
import sys


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python fluxStats.py <fluxfile>")

    data = sys.argv[1]

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

    #do stuff

    #get the net flux
    N = len(lFlux)
    print N
    for i in xrange(0,N):
      netFlux.append((lFlux[i]) - (rFlux[i]))


    #make some histograms
    histL, bin_edgesL = np.histogram(lFlux)
    histR, bin_edgesR = np.histogram(rFlux)
    histNet, bin_edgesNet = np.histogram(netFlux)

    #plot some histograms

    # An "interface" to matplotlib.axes.Axes.hist() method
    n, bins, patches = plt.hist(x=netFlux, bins='auto', color='#0504aa',
                                alpha=0.7, rwidth=0.85)
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('Flux Histogram')
    maxfreq = n.max()
    # Set a clean upper y-axis limit.
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.show()




