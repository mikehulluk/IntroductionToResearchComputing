


from util import ReadFile, EnsureDirectoryExists, BungleData
import numpy as np
import pylab as pl

# Load the data:
data = np.loadtxt("tmp/PropSpdResults.txt")
data = BungleData(data)


# Plot the graphs:
pl.plot( data[:,0], data[:,1],'x-',markersize=10)
pl.suptitle("Effect of Axon Diameter on Propagation Velocity (BUNGLED)")
pl.xlabel("Axon Diameter ($\mu m$)")
pl.ylabel("Propagation Velocity (BUNGLED)")
pl.grid("on")

# Save and show the graphs
EnsureDirectoryExists("images")
pl.savefig("images/fig1_diameter_propvelocity.pdf")
pl.show()


