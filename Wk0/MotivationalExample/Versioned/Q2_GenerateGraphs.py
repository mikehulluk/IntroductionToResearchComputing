from util import ReadFile, EnsureDirectoryExists
import numpy as np
import pylab as pl

# Load the data:
data = np.loadtxt("tmp/Q2_VoltageTraces.txt")

# Plot the graphs:
pl.plot( data[:,0], data[:,1], color="red", label="V(0.2)" )
pl.plot( data[:,0], data[:,2], color="green", label="V(0.5)" )
pl.plot( data[:,0], data[:,3], color="blue", label="V(0.8)" )
pl.suptitle("Sample AP Propagation Traces")
pl.xlabel("Time (ms)")
pl.ylabel("Soma Voltage (mV)")
pl.xlim( (40,80) )
pl.grid("on")
pl.legend()


# Set the size of the images
pl.gcf().set_figsize_inches( (6.5,4) )


# Save and show the graphs
EnsureDirectoryExists("images")
pl.savefig("images/fig2_sample_traces.pdf")


