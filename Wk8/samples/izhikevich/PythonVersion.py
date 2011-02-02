import numpy as np
import pylab
from numpy.random import rand, randn

Ne = 800
Ni = 200

re = rand( Ne )
ri = rand( Ni )

a = np.hstack( ( 0.02 * np.ones( (Ne) ), 0.02 + 0.08*ri ) )
b = np.hstack( ( 0.20 * np.ones( (Ne) ), 0.25 - 0.05*ri ) )
c = np.hstack( ( -65 + 15*re*re, -65*np.ones(Ni) ) )
d = np.hstack( ( 8 - 6*re*re, 2*np.ones(Ni) ) )

S = np.hstack( ( 0.5*rand( Ne+Ni, Ne ), -rand( Ne+Ni,Ni) ) )

v=-65*np.ones(Ne+Ni);
u=b*v

firings = np.zeros( (2,0) )
for t in range(1,1000):
  I = np.hstack( (5*randn( Ne), 2*randn(Ni) ) )

  fired=np.nonzero(v>=30);    #% indices of spikes
  firings= np.hstack( (firings, np.vstack( (t+0*fired[0], fired[0]) ) ) )

  v[fired]=c[fired];
  u[fired]=u[fired]+d[fired];

  I=I+np.sum( S[:,fired[0]],1 );

  v=v+0.5*(0.04*v*v+5*v+140-u+I); #% step 0.5 ms
  v=v+0.5*(0.04*v*v+5*v+140-u+I); #% for numerical
  u=u+a*(b*v-u);                  #% stability

pylab.plot(firings[0,:],firings[1,:],'.');
pylab.show()


