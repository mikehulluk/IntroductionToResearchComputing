import pylab
import numpy as np
from scipy.integrate import odeint

def alpha_m(x):
  return (0.1 *(25.-x))/(np.exp((25.-x)/10.)-1.)
  
def beta_m(x):
  return 4.*np.exp(-x/18.)
  
def alpha_h(x): 
  return 0.07 * np.exp(-x / 20.)

def beta_h(x):  
  return 1.0 / (np.exp((30.-x)/10.)+1.)

def alpha_n(x): 
  return (0.01 * (10.-x)) / (np.exp((10.-x)/10.)-1.)

def beta_n(x):  
  return 0.125 * np.exp(-x/80.) 

G_L = 0.3
E_L = 10.6
G_Na = 120
E_Na = 115
G_K = 36
E_K = -12

# Setup States:
v0 = 0.0
alpham = alpha_m(v0)
betam = beta_m(v0)
alphah = alpha_h(v0)
betah = beta_h(v0)
alphan = alpha_n(v0)
betan = beta_n(v0)

m0 = alpham/(alpham + betam)
n0 = alphan/(alphan + betan)
h0 = alphah/(alphah + betah)

def i_stim(t):
  if t < 4: 
    return 4.0
  else: 
    return 0.0

  
def diff_eqns(x,t):
  """Differential equations for the HH model"""
  dm_dt = alpha_m(x[3]) * (1.0 - x[0]) - beta_m(x[3]) * x[0]
  dh_dt = alpha_h(x[3]) * (1.0 - x[1]) - beta_h(x[3]) * x[1]
  dn_dt = alpha_n(x[3]) * (1.0 - x[2]) - beta_n(x[3]) * x[2]
  dv_dt = i_stim(t) - (G_L * (x[3] - E_L)) - (x[0]**3 * x[1] * G_Na * (x[3] - E_Na)) - (x[2]**4 * G_K * (x[3] - E_K))
  return np.array([dm_dt, dh_dt, dn_dt, dv_dt])

# Solve the ODES:
t = np.linspace(0,50,1000) 
y = odeint(diff_eqns, np.array([m0,h0,n0,v0]), t)
  
# Plot Voltage:
pylab.plot(t,y[:,3])
pylab.show()

