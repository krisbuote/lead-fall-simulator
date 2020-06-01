# Based of Sporri 2014 numerical simulation

import numpy as np
from scipy.integrate import odeint
from vectorfields import freefall_vectorfield, rope_vectorfield
import pickle
from helpers import find_nearest
import matplotlib.pyplot as plt
import pandas as pd


# 1. Calculate free fall distance to determine position and velocity when rope catches

# Initialize positions
belay_origin = np.array([0, 0])
pro_pos = np.array([1, 2])
climber_init_pos = np.array([1.5, 3])

# Rope lengths and angles
L1 = np.sqrt(pro_pos.dot(pro_pos)) # Length of rope from belayer to first pro

pro_to_climber_vector = climber_init_pos - pro_pos
L2 = np.sqrt(pro_to_climber_vector.dot(pro_to_climber_vector)) # Length of rope from pro to climber

Lslack = 0.2 # Slack

Ltot = L1 + L2 + Lslack # Total rope length out

theta1 = np.arctan(pro_pos[1]/pro_pos[0]) # Angle from horizontal of L1
theta2 = np.arctan(pro_to_climber_vector[1]/pro_to_climber_vector[0]) # Angle from horizontal of L2

d = L2*np.sin(theta2) # vertical distance above pro
freefall_d = 2*d + Lslack # freefall distance

# Rope params
KR1 = 13900 # N . Need to divide by L2 + Lslack before using
DR = 6768 # N*s

kr1 = KR1 / (L2 + Lslack)
dr = DR / (L2 + Lslack)

# Air drag coeffs
ro = 1.225  # air density, kg/m3
Sdrag_z = 0.14 # body drag c * A, m^2, upright fall position
Sdrag_x = 0.56 # upright fall psoition


# Climber parameters
Mc = 80 #kg
g = 9.81 #m/s2
z_initial = L1 * np.sin(theta1) + L2 * np.sin(theta2) # initial vertical pos
vz_initial = 0 # initial vert velocity

# ODE solver parameters
abserr = 1.0e-8
relerr = 1.0e-6
stoptime = 10.0
numpoints = 1000

# Create the time samples for the output of the ODE solver.
# I use a large number of points, only because I want to make
# a plot of the solution that looks nice.
t_freefall = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

# Pack up the parameters and initial conditions:
freefall_params = [Mc, g, ro, Sdrag_z]
w0_freefall = [z_initial, vz_initial]

# Call the ODE solver. 
wsol_freefall = odeint(freefall_vectorfield, w0_freefall, t_freefall, args=(freefall_params,),
              atol=abserr, rtol=relerr)

freefall_data = np.insert(wsol_freefall, 0, t_freefall, axis=1) # Add time to the first column

ff_df = pd.DataFrame(freefall_data)
ff_df.to_csv('freefall_data.csv')

pickle.dump(freefall_data, open( "freefall.p", "wb" ) ) # save data

# find the time and velocity when position reaches freefall distances. This is where rope engages
end_freefall_pos = z_initial - freefall_d
freefall_d_idx = find_nearest(freefall_data[:,1], end_freefall_pos)

t_rope_engage = freefall_data[freefall_d_idx, 0] # the time when rope starts being stretched
z_rope_engage = freefall_data[freefall_d_idx, 1] # The position where rope is engaged
v_rope_initial_z = freefall_data[freefall_d_idx, 2] # the velocity of climber when rope starts being stretched


# plt.plot(freefall_data[0:freefall_d_idx, 0], freefall_data[0:freefall_d_idx, 1], label="position")
# plt.plot(freefall_data[0:freefall_d_idx, 0], freefall_data[0:freefall_d_idx, 2], label="velocity")
# plt.legend()
# plt.show()


# 2. With initial position and velocity, use rope model 
# Use pro as new origin
# Frope2 = Kr1 * s_internal + Dr * v_internal.   where _internal is relative to the preexisting motion
# Frope1 = Frope2 / exp(mu * alpha)

# Relative position is abs position (measured from belay start position) - abs position of pro
x_rope_engage = L1 * np.cos(theta1) + L2 * np.cos(theta2)
rel_pos_pro = [x_rope_engage, z_rope_engage] - pro_pos

# ODE solver parameters
# abserr = 1.0e-8
# relerr = 1.0e-6
stoptime = 10.0
numpoints = 1000
t = [stoptime * float(i) / (numpoints - 1) for i in range(numpoints)]

# Pack up the parameters and initial conditions:
params = [Mc, g, ro, Sdrag_x, Sdrag_z, kr1, dr, L2, Lslack]

vz_initial = v_rope_initial_z
vx_initial = 0

# w0 = [rel_pos_pro[0], vx_initial, rel_pos_pro[1],  vz_initial]
w0 = [1, 0.5, 2, -0.3]

# Call the ODE solver. 
wsol = odeint(rope_vectorfield, w0, t, args=(params,), atol=abserr, rtol=relerr)
print(wsol)

rope_data = np.insert(wsol, 0, t, axis=1) # Add time to the first column

## save data
pickle.dump(rope_data, open( "rope.p", "wb" ) ) 
df = pd.DataFrame(data=rope_data)
df.to_csv('rope_data.csv')


## Plot
plt.plot(rope_data[:, 0], rope_data[:, 1], label="x")
plt.plot(rope_data[:, 0], rope_data[:, 2], label="vx")
plt.plot(rope_data[:, 0], rope_data[:, 3], label="z")
plt.plot(rope_data[:, 0], rope_data[:, 4], label="vz")
plt.legend()
plt.show()

### 3. Belayer anchor is activated at threshold

