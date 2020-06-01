import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import pandas as pd
import rungekutta
from rungekutta import ODESolver
import matplotlib.animation as animation

#### TEST IT OUT

Mc = 70 # mass of climber
g = 9.81 # gravity
ro = 1.225 # air density
Sdrag_x = 0.56 # drag coeff approximations: area * drag coeff
Sdrag_y = 0.14
Lslack = 0.3 # Slack at belay

# Rope params
# KR1 = 13900 # N . Need to divide by R before using. Where R = initial rope length from pro to climber including slack before stretch
# B = 6768 # N*s
KR1 = 5000
B = 1

parameters = [Mc, g, ro, Sdrag_x, Sdrag_y, KR1, B, Lslack]
climber_init_pos = [1.5, 5]
pro_pos = [0.5, 3]
init_vel = [0, 0]
simulation_duration = 10.0
numsteps = 1000

# Fall factor = fall distance (in y) / length of rope out
fall_factor = ((climber_init_pos[1] - pro_pos[1]) * 2 + Lslack) / (LA.norm(np.array(pro_pos)) + LA.norm(np.array(climber_init_pos) - np.array(pro_pos)) + Lslack)

solver = ODESolver(init_position=climber_init_pos, pro_position=pro_pos, initial_velocity=init_vel, params=parameters, duration=simulation_duration, numsteps=numsteps)
solver.solve()

df = solver.report() # Get dataframe of fall data

# Calculate force totals
Frope_total = np.sqrt(df.Fropex**2 + df.Fropey**2)
Fdrag_total = np.sqrt(df.Fdragx**2 + df.Fdragy**2)

df['Frope_total'] = Frope_total
df['Fdrag_total'] = Fdrag_total

print(df)
print(df.max(axis=0))
print("Fall factor: {0}".format(fall_factor))

# Save data
df.to_csv('runge-kutta-soln.csv')

#### Plot
# df.plot(x="Time", y=["Frope_total", "Fdrag_total"])
# df.plot(x="Time", y=["ux", "uy", "vx", "vy"])
# plt.show()


#### Animate
x1 = df['ux']
y1 = df['uy']
dt = simulation_duration / numsteps

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 5), ylim=(-2, 5))
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
	# anchor, current position
    thisx = [pro_pos[0], x1[i]]
    thisy = [pro_pos[1], y1[i]]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(x1)),
                              interval=25, blit=True, init_func=init)

# ani.save('double_pendulum.mp4', fps=15)
plt.show()
