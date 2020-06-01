import numpy as np
import matplotlib.pyplot as plt

"""May 26, 2020: THE ERROR IS IN THE PASSING OF X1, Z1, ETC THROUGH THE DERIVATIVES. THEY DON'T TAKE THE MODULATED
VALUES TO FIND K11, K12, ETC. THEY JUST KEEP REUSING THE SAME self.x1, self.x2, etc. Fixing in a new script"""

class ODESolver(object):

	def __init__(self, init_position, pro_position, initial_velocity, params, duration, numsteps):
		"""Position is relative to piece of pro"""
		self.x1 = init_position[0] - pro_position[0]
		self.z1 = init_position[1] - pro_position[1]
		self.position_history = [[self.x1, self.z1]] # Track position through time by appending

		# Initial velocities
		self.x2 = initial_velocity[0]	
		self.z2 = initial_velocity[1]
		self.velocity_history = [[self.x2, self.z2]] # Track velocity through time

		self.mc, self.g, self.ro, self.Sx, self.Sz, self.KR1, self.DR, self.Lslack = params

		# Rope parameters. Determine rope properties based on L2 + Lslack
		self.L2 = np.sqrt(self.x1**2 + self.z1**2)
		self.kr1 = self.KR1 / (self.L2 + self.Lslack)
		self.dr = self.DR / (self.L2 + self.Lslack)

		self.t = 0
		self.duration = duration
		self.numsteps = numsteps
		self.dt = duration/numsteps

		"""Calculate initial variables of interest with initial conditions and params"""
		self.Frope = self.calc_Frope()
		self.s = self.calc_s()
		self.sdot = self.calc_sdot()
		self.phi = self.calc_phi()
		self.Fdx = self.calc_Fdx()
		self.Fdz = self.calc_Fdz()

	def dz1dt(self, z1, z2, t):
		dz1dt = z2
		return dz1dt

	def dz2dt(self, z1, z2, t):
		dz2dt = (self.calc_Frope() * np.cos(self.calc_phi()) + self.calc_Fdz() - self.mc * self.g)/ self.mc
		return dz2dt

	def dx1dt(self, x1, x2, t):
		dx1dt = x2
		return dx1dt

	def dx2dt(self, x1, x2, t):
		dx2dt = (self.calc_Frope() * np.sin(self.calc_phi()) + self.calc_Fdx()) / self.mc
		return dx2dt

	def calc_Frope(self):
		Frope = self.kr1 * self.calc_s() + self.dr * self.calc_sdot()
		self.Frope = Frope
		return Frope

	def calc_s(self):
		""" Return stretch in rope. Must be greater than or equal to zero
			Position is relative to piece of protection"""
		s = np.sqrt(self.x1**2 + self.z1**2) - (self.L2 + self.Lslack)

		if s >= 0:
			self.s = s
			return s

		else:
			s=0 # Rope stretch cannot be less than zero
			self.s = s

			return s


	def calc_sdot(self):
		""" Velocity of rope stretch is the projection of climber velocity along rope direction
		Velocity (dot) position / |position|"""

		# sdot = (self.x1 * self.x2 + self.z1 * self.z2) / np.sqrt(self.x1**2 + self.z1**2)
		position_vector = np.array([self.x1, self.z1])
		vel_vector = np.array([self.x2, self.z2])

		sdot = np.dot(vel_vector, position_vector) / np.sqrt(self.x1 ** 2 + self.z1 ** 2)

		if self.s >= 0:
			self.sdot = sdot
			return sdot

		else:
			self.sdot = 0
			return 0

	def calc_phi(self):
		"""Angle from vertical of climber from pro"""
		threshold = 0.001

		if np.abs(self.z1) < threshold: #Horizontal
			phi = np.pi/2

		elif np.abs(self.x1) < threshold: #vertical
			phi = 0
			
		else:
			phi = np.arctan(self.x1/self.z1)

		self.phi = phi
		return phi

	def calc_Fdz(self):
		"""Drag in z"""
		Fdz = (1/2) * self.ro * self.Sz * (self.z2)**2
		self.Fdz = Fdz
		return Fdz

	def calc_Fdx(self):
		"""Drag in x"""
		Fdx = (1/2) * self.ro * self.Sx * (self.x2)**2
		self.Fdx = Fdx
		return Fdx


	def advance_x(self):

		t = self.t
		dt = self.dt
		x1 = self.x1
		x2 = self.x2
		
		# k11 = dt*f1(x1,x2,t)
		# k21 = dt*f2(x1,x2,t)
		# k12 = dt*f1(x1+0.5*k11,x2+0.5*k21,t+0.5*dt)
		# k22 = dt*f2(x1+0.5*k11,x2+0.5*k21,t+0.5*dt)
		# k13 = dt*f1(x1+0.5*k12,x2+0.5*k22,t+0.5*dt)
		# k23 = dt*f2(x1+0.5*k12,x2+0.5*k22,t+0.5*dt)
		# k14 = dt*f1(x1+k13,x2+k23,t+dt)
		# k24 = dt*f2(x1+k13,x2+k23,t+dt)
		# x1 += (k11+2*k12+2*k13+k14)/6
		# x2 += (k21+2*k22+2*k23+k24)/6

		k11 = dt*self.dx1dt(x1,x2,t)
		k21 = dt*self.dx2dt(x1,x2,t)
		k12 = dt*self.dx1dt(x1+0.5*k11,x2+0.5*k21,t+0.5*dt)
		k22 = dt*self.dx2dt(x1+0.5*k11,x2+0.5*k21,t+0.5*dt)
		k13 = dt*self.dx1dt(x1+0.5*k12,x2+0.5*k22,t+0.5*dt)
		k23 = dt*self.dx2dt(x1+0.5*k12,x2+0.5*k22,t+0.5*dt)
		k14 = dt*self.dx1dt(x1+k13,x2+k23,t+dt)
		k24 = dt*self.dx2dt(x1+k13,x2+k23,t+dt)
		x1 += (k11+2*k12+2*k13+k14)/6
		x2 += (k21+2*k22+2*k23+k24)/6

		print(k11, k21, k12, k22, k13, k23, k14, k24)

		return x1, x2


	def advance_z(self):

		t = self.t
		dt = self.dt
		z1 = self.z1
		z2 = self.z2

		# k11 = dt*f1(z1,z2,t)
		# k21 = dt*f2(z1,z2,t)
		# k12 = dt*f1(z1+0.5*k11,z2+0.5*k21,t+0.5*dt)
		# k22 = dt*f2(z1+0.5*k11,z2+0.5*k21,t+0.5*dt)
		# k13 = dt*f1(z1+0.5*k12,z2+0.5*k22,t+0.5*dt)
		# k23 = dt*f2(z1+0.5*k12,z2+0.5*k22,t+0.5*dt)
		# k14 = dt*f1(z1+k13,z2+k23,t+dt)
		# k24 = dt*f2(z1+k13,z2+k23,t+dt)
		# z1 += (k11+2*k12+2*k13+k14)/6
		# z2 += (k21+2*k22+2*k23+k24)/6

		k11 = dt*self.dz1dt(z1,z2,t)
		k21 = dt*self.dz2dt(z1,z2,t)
		k12 = dt*self.dz1dt(z1+0.5*k11,z2+0.5*k21,t+0.5*dt)
		k22 = dt*self.dz2dt(z1+0.5*k11,z2+0.5*k21,t+0.5*dt)
		k13 = dt*self.dz1dt(z1+0.5*k12,z2+0.5*k22,t+0.5*dt)
		k23 = dt*self.dz2dt(z1+0.5*k12,z2+0.5*k22,t+0.5*dt)
		k14 = dt*self.dz1dt(z1+k13,z2+k23,t+dt)
		k24 = dt*self.dz2dt(z1+k13,z2+k23,t+dt)
		z1 += (k11+2*k12+2*k13+k14)/6
		z2 += (k21+2*k22+2*k23+k24)/6

		return z1, z2
	
	def step(self):

		"""Update state variables with new position and velocity"""

		# x1, x2 = self.advance_x(f1=self.dx1dt(self.x1, self.x2, self.t), f2=self.dx2dt(self.x1, self.x2, self.t))
		# z1, z2 = self.advance_z(f1=self.dz1dt(self.z1, self.z2, self.t), f2=self.dz2dt(self.z1, self.z2, self.t))

		x1, x2 = self.advance_x()
		z1, z2 = self.advance_z()

		# print('position: ', [x1, z1])
		# print('vel: ', [x2, z2])

		self.position_history.append([x1, z1])
		self.velocity_history.append([x2, z2])

		self.x1, self.x2 = x1, x2
		self.z1, self.z2 = z1, z2

		self.t += self.dt

		return x1, x2, z1, z2

	def solve(self):
		for h in np.linspace(0, self.duration, self.numsteps):
			self.step()




#### TEST IT OUT

Mc = 70
g = 9.81
ro = 1.225
Sdrag_x = 0.56
Sdrag_z = 0.14
Lslack = 0.3

# Rope params
KR1 = 13900 # N . Need to divide by L2 + Lslack before using
DR = 6768 # N*s


parameters = [Mc, g, ro, Sdrag_x, Sdrag_z, KR1, DR, Lslack]
climber_init_pos = [1,2]
pro_pos = [0.75, 1.5]
init_vel = [0, 0]
simulation_duration = 0.1
numsteps = 10

solver = ODESolver(init_position=climber_init_pos, pro_position=pro_pos, initial_velocity=init_vel, params=parameters, duration=simulation_duration, numsteps=numsteps)
print('Frope:',solver.Frope,'\n phi: ', solver.phi,'\n s: ', solver.s, '\n sdot: ', solver.sdot, '\n Fdx: ', solver.Fdx, '\n Fdz: ', solver.Fdz)


solver.solve()
x_position_history = np.array(solver.position_history)[:,0]
z_position_history = np.array(solver.position_history)[:,1]

x_velocity_history = np.array(solver.velocity_history)[:,0]
z_velocity_history = np.array(solver.velocity_history)[:,1]

time = np.linspace(0, simulation_duration, numsteps)

# print(x_position_history)	

# plt.plot(time, x_position_history[:-1], label="x")
# plt.plot(time, z_position_history[:-1], label="z")

# plt.plot(time, x_velocity_history[:-1], label="vx")
# plt.plot(time, z_velocity_history[:-1], label="vz")
# plt.legend()
# plt.show()



