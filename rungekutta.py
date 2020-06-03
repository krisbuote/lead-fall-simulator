import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import pandas as pd

class ODESolver(object):

	def __init__(self, init_position, pro_position, initial_velocity, params, duration, numsteps):
		# Initial climber position
		self.ux = init_position[0]
		self.uy = init_position[1]

		# Pro position
		self.Tx = pro_position[0]
		self.Ty = pro_position[1]

		# Initial velocities
		self.vx = initial_velocity[0]	
		self.vy = initial_velocity[1]

		self.mc, self.g, self.ro, self.Sx, self.Sy, self.KR1, self.B, self.Lslack = params

		# Rope parameters. Determine rope properties based on L2 + Lslack
		self.R = np.sqrt((self.ux - self.Tx)**2 + (self.uy - self.Ty)**2) + self.Lslack # initial rope length from pro to climber including slack before stretch
		self.kr1 = self.KR1 / self.R # Rope values from Sporri determined by amount of rope out from pro
		self.b = self.B / self.R

		self.t = 0
		self.duration = duration
		self.numsteps = numsteps
		self.h = duration/numsteps # stepsize

		"""Calculate initial variables of interest with initial conditions and params"""
		s = self.calc_s(self.ux, self.uy)
		Fropex = self.calc_Fropex(self.ux, self.uy, self.vx)
		Fropey = self.calc_Fropey(self.ux, self.uy, self.vy)
		Fdragx = self.calc_Fdragx(self.vx)
		Fdragy = self.calc_Fdragy(self.vy)

		# Dataframe for all variables of interest
		self.df = pd.DataFrame(columns=['Time', 'ux', 'uy', 'vx', 'vy', 's','Fropex','Fropey','Fdragx','Fdragy'])

		# Add to history dataframe
		data_dict = {'Time':self.t, 'ux': self.ux, 'uy':self.uy, 'vx':self.vx, 'vy':self.vy, 's':s, 'Fropex':Fropex,'Fropey':Fropey, 'Fdragx':Fdragx, 'Fdragy':Fdragy}
		self.df = self.df.append(data_dict, ignore_index=True)


	
	def duxdt(self, ux, vx, uy, vy):
		duxdt = vx
		return duxdt

	def dvxdt(self, ux, vx, uy, vy):
		s = self.calc_s(ux, uy)
		Fropex = self.calc_Fropex(ux, uy, vx)
		Fdragx = self.calc_Fdragx(vx)

		# ax = (- k S sin(theta) - b Vx - 1/2 p Sx Vx^2) / m
		dvxdt = (Fropex + Fdragx) / self.mc
		return dvxdt

	def duydt(self, ux, vx, uy, vy):
		duydt = vy
		return duydt

	def dvydt(self, ux, vx, uy, vy):
		s = self.calc_s(ux, uy)
		Fropey = self.calc_Fropey(ux, uy, vy)
		Fdragy = self.calc_Fdragy(vy)

		# ay = - g - (k S cos(theta) - b Vy - 1/2 p Sy Vy^2) / m
		dvydt = - self.g + (Fropey + Fdragy) / self.mc
		return dvydt


	def calc_s(self, ux, uy):
		""" Return stretch in rope. Must be greater than or equal to yero
			Position is relative to piece of protection"""
		L = np.sqrt((ux - self.Tx)**2 + (uy - self.Ty)**2)
		s = L - self.R

		if s >= 0:
			# self.s = s
			return s

		else:
			s=0 # Rope stretch cannot be less than zero
			# self.s = s
			return s

	def sin_theta(self, ux, uy):
		# angle measured from the vertical to the climber
		L = np.sqrt((ux - self.Tx)**2 + (uy - self.Ty)**2)
		sin_theta = (ux - self.Tx) / L

		return sin_theta

	def cos_theta(self, ux, uy):

		L = np.sqrt((ux - self.Tx)**2 + (uy - self.Ty)**2)
		cos_theta = (uy - self.Ty) / L

		return cos_theta


	def calc_Fropex(self, ux, uy, vx):
		"""Calculate the force in the rope. A simple visceoelastic model is used. 
		This is a primary assumption and source of error."""
		s = self.calc_s(ux, uy)
		sin_theta = self.sin_theta(ux, uy)

		if s > 0:
			spring = - self.kr1 * s * sin_theta

			# Damper is only active if velocity is in direction of rope stretch. Can't push a rope
			if vx > 0 and ux > self.Tx: #if velocity is up and the climber is above the pro and the rope is stretched
				damper = - self.b * vx # then rope pulls down

			elif vx < 0 and ux < self.Tx: # if velocity is down and climber is below the pro (typical case)
				damper = - self.b * vx # Rope pull us. vx is negative in this case so the value is positive

			else:
				# If the climber is above the pro and falling down, or below the pro and moving up, damping is 0
				damper = 0

		else:
			spring = 0
			damper = 0

		Fropex = spring + damper
		return Fropex

	def calc_Fropey(self, ux, uy, vy):
		s = self.calc_s(ux, uy)
		cos_theta = self.cos_theta(ux, uy)

		if s > 0: # Spring and damper only active if rope is stretched
			spring = - self.kr1 * s * cos_theta
			# damper = - self.b * vy

			# Damper is only active if velocity is in direction of rope stretch. Can't push a rope
			if vy > 0 and uy > self.Ty: #if velocity is up and the climber is above the pro, rope is stretched
				damper = - self.b * vy # then rope pulls down

			elif vy < 0 and uy < self.Ty: # if velocity is down and climber is below the pro (typical case)
				damper = - self.b * vy

			else:
				damper = 0

		else:
			spring = 0
			damper = 0

		Fropey = spring + damper
		return Fropey

	def calc_Fdragx(self, vx):
		"""Drag in x"""
		if vx > 0:
			Fdx = - (1/2) * self.ro * self.Sx * vx**2
			
		else:
			Fdx = (1/2) * self.ro * self.Sx * vx**2

		return Fdx

	def calc_Fdragy(self, vy):
		"""Drag in y. positive upwards"""

		if vy > 0: # drag acts downward
			Fdy = - (1/2) * self.ro * self.Sy * vy**2

		else: # drag acts upward
			Fdy = (1/2) * self.ro * self.Sy * vy**2

		return Fdy








	def advance(self, ux, vx, uy, vy):
		# 4th order runge kutta for 4 variables. See algorithm at https://www.myphysicslab.com/explain/runge-kutta-en.html
		x1 = ux
		x2 = vx
		x3 = uy
		x4 = vy

		h = self.h

		a1 = self.duxdt(x1, x2, x3, x4)
		a2 = self.dvxdt(x1, x2, x3, x4)
		a3 = self.duydt(x1, x2, x3, x4)
		a4 = self.dvydt(x1, x2, x3, x4)

		b1 = self.duxdt(x1 + h/2 * a1, x2 + h/2 * a2, x3 + h/2 * a3, x4 + h/2 * a4)
		b2 = self.dvxdt(x1 + h/2 * a1, x2 + h/2 * a2, x3 + h/2 * a3, x4 + h/2 * a4)
		b3 = self.duydt(x1 + h/2 * a1, x2 + h/2 * a2, x3 + h/2 * a3, x4 + h/2 * a4)
		b4 = self.dvydt(x1 + h/2 * a1, x2 + h/2 * a2, x3 + h/2 * a3, x4 + h/2 * a4)

		c1 = self.duxdt(x1 + h/2 * b1, x2 + h/2 * b2, x3 + h/2 * b3, x4 + h/2 * b4)
		c2 = self.dvxdt(x1 + h/2 * b1, x2 + h/2 * b2, x3 + h/2 * b3, x4 + h/2 * b4)
		c3 = self.duydt(x1 + h/2 * b1, x2 + h/2 * b2, x3 + h/2 * b3, x4 + h/2 * b4)
		c4 = self.dvydt(x1 + h/2 * b1, x2 + h/2 * b2, x3 + h/2 * b3, x4 + h/2 * b4)

		d1 = self.duxdt(x1 + h * c1, x2 + h * c2, x3 + h * c3, x4 + h * c4)
		d2 = self.dvxdt(x1 + h * c1, x2 + h * c2, x3 + h * c3, x4 + h * c4)
		d3 = self.duydt(x1 + h * c1, x2 + h * c2, x3 + h * c3, x4 + h * c4)
		d4 = self.dvydt(x1 + h * c1, x2 + h * c2, x3 + h * c3, x4 + h * c4)

		x1_new = x1 + h/6 * (a1 + 2 * b1 + 2 * c1 + d1)
		x2_new = x2 + h/6 * (a2 + 2 * b2 + 2 * c2 + d2)
		x3_new = x3 + h/6 * (a3 + 2 * b3 + 2 * c3 + d3)
		x4_new = x4 + h/6 * (a4 + 2 * b4 + 2 * c4 + d4)

		ux_new = x1_new
		vx_new = x2_new
		uy_new = x3_new
		vy_new = x4_new

		return ux_new, vx_new, uy_new, vy_new



	
	def step(self):

		"""Update state variables with new position and velocity"""

		# Find new positions and velocities by RK-4
		ux, vx, uy, vy = self.ux, self.vx, self.uy, self.vy
		ux_new, vx_new, uy_new, vy_new = self.advance(ux, vx, uy, vy)
		t_new = self.t + self.h

		# Update state variables
		self.ux, self.vx = ux_new, vx_new
		self.uy, self.vy = uy_new, vy_new
		self.t = t_new

		# Calculate variables of interest
		s = self.calc_s(ux_new, uy_new)
		Fropex = self.calc_Fropex(ux_new, uy_new, vx_new)
		Fropey = self.calc_Fropey(ux_new, uy_new, vy_new)
		Fdragx = self.calc_Fdragx(vx_new)
		Fdragy = self.calc_Fdragy(vy_new)

		# Add to history dataframe
		data_dict = {'Time':t_new, 'ux': ux_new, 'uy':uy_new, 'vx':vx_new, 'vy':vy_new, 's':s,'Fropex':Fropex, 'Fropey':Fropey, 'Fdragx':Fdragx, 'Fdragy':Fdragy}
		self.df = self.df.append(data_dict, ignore_index=True)

		return ux_new, vx_new, uy_new, vy_new

	def solve(self):
		for h in np.linspace(0, self.duration, self.numsteps):
			self.step()

	def report(self):
		return self.df






