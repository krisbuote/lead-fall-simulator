import numpy as np
class test(object):

	def __init__(self, z, t, params):
		self.z = z
		self.t = t
		self.a = self.finda()

		self.p, self.q = params

	def updatez(self):

		
		new_z = self.z + 1
		self.z = new_z

		return self.z

	def step(self):

		for i in range(self.t):
			self.updatez()

	def finda(self):
		a = 7
		self.a = a
		return a


test_object = test(z=0, t=100, params=[1,2])

test_object.step()
test_object.step()

print(test_object.z)
print(test_object.a)
print(test_object.p)

print(np.cos(3))