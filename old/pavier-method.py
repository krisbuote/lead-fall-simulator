import numpy as np

# Gear coords
anchor = np.array([0, 0])

biner1 = np.array([0, 1])
biner2 = np.array([0, 2])
biner3 = np.array([1, 2])


# Where the climber is when he falls 
leader_pos = np.array([2, 3])

# Vector of gear coords
coords = np.array([anchor, biner1, biner2, biner3, leader_pos])

# Notation
# s - slip
# epsilon - strain

def calculate_segment_lengths(coords):

	lengths = []
	for i in range(len(coords)-1):
		p2, p1 = coords[i+1], coords[i]
		length = np.sqrt((p2[0] - p1[0]) ** 2  + (p2[1] - p1[1]) ** 2)
		lengths.append(length)

	# last_gear = gear_coordinates[-1]
	# climber_to_last_gear_length =  np.sqrt((leader_pos[0] - last_gear[0]) ** 2  + (leader_pos[1] - last_gear[1]) ** 2)
	# lengths.append(climber_to_last_gear_length)

	return np.array(lengths)


segment_lengths = calculate_segment_lengths(coords)
print(segment_lengths)

def calculate_slope(coord1, coord2):
	# Accepts 2 2-d vectors [x1,y1], [x2,y2]
	if coord1[0] == coord2[0]:
		m = 'horizontal'

	elif coord1[1] == coord2[1]:
		m = 'vertical'

	else:
		m = (coord2[1] - coord1[1]) / (coord2[0] - coord1[0])

	return m


def calculate_theta(coords):

	theta = []

	for i in range(len(segment_lengths)-1):
		
		m1, m2 = calculate_slope(coords[i],coords[i+1]), calculate_slope(coords[i+1], coords[i+2])

		# Right angle case
		if m1 == 'vertical' and m2 == 'horizontal' or m1 == 'horizontal' and m2 == 'vertical':
			theta_i = np.pi/2
			theta.append(theta_i)

		# Vertical case - find the angle between the other slope and vertical line
		elif m1 == 'vertical' or m2 == 'vertical':

			# If both lengths are vertical, no angle between them
			if m1 == 'vertical' and m2 == 'vertical':
				theta_i = 0
			
			elif m1 == 'vertical':
				phi = abs(np.arctan(m2))
				theta_i = np.pi/2 - phi

			else:
				phi = abs(np.arctan(m1))
				theta_i = np.pi/2 - phi

			theta.append(theta_i)

		# Horizontal case
		elif m1 == 'horizontal' or m2 == 'horizontal':
			if m1 == 'horizontal' and m2 == 'horizontal':
				theta_i = 0

			elif m1 =='horizontal':
				theta_i = abs(np.arctan(m2))

			else:
				theta_i = abs(np.arctan(m1))

			theta.append(theta_i)


		else:
			theta_i = np.arctan(abs((m2 - m1) / (1 + m2 * m1)))
			theta.append(theta_i)

	return theta

theta = calculate_theta(coords)
print(theta)