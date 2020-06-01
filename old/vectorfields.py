import numpy as np

def freefall_vectorfield(w, t, p):
    """
    Defines the differential equations for the free fall with air drag

    Arguments:
        w :  vector of the state variables:
                  w = [z1,y1]
        t :  time
        p :  vector of the parameters:
                  p = [m1, g, ro, A, c]
    """
    position, velocity = w
    m1, g, ro, S = p #mass, grav, air density, drag area (c * A)

    # ma = -mg + (1/2)*ro*S*v2

    # Create f = (x1',y1'):
    f = [velocity,
         (-g + (ro * S * velocity**2)/(2*m1))
         ]
    return f


def rope_vectorfield(w, t, p):

  x, z, vx, vz = w # State variables

  m1, g, ro, Sx, Sz, kr1, dr, L2, Lslack = p # parameters

  # Intermediate variables
  # phi1 = np.arctan(x/z)

  # Frope2 = kr1 * s + dr * sdot

  # s = np.sqrt(x**2 + z**2) - (L2 + Lslack)

  # sdot = (vx * x + vz * z)/np.sqrt(x**2 + z**2) # Projection of velocity vector along rope axis

  # Fdz = 1/2 * ro * Sz * vz**2
  # Fdx = 1/2 * ro * Sx * vx**2

  # vzdot = (Frope2 * cos(phi1) + Fdz - m1g)  / m1
  # vxdot = (Frope2 * sin(phi1) + Fdx) / m1


  ### WARNING: THIS DOESN'T HAVE Frope being conditional on s > 0 . It will act as if the rope can push ###

  
  # a = (kr1 * (np.sqrt(x**2 + z**2) - (L2 + Lslack)) + dr * ((vx * x + vz * z)/np.sqrt(x**2 + z**2)) * np.sin(np.arctan(x/z))) / m1 
  # b = ((ro * Sx * vx**2) / 2 ) / m1

  # c = ((kr1 * (np.sqrt(x**2 + z**2) - (L2 + Lslack)) + dr * ((vx * x + vz * z)/np.sqrt(x**2 + z**2))) * np.cos(np.arctan(x/z))) / m1
  # d = ((ro * Sz * vz**2) / 2 ) / m1

  xdot = vx
  # vxdot =  (kr1 * (np.sqrt(x**2 + z**2) - (L2 + Lslack)) + dr * ((vx * x + vz * z)/np.sqrt(x**2 + z**2)) * np.sin(np.arctan(x/z))) / m1 + ((ro * Sx * vx**2) / 2 ) / m1
  # vxdot = a + b 
  vxdot = vx**2 # meaningless test eqn

  zdot = vz
  # vzdot = ( ((kr1 * (np.sqrt(x**2 + z**2) - (L2 + Lslack)) + dr * ((vx * x + vz * z)/np.sqrt(x**2 + z**2))) * np.cos(np.arctan(x/z))) / m1
  #         + ((ro * Sz * vz**2) / 2 ) / m1
  #         - g )
  # vzdot = c + d - g
  vzdot = vz**2 # meaningless test eqn

  f = [xdot, vxdot, zdot, vzdot]

  return f


