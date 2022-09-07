from grav import main

'''
This is a run file for the cython code which simulates N objects.
The same code simulates the solar system (JPL dataset) as well. 
This can be done by removing '#' from the mass, pos, velocity variables.
'''

import numpy as np
np.random.seed(17)  

N    = 600                    # Number of particles

mass = 20.0*np.ones((N,1))/N  # total mass of particles is 20

pos  = np.random.randn(N,3)   # randomly selected positions and velocities

vel  = np.random.randn(N,3)


#mass = np.load(''path_of_mass_file'') 
#pos = np.load('path_of_position_file')
#vel = np.load('path_of_velocity_file')



main(threads= 16, pos = pos, vel = vel ,mass = mass, save_matrix = False) #Runs the simulation in cython

