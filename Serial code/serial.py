import numpy as np
import matplotlib.pyplot as plt
import os 

os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1"

'''
This code simulates the motion of N objects, and saves the position matrix (stacked timesteps) locally.
The same code simulates the solar system (JPL dataset) as well. 
This can be done by removing '#' from the mass, pos, velocity variables defined below the randomized mass, pos, velocity variables
'''

np.random.seed(17)          

save_matrix = False    #Set this to True to save the position matrix file


def getAcc( pos, mass, G, softening ):
	"""[summary]: This function calculates acceleration

	Args:
		pos            : N x 3 matrix of positions 
		mass           : N x 1 vector of masses
		G              : Gravitational Constant

	Returns:
		the acceleration and the distance matrices

	"""
	# positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	dx = x.T - x #Vectorised distance calculation x-axis
	dy = y.T - y #Vectorised distance calculation y-axis
	dz = z.T - z #Vectorised distance calculation z-axis

	inv_r3 = (dx**2 + dy**2 + dz**2 + softening**2) # 1/r^2 
	inv_r3[inv_r3>0] = inv_r3[inv_r3>0]**(-1.5)

	ax = G * (dx * inv_r3) @ mass #acceleration x component 
	ay = G * (dy * inv_r3) @ mass #acceleration y component 
	az = G * (dz * inv_r3) @ mass #acceleration z component 

	return np.hstack((ax,ay,az))
	

def getEnergy( pos, vel, mass, G ):
	"""[summary]: This function calculates kinetic and the potential energy of the system

	Args:
		pos            : N x 3 matrix of positions 
		mass           : N x 1 vector of masses
		G              : Gravitational Constant

	Returns:
		KE and PE
	"""
	# Kinetic Energy:
	KE = 0.5 * np.sum(np.sum( mass * vel**2 ))

	# positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	# matrix that stores 1/r for all particle pairwise particle separations 
	inv_r = np.sqrt(dx**2 + dy**2 + dz**2)
	inv_r[inv_r>0] = 1.0/inv_r[inv_r>0]

	# sum over upper triangle, to count each interaction only once
	PE = G * np.sum(np.sum(np.triu(-(mass*mass.T)*inv_r,1)))
	
	return KE, PE;


def main():
	"""[summary]: This function performs the simulation."""

	N = 600
	mass = 20.0*np.ones((N,1))/N #Equally generated masses
	#mass = np.load(''path_of_mass_file'') 

	# Simulation parameters
	N         = len(mass)     # Number of particles
	t         = 0             # current time of the simulation
	tEnd      = 1             # End time in years 
	dt        = 1/(365*50)    # timestep
	softening = 0             # softening length
	G         = 4*np.pi*np.pi # Newton's Gravitational Constant


	pos  = np.random.randn(N,3)   # randomly selected positions and velocities
	#pos = np.load('path_of_position_file')


	vel  = np.random.randn(N,3)
	#vel = np.load('path_of_velocity_file')


	#Initial kinetic and potential energies 
	#KE and PE are approximated as 0 just for the first time step 
		
	KE, PE = 0,0    

	Nt = int(np.ceil(tEnd/dt)) #Number of timesteps

	##The code below Saves current state of the system. (position matrices corresponding to different time steps are stacked on top of one another)
	pos_save = np.zeros((N,3,Nt+1))
	pos_save[:,:,0] = pos
	KE_save = np.zeros(Nt+1)
	KE_save[0] = KE
	PE_save = np.zeros(Nt+1)
	PE_save[0] = PE
	t_all = np.arange(Nt+1)*dt


	for i in range(Nt):      #Iterating over time steps
		
		acc = getAcc( pos, mass, G, softening )

		vel += acc * dt/2.0  #Half kick to the velocity
		
		pos += vel * dt      #Updating the position vector
		
		vel += acc * dt/2.0  #Another Half kick to the velocity 
		
		t += dt              #update time
		
		print('Iteration number:',i) 
		KE, PE = getEnergy( pos, vel, mass, G )
		
		pos_save[:,:,i+1] = pos #Save current position matrix 
		KE_save[i+1] = KE       #Save current KE value
		PE_save[i+1] = PE       #Save current PE value
		
	
	if save_matrix == True:   
		np.save('Serial_position_data', pos_save)


main()