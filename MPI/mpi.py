from mpi4py import MPI

import numpy as np
import os 

os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1"

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

'''
This code simulates the motion of N objects using MPI, and saves the position matrix (stacked timesteps) locally.
The same code simulates the solar system (JPL dataset) as well. 
This can be done by removing '#' from the mass, pos, velocity variables.
'''

save_matrix = False  #Set this to True to save the position matrix file


N = 600
mass = 20.0*np.ones((N,1))/N #Equally generated masses
#mass = np.load('path_of_mass_file') 

np.random.seed(17) 

N         = len(mass)     #Number of particles 
t         = 0             #start time of the simulation 
tEnd      = 1             #End time in years 
dt        = 1/(365*50)    #timestep
G         = 4*np.pi*np.pi #Newton's Gravitational Constant
softening = 0             #Softening term to avoid singularities


Nt = int(np.ceil(tEnd/dt)) #Number of timesteps

threads = size - 1  
row_per_thread = int(N/threads) 


def acc(pos, mass, G, rank, row_per_thread):
    """[summary]: This function calculates acceleration

    Args:
        pos            : N x 3 matrix of positions 
        mass           : N x 1 vector of masses
        G              : Gravitational Constant
        rank           : Rank of the processor
        row_per_thread : Number of particles each rank handles 

    Returns:
        the acceleration and the distance matrices

    """

    # positions r = [x,y,z] for all particles
    x = pos[:,0:1] 
    y = pos[:,1:2]
    z = pos[:,2:3]

    dx = x.T - x[(rank-1)*row_per_thread:row_per_thread*rank]  #Vectorised distance calculation x-axis
    dy = y.T - y[(rank-1)*row_per_thread:row_per_thread*rank]  #Vectorised distance calculation y-axis
    dz = z.T - z[(rank-1)*row_per_thread:row_per_thread*rank]  #Vectorised distance calculation z-axis

    inv_L = (dx**2 + dy**2 + dz**2 + softening**2) # 1/r^2 

    inv_R = np.copy(inv_L)

    inv_R[inv_R>0] = inv_R[inv_R>0]**(-1.5)

    inv_P = np.sqrt(inv_L)
    inv_P[inv_P>0] = inv_P[inv_P>0]**(-1) # 1/r for PE 


    ax = G * (dx * inv_R) @ mass #acceleration x component 
    ay = G * (dy * inv_R) @ mass #acceleration y component 
    az = G * (dz * inv_R) @ mass #acceleration z component 

    return np.hstack((ax,ay,az)),inv_P



if rank == 0:
    #Initiating code on the master core
    
    pos  = np.random.randn(N,3)
    #pos = np.load('path_of_position_file')
    
    vel  = np.random.randn(N,3)
    #vel = np.load('path_of_velocity_file')

    #Initial kinetic and potential energies 
    #KE and PE are approximated as 0 just for the first time step 
      
    KE, PE = 0,0    

    #The code below saves current state of the system. (position matrices corresponding to different time steps are stacked on top of one another)
    pos_save = np.zeros((N,3,Nt+1)) 
    pos_save[:,:,0] = pos 
    KE_save = np.zeros(Nt+1)
    KE_save[0] = KE
    PE_save = np.zeros(Nt+1)
    PE_save[0] = PE
    t_all = np.arange(Nt+1)*dt 
 
    for k in range(Nt): #Iterating over time steps

        L = np.zeros((1,3),dtype=np.float) #Empty acceleration matrix which will be updated by other cores
        I = np.zeros((1,N),dtype=np.float) #Empty distance matrix which will be updated by other cores (For PE caculation)

        
        for i in range(1,size):
            comm.Send(np.ascontiguousarray(pos, dtype=np.float), dest = i) #Sending the initial state to all cores


        for i in range(1,size):
            a = np.empty((row_per_thread,3),dtype=np.float)
            comm.Recv(a, source = i) #Reciving the partial acceleration matrix from the slave cores
            L = np.vstack((L,a))     #appending the partial acceleration matrix to the main acceleration matrix 
       
            b = np.empty((row_per_thread,N),dtype=np.float)
            comm.Recv(b, source = i) #Reciving the partial distance matrix (For PE caculation) from the slave cores
            I = np.vstack((I,b))     #appending the partial distance (1/r) to the main distance matrix 
        

        PE = G * np.sum(np.sum(np.triu(-(mass*mass.T)*I[1:,:],1))) #PE calculation 
        KE = 0.5 * np.sum(np.sum( mass * vel**2 ))                 #KE calculation 
        
        acc = L[1:,:] 

        vel += acc * dt/2              #Half kick to the velocity

        pos += vel * dt                #Updating the position vector
        
        vel += acc * dt/2              #Another Half kick to the velocity 
    
        t += dt                        #update time

        print('Iteration number:',k)   #Prints current itretation label

        pos_save[:,:,k+1] = pos        #Save current position matrix 
        KE_save[k+1] = KE              #Save current KE value 
        PE_save[k+1] = PE              #Save current PE value 

    if save_matrix == True:   
        np.save('MPI_position_data', pos_save)
    


for k in range(Nt): 
    if rank != 0: #Code running on the slave cores 
        pos = np.empty((N,3),dtype=np.float) 

        comm.Recv(pos, source = 0) #Reciving the current state of the system from the master core

        a,inv_P = acc(pos, mass, G, rank, row_per_thread) #Each core computes the acceleration of every particle assigned to it
        comm.Send(a, dest = 0) #Sending this partial acceleration matrix to the master core
        comm.Send(inv_P, dest = 0) #Sending this partial distance matrix to the master core (KE calculation)


        
