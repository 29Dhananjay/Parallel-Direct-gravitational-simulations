'''
This code simulates the motion of N objects using OpenMP, and saves the position matrix (stacked timesteps) locally.
'''


from cython.parallel cimport parallel
cimport openmp

import os 

os.environ["MKL_NUM_THREADS"] = "1" 
os.environ["NUMEXPR_NUM_THREADS"] = "1" 
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np 
cimport numpy as np


softening = 0 #Softening term to avoid singularities


def acc(np.ndarray pos, np.ndarray mass, double G, double rank, double row_per_thread):
    """[summary]: This function calculates the acceleration matrix

    Args:
        pos            : N x 3 matrix of positions 
        mass           : N x 1 vector of masses
        G              : Gravitational Constant
        rank           : Rank of the processor
        row_per_thread : Number of particles each rank handles 

    Returns:
        the acceleration and the distance matrices (PE calculations)

    """

    cdef:
        np.ndarray x, y, z, dx, dy, dz, inv_L, inv_R, ax, ay, az

    # positions r = [x,y,z] for all particles
    x = pos[:,0:1]
    y = pos[:,1:2]
    z = pos[:,2:3]
    

    dx = x.T - x[(rank-1)*row_per_thread:row_per_thread*rank]  #Vectorised distance calculation x-axis
    dy = y.T - y[(rank-1)*row_per_thread:row_per_thread*rank]  #Vectorised distance calculation y-axis
    dz = z.T - z[(rank-1)*row_per_thread:row_per_thread*rank]  #Vectorised distance calculation z-axis


    inv_L = (dx**2 + dy**2 + dz**2 + softening**2) # 1/r^2 distance

    inv_R = np.copy(inv_L)

    inv_R[inv_R>0] = inv_R[inv_R>0]**(-1.5)

    inv_P = np.sqrt(inv_L)

    inv_P[inv_P>0] = inv_P[inv_P>0]**(-1)

    ax = G * (dx * inv_R) @ mass #acceleration x component 
    ay = G * (dy * inv_R) @ mass #acceleration y component 
    az = G * (dz * inv_R) @ mass #acceleration z component 

    return np.hstack((ax,ay,az)),inv_P
    

def get_acc(int threads, np.ndarray pos, np.ndarray mass, double G):
    """[summary]: This function calculates acceleration in parallel  

    Args:
        threads    : Number of threads made available. 
        pos        : N x 3 matrix of positions 
        mass       : N x 1 vector of masses
        G          : Gravitational Constant

    Returns:
        the acceleration, and the distance matrix computed in parallel
    """

    cdef:
        double row_per_thread
        np.ndarray stack, I, a, b, rank_sort, value1, value2
        int rank, i, index


    row_per_thread = int(len(mass)/(threads-1)) 
    
    I_stack = []     #Empty distance list which will be updated by other cores (PE calculation)
    acc_stack = []   #Empty acceleration list which will be updated by other cores
    rank_stack = [] 

    stack = np.zeros((1,3))
    I = np.zeros((1,len(mass)),dtype=np.float)    


    with nogil, parallel(num_threads=threads): #Initiating parallel threads 
        rank = openmp.omp_get_thread_num()
        with gil:
            for i in range(1,threads):
                if rank == i:
                    a,b = acc(pos, mass, G, i, row_per_thread) #Each core computes the acceleration of every particle assigned to it
                    acc_stack.append(a)                        #Appending the acceleration calculation to the acceleration list
                    I_stack.append(b)                          #Appending the distance calculation to the distance list
                    rank_stack.append(rank) 
                   

    #the code below sorts the acceleration matrix in the correct order
    rank_sort = np.array(rank_stack).argsort()
    for i in range(len(rank_sort)):
        index = rank_sort[i]
        value1 = acc_stack[index]
        stack = np.vstack((stack,value1))
        value2 = I_stack[index]
        I = np.vstack((I,value2))
    
    return stack[1:,:],I[1:,:] #Returning the acceleration (force) and the distance (PE) matrices


        


def main(int threads, np.ndarray pos, np.ndarray vel, np.ndarray mass, save_matrix):
    """[summary]: This function performs the simulation. 

    Args:
        threads  : Number of threads made available. 
        pos      : N x 3 matrix of positions 
        vel      : N x 3 vector of velocities 
        mass     : N x 1 vector of masses
       
    """


    cdef:
        int N, i, rank, Nt
        double row_per_thread, G, i_year, KE, PE, dt, tEnd
        np.ndarray pos_save, KE_save, PE_save, acceleration

    
    N         = len(mass)     # Number of particles
    t         = 0             # current time of the simulation
    tEnd      = 1             # End time in years 
    dt        = 0.0000547945  # timestep
    G         = 4*np.pi*np.pi # Newton's Gravitational Constant

    #Initial kinetic and potential energies 
	#KE and PE are approximated as 0 just for the first time step 
    KE, PE  = 0,0 
    

    Nt = int(np.ceil(tEnd/dt)) #Number of timesteps

    #The code below Saves current state of the system. (position matrices corresponding to different time steps are stacked on top of one another)
    pos_save = np.zeros((N,3,Nt+1))
    pos_save[:,:,0] = pos
    KE_save = np.zeros(Nt+1)
    KE_save[0] = KE
    PE_save = np.zeros(Nt+1)
    PE_save[0] = PE
    t_all = np.arange(Nt+1)*dt


    for i in range(Nt): #Iterating over time steps

        acceleration,I = get_acc(threads, pos, mass, G) #Calculating acceleration in parallel
        
        vel += acceleration * dt/2   #Half kick to the velocity

        pos += vel * dt              #Updating the position vector

        vel += acceleration * dt/2   #Another Half kick to the velocity 

        t += dt                      #update time
        
        print('Iteration number:',i) #Prints current itretation label
        
        KE = 0.5 * np.sum(np.sum(mass*vel**2))
        PE = G * np.sum(np.sum(np.triu(-(mass*mass.T)*I,1)))
        

        pos_save[:,:,i+1] = pos      #Save current position matrix
        KE_save[i+1] = KE            #Save current KE value 
        PE_save[i+1] = PE            #Save current PE value 
       
    '''remove the '#' to save the position matrix locally'''
    if save_matrix == True:
        np.save('OpenMP_position_data', pos_save)










