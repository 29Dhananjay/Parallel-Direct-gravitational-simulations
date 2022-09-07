import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astroquery.jplhorizons import Horizons
import pandas as pd

startdate = "01/01/2015"


#Sizes of planets
rad = np.array([0.5,0.0035068,0.0086989,0.0091577,0.0048721,
0.20049,0.113703,0.036455,0.035392,0.0017081,
0.0024969,0.0026184,0.0022435])


moon_size = 10
asteroid_size = 10

#Assigning an equal size to all natural satellites and asteroids 
jovian_rad = np.empty(67)
jovian_rad.fill(0.0017081/moon_size)

saturn_rad = np.empty(80)
saturn_rad.fill(0.0017081/moon_size)

uranus_rad = np.empty(27)
uranus_rad.fill(0.0017081/moon_size)

neptune_rad = np.empty(13)
neptune_rad.fill(0.0017081/moon_size)

asteroid_belt = np.empty(280)
asteroid_belt.fill(0.0017081/asteroid_size)

rad = np.concatenate((rad,jovian_rad,saturn_rad,uranus_rad,neptune_rad,asteroid_belt))


#Color of the planets
color = ['darkorange', 'lightgray', 'goldenrod', 'royalblue', 'firebrick', 
'burlywood', 'wheat', 'cornflowerblue', 'mediumblue', 'ghostwhite',
'ghostwhite','ghostwhite','ghostwhite']

jovian_color = [
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite']

saturn_color = [
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite']

uranus_color = [
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite']

neptune_color = [
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite','ghostwhite','ghostwhite',
'ghostwhite','ghostwhite','ghostwhite']#,'ghostwhite']


asteroid_belt_color = ['ghostwhite' for i in range(280)]


color = color+jovian_color+saturn_color+uranus_color+neptune_color+asteroid_belt_color


pos_save = np.load('path_to_position_matrices')

#pos_save = np.load('2015_simulation.npy')   #use this to plot the animation for the year 2015. 

length = np.shape(pos_save)[2] #Number of iterations 

fig = plt.figure(figsize=(8,8),dpi=125)
ax1 = plt.subplot()
fig.set_facecolor('black') 
ax1.set_facecolor('black') 
ax1.set_facecolor((0, 0, 0))

time_step = 30
#Plotting every 30th timestep
for i in range(0,length,time_step):
    
	plt.sca(ax1)
	plt.cla()
	pos = pos_save[:,:,i] 
	x = pos[:,0]
	y = pos[:,1]
	
	enddate = pd.to_datetime(startdate) + pd.DateOffset(days=i/50) #Date 
	enddate = enddate.date()


	ax1.scatter(x,y,s=65*rad,color = color) #Plotting all bodies 
	ax1.text(-2.6,5.9,'Date: '+ str(enddate),color = 'w',fontsize=12)    
	font_size = 6 
	#labeling major bodies 
	ax1.text(x[0]+0.1,y[0]+0.2,'Sun',color = 'w',fontsize=5)          
	ax1.text(x[1]+0.1,y[1]+0.1,'Mercury',color = 'w',fontsize=5)     
	ax1.text(x[2]+0.1,y[2]+0.1,'Venus',color = 'w',fontsize=5)        
	ax1.text(x[3]+0.1,y[3]+0.1,'Earth',color = 'w',fontsize=5)
	ax1.text(x[4]+0.1,y[4]+0.1,'Mars',color = 'w',fontsize=font_size)
	ax1.text(x[5]+0.1,y[5]+0.1,'Jupiter',color = 'w',fontsize=font_size)
	ax1.text(x[6]+0.1,y[6]+0.1,'Saturn',color = 'w',fontsize=font_size)
	ax1.text(x[7]+0.1,y[7]+0.1,'Uranus',color = 'w',fontsize=font_size)
	ax1.text(x[8]+0.1,y[8]+0.1,'Neptune',color = 'w',fontsize=font_size)
	ax1.text(x[9]+0.1,y[9]+0.1,'Pluto',color = 'w',fontsize=font_size)
 

	#Limits of the plot -10.3 AU < X < 10.3 AU and -10.3 AU < Y < 10.3 AU
	ax1.set(xlim=(-10.3, 10.3), ylim=(-10.3, 10.1))
	ax1.grid(False)
	#Pausing between frames
	plt.pause(0.000000000001)

plt.show()
