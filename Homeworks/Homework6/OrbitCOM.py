

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass




def OrbitCOM(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
          galaxy: 'string'
              galaxy name
        start: 'int'
            the number of first snapshot
        end: 'int'
            the number of the last snapshot
        n : 'int'
            intervels of COM calculation
    outputs: 
        
    
    """
    
    # compose the filename for output
    fileout = "Orbit_" + galaxy+ ".txt"
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1 
    if galaxy == "M33": 
        volDec = 4
    else: 
        volDec = 2 
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    snap_ids = np.arange(start, end, n)
    if snap_ids.size == 0: 
          print("no elements")

    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros((snap_ids.size, 7))
    
    # a for loop 
    for  i, snap_id in enumerate(snap_ids):
        
        # compose the data filename (be careful about the folder)
        ilbl = '000' + str(i)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        filename="%s_"%(galaxy) + ilbl + '.txt'
        
        # Initialize an instance of CenterOfMass class, using disk particles ptype == 2 
        COM = CenterOfMass(filename, 2)
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        COM_pos = COM.COM_P(0.1 ,volDec)
        COM_vel = COM.COM_V( COM_pos[0], COM_pos[1], COM_pos[2])
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        
        time = COM.time.value/ 1000
        orbit[i] = np.array([time,  COM_pos[0].value, COM_pos[1].value, COM_pos[2].value, COM_vel[0].value, COM_vel[1].value, COM_vel[2].value])
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 




# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt




# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  




# Determine the magnitude of the relative position and velocities 

# of MW and M31

# of M33 and M31




# Plot the Orbit of the galaxies 
#################################




# Plot the orbital velocities of the galaxies 
#################################

