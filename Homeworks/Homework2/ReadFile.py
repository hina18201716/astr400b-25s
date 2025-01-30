#!/usr/bin/env python
# coding: utf-8

# In[9]:


# Import All  Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants


# In[10]:


# Function to read file 
def Read(filename): 
    """ This function read a file formatted as 
            Time (Myr)
            # of Particle
            units
            Header names
            - type Mass(10e M_sun) Position(kpc) Velocity(km/s)
            Coordinate is layed s.t. origin in at the Galactic Center. 
        Input: filename
        Output: Time(Myr)
    """
    # open file 
    file = open(filename, 'r')

    line1 = file.readline() # read a line from the file
    label, value = line1.split() # convert from one string to two 
    time = float(value) * u.Myr # convert unit
    line2 = file.readline() # read a line from the file
    label2, n = line2.split() # convert from one string to two
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    # close file
    file.close()
    
    return time, n, data

