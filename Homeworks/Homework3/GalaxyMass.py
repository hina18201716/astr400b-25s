#!/usr/bin/env python
# coding: utf-8

# In[5]:


#import all library needed 

import numpy as np
import astropy.units as u
from ReadFile import Read
from tabulate import tabulate


# In[6]:


def ComponentMass(filename, Ptype): 
    """
    This function calculate the total mass of one component in the galaxy.
    Input: 
        Filename 
        Particle name: 1-Halo type, 2-Disk type, 3-Bulge type
    Output: Mass of the componets in the galaxy in 3 decimal place (10e12 M_sun)
    """
    time, Pn, data = Read(filename) # read file 
    
    Tmass = 0 # set total mass to be zero 
    for i in data: 
        if i['type'] == Ptype: 
            Tmass += i['m'] *1e10 * u.M_sun   # Total mass in 1e10 solar mass
    
    Tmass =Tmass / 1e12 #convert into 1e12 solar mass 
     
    return np.round( Tmass, 3) # round up to 3 significant digits and return


# In[7]:


##### optional make table in python
header = ['Galaxy Name', 'Halo Mass', 'Disk Mass', 'Bulge Mass', 'Total Mass']

MW = ['MW']
M31 = ['M31']
M33 = ['M33']

tot_MW = 0
tot_M31 = 0
tot_M33 = 0


# In[8]:


for i in range(0,3):
    MW_mass = ComponentMass("MW_000.txt", i+1)
    M31_mass = ComponentMass("M31_000.txt", i+1)
    M33_mass = ComponentMass("M33_000.txt", i+1)    
    
    MW.append(MW_mass)
    M31.append(M31_mass)
    M33.append(M33_mass)
    
    tot_MW += MW_mass 
    tot_M31 += M31_mass
    tot_M33 += M33_mass
    


# In[9]:


MW.append(np.round(tot_MW,3))
M31.append(np.round(tot_M31,3))
M33.append(np.round(tot_M33,3))


# In[10]:


table = np.array([header, MW, M31, M33])


# In[11]:


print(tabulate(table))

