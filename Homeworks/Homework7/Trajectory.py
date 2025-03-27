
# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 




# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass2 import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
import GalaxyMass

# # M33AnalyticOrbit




class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, galaxy):
        """ determine the acceleration M33 feels from M31 
        and integrate its current position and velocity forwards in time 
        PRAMETER
        ---------
            galaxy (string)
                galaxy name
        """
        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        fileout = "Orbit_" + galaxy+ ".txt"

        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        COM_M33 = CenterOfMass("M33_000.txt", 2)

        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        COMP_M33_unit = COM_M33.COM_P(0.1 ,4)
        COMP_M33 = COMP_M33_unit.value # 4 because M33, for MW M31, use 2

        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        COMV_M33 = COM_M33.COM_V( COMP_M33_unit[0], COMP_M33_unit[1], COMP_M33_unit[2] ).value


        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        COM_M31 = CenterOfMass("M31_000.txt", 2)

        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        COMP_M31_unit = COM_M31.COM_P(0.1 ,2) # 4 because M33, for MW M31, use 2
        COMP_M31 = COMP_M31_unit.value # 4 because M33, for MW M31, use 2
        
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        COMV_M31 = COM_M31.COM_V( COMP_M31_unit[0], COMP_M31_unit[1], COMP_M31_unit[2] ).value


        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = COMP_M33 - COMP_M31
        self.v0 = COMV_M33 - COMV_M31

        ### get the mass of each component in M31 
        # ptype 1 - Halo, 2 - disk, 3 - bulge 
        ### disk
        # **** self.rdisk = scale length (no units)
        self.rdisk = 5

        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = GalaxyMass.ComponentMass("M31_000.txt", 2) *1e12

        ### bulge
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1

        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = GalaxyMass.ComponentMass("M31_000.txt", 3) *1e12

        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 60
        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = GalaxyMass.ComponentMass("M31_000.txt", 1) *1e12

        self.filename = "output.txt"
    
    def HernquistAccel(self, M, r_a, r): # it is easiest if you take as an input the position VECTOR 
        """ Function to calculate Halo and Bulge Hernquist Acceleration defined as -G*M/(rmag *(ra + rmag)**2) * r, 
        where the last r is a vector. 
        
        Input 
            M (float):
                mass of componets to calcualte accerelation
            ra (float):
                scale length
            r (np.array): 
                position vector of components center of mass
        Output 
            Hern (np.array): 
                vector of accerelation of particle 
        
        """
        
        ### **** Store the magnitude of the position vector
        mag2 = np.sum(r**2)
        rmag = np.sqrt(mag2)
        
        ### *** Store the Acceleration
        Hern =  - (self.G * M ) * r/ (rmag * (r_a + rmag)**2 ) 
        #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern.value
    
    
    
    def MiyamotoNagaiAccel(self, M, rd, r):# it is easiest if you take as an input a position VECTOR  r 
        """Function to calculate Disk Acceleration defined as -G*M/(R**2 + B**2) * r * [1, 1, B/np.sqrt(z**2 + zd**2 )], 
        where the last r is a vector , R = sqrt(x**2+y**2), B = rd+sqrt(z**2+zd**2). 
        
        Input 
            M (float):
                mass of componets to calcualte accerelation
            rd (float):
                scale length
            r (np.array): 
                position vector of components center of mass
        Output 
            MNAccel (np.array): 
                vector of accerelation of particle 
         """

        
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whole thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        
        zd = rd /5.0
        R = np.sqrt(r[0]**2 + r[1]**2)
        B = rd + np.sqrt( r[2]**2 + zd**2 )
        z = r[2]
        
        c = np.array( [1, 1, B/np.sqrt(z**2 + zd**2)] )
        a = - (self.G * M) / (R**2 + B**2 )**(3/2)
        
        MNAccel = a* r *c
                       
        return MNAccel.value
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r): # input should include the position vector, r
        """ Function to calculate acceletation of M31. 
        PARAMETER
        ---------
            r (np.array of float)
                position vector
        """

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
                       
        halo    = self.HernquistAccel        (self.Mhalo, self.rhalo, r)
        disk    = self.MiyamotoNagaiAccel    (self.Mdisk, self.rdisk, r)
        bulge   = self.HernquistAccel        (self.Mbulge, self.rbulge, r)
        
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return halo + disk + bulge
    
    
    
    def LeapFrog(self, dt, r, v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """ One step in LeapFrog method. 
        Function to update the position and velocity by using kinemtatic equation. 
            x = x0 + vt + 1/2 at^2
            v = v0 + at 
            
        PAPAMETER
        ---------
            dt (float)
                step size
            r (np array of float)
                position vector 
            v (np array of float)
                velocity vector
        
        OUTPUT 
        -------
            rnew (np array of float)
                updated position vector 
                
            vnew (np array of float)
                updated velocity vector 
        """
        
        # predict the position at the next half timestep
         
        rhalf = r + 0.5 * dt * v
        
        # predict the final velocity at the next tip using the acceleration field at the rhalf position 
        ahalf = self.M31Accel(rhalf)
        
        # vnew = np.array([])

        vnew = v + dt * ahalf
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew * 0.5 * dt   
        return rnew, vnew # **** return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """FUnction to solve equation of motion by running LeapFrog iterratively. 
        This will create text file that contains position and velocity vector at each time step. 
        PARAMETER
        -----------
            t0 (float)
                initial time 
            dt (float)
                step size
            tmax (float)
                end time
        
        """
        # initialize the time to the input starting time
        t = t0
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros( (int(tmax/dt)+2 , 7) )
    
        r = self.r0 
        v = self.v0 
        # # initialize the first row of the orbit
        # orbit[0] = t0, *tuple(r), *tuple(v)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
      
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while ( t < tmax):  # as long as t has not exceeded the maximal time 
            print("Step: ", i)
            # **** advance the time by one timestep, dt
            t += dt
            # **** store the new time in the first column of the ith row
            orbit[i,0] = t 
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            r = orbit[i-1, 1:4]
            v = orbit[i-1, 4:7]
            
            rnew, vnew = self.LeapFrog(dt, r, v)
         
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            
            orbit[i, 1:4] = rnew[0], rnew[1], rnew[2]
            orbit[i, 4:7] = vnew[0], vnew[1], vnew[2]
            
    
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function

