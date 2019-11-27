# -*- coding: utf-8 -*-
#This code was written to orient olivine diffusion profiles relative to their crystal orientations
#Currently it is set up for the EDAX EBSD at the American Museum of Natural History.
#Henry Towbin
# Rotation Matrix Code modified from Britton et al 2016 "Tutorial: Crystal orientations and EBSD - Or which way is up?"
 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math as m
import matplotlib.image as image
#import mplstereonet # this is commented out becaus it requires a seperate module 

#alpha = np.deg2rad(alpha) # convert angles into radians
#beta  = np.deg2rad(beta)
#gamma = np.deg2rad(gamma)
#
## Euler Angles in Bunge convention
#phi1=float(BRKR_phi1)*m.pi/180.0
#Phi =float(BRKR_Phi) *m.pi/180.0
#phi2=float(BRKR_phi2)*m.pi/180.0


#------------------------------------------------------------------------------
# Rotation matrices 
#------------------------------------------------------------------------------
#  see e.g.
#  J. B. Kuipers "Quaternions and Rotation Sequences", 
#  Princeton University Press, 1999

def Rx(RotAngle):
    '''
    provides the X axis (e1) rotation matrix in cartesian systems,
	input "RotAngle" in radians

    meaning for a COLUMN of VECTORS: transformation matrix U for VECTORS
    this matrix rotates a set of "old" basis vectors 
    $(\vec{e1},\vec{e2},\vec{e3})^T$ (column) by +RotAngle (right hand rule) 
    to a new set of basis vectors $(\vec{e1}',\vec{e2}',\vec{e3}')^T$ (column)

    meaning for a COLUMN of vector COORDINATES:
    (1) (N_P_O):    coordinates of fixed vector in a "new" basis that is 
                    rotated by +RotAngle (passive rotation)
    (2) (O_P_O):    active rotation of vector coordinates in the same 
                    "old" basis by -RotAngle
    '''
    mat=np.matrix([ [1,                 0,               0 ],
                    [0,  np.cos(RotAngle), np.sin(RotAngle)],
                    [0, -np.sin(RotAngle), np.cos(RotAngle)]])
    return mat
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def Rz(RotAngle):
    '''
    provides the Z axis (e3) rotation matrix in cartesian systems,
	input "RotAngle" in radians

    meaning for a COLUMN of VECTORS: transformation matrix U for VECTORS
    this matrix rotates a set of "old" basis vectors 
    $(\vec{e1},\vec{e2},\vec{e3})^T$ (column) by +RotAngle (right hand rule) 
    to a new set of basis vectors $(\vec{e1}',\vec{e2}',\vec{e3}')^T$ (column)

    meaning for a COLUMN of vector COORDINATES:
    (1) (N_P_O):    coordinates of fixed vector in a "new" basis that is 
                    rotated by +RotAngle (passive rotation)
    (2) (O_P_O):    active rotation of vector coordinates in the same 
                    "old" basis by -RotAngle
    '''
    mat=np.matrix([[ np.cos(RotAngle) , np.sin(RotAngle), 0 ],
                   [-np.sin(RotAngle) , np.cos(RotAngle), 0 ],
                   [                0 ,                0, 1 ]])
    return mat
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
def CalcStructureMatrix(a=1.0,b=1.0,c=1.0,
                        alpha=np.deg2rad(90.0),
                        beta=np.deg2rad(90.0),
                        gamma=np.deg2rad(90.0)):
    '''
    computes the structure matrix from lattice parameters
    input angles in RADIANS    
    convention: McKie&McKie, "Essentials of Crystallography", 1986
    '''

    ca = m.cos(alpha)
    sa = m.sin(alpha)
    cb = m.cos(beta)
    cg = m.cos(gamma)
    
    ax = a * m.sqrt(1.0+2.0*ca*cb*cg-(ca*ca+cb*cb+cg*cg))/sa
    ay = a * (cg-ca*cb)/sa
    az = a * cb
    
    by = b * sa
    bz = b * ca
    
    cz = c
    
    StructureMat=np.matrix([[ax , 0,  0 ],
                            [ay , by, 0 ],
                            [az , bz, cz]])
    return StructureMat
    
    
#------------------------------------------------------------------------------
# STRUCTURE MATRICES
#------------------------------------------------------------------------------
# direct structure matrix A
# A gives the direct lattice basis vector coordinates (column vectors)
# in the Cartesian system
#A = CalcStructureMatrix(a,b,c,alpha,beta,gamma)
# reciprocal structure matrix A+: transpose of inverse of A
# A+ gives to reciprocal lattice basis vector coordinates (column vectors)
# in the Cartesian system
#Aplus = (A.I).T
#------------------------------------------------------------------------------

# rotation of Sample System to Cartesian Crystal System = ORIENTATION OF GRAIN
#

#I am using Row Vectors in accoradnace  with Britton et al 2016 
#it is important to remember the order of multiplcation row row vectors xA=B not Ax=b(column vectors)

#U_O = Rz(phi2)  * Rx(Phi) * Rz(phi1)


def Crystal2Sample(uvw,phi1,Phi,phi2):
    uvw = np.matrix(uvw)
    
    U_O = Rz(np.deg2rad(phi2))  * Rx(np.deg2rad(Phi)) * Rz(np.deg2rad(phi1))
    
    Structmat = CalcStructureMatrix(a=1,b=1,c=1,  # right now its for the cubic system. I need to determine if we need vectors in the basis of for the unit cell of its just the raw angle 
                        alpha=np.deg2rad(90.0),     #olivine a=4.78,b=10.25,c=6.3
                        beta=np.deg2rad(90.0),
                        gamma=np.deg2rad(90.0))
                        
    #Sample_xyz= uvw * Structmat.T * U_O # commented out to remove struct matrix
    Sample_xyz= uvw * U_O
    
    return Sample_xyz


    
    
def Sample2Crystal(xyz,phi1,Phi,phi2):
    uvw = np.matrix(xyz)
    
    U_O = Rz(np.deg2rad(phi2))  * Rx(np.deg2rad(Phi)) * Rz(np.deg2rad(phi1))
    
    Structmat = CalcStructureMatrix(a=1,b=1,c=1,  # right now its for the cubic system. I need to determine if we need vectors in the basis of for the unit cell of its just the raw angle 
                        alpha=np.deg2rad(90.0),     #olivine a=4.78,b=10.25,c=6.3
                        beta=np.deg2rad(90.0),
                        gamma=np.deg2rad(90.0))
                        
    Sample_xyz= uvw * U_O.T #* (Structmat.T).I #commented out to remove structure matrix 
     
    return Sample_xyz
   
   
   
def Crystal2Sample2(uvw,phi1,Phi,phi2):
    uvw = np.matrix(uvw)
    
    U_O = Rz(np.deg2rad(phi2))  * Rx(np.deg2rad(Phi)) * Rz(np.deg2rad(phi1))
    
    Structmat = CalcStructureMatrix(a=1,b=1,c=1,  # right now its for the cubic system. I need to determine if we need vectors in the basis of for the unit cell of its just the raw angle 
                        alpha=np.deg2rad(90.0),     #olivine a=4.78,b=10.25,c=6.3
                        beta=np.deg2rad(90.0),
                        gamma=np.deg2rad(90.0))
                        
    #Sample_xyz= uvw * Structmat.T * U_O # commented out to remove struct matrix
    Sample_xyz=  U_O * uvw.T
    
    return Sample_xyz   
   
def Sample2Crystal2(xyz,phi1,Phi,phi2):
    uvw = np.matrix(xyz)
    
    U_O = Rz(np.deg2rad(phi2))  * Rx(np.deg2rad(Phi)) * Rz(np.deg2rad(phi1))
    
    Structmat = CalcStructureMatrix(a=1,b=1,c=1,  # right now its for the cubic system. I need to determine if we need vectors in the basis of for the unit cell of its just the raw angle 
                        alpha=np.deg2rad(90.0),     #olivine a=4.78,b=10.25,c=6.3
                        beta=np.deg2rad(90.0),
                        gamma=np.deg2rad(90.0))
                        
    Sample_xyz=  U_O.T * uvw.T #* (Structmat.T).I #commented out to remove structure matrix 
     
    return Sample_xyz
    
def profile_angle(theta): # Theta in degrees from the positive X axis
    x = np.sin(np.deg2rad(theta)) # this configuration has the x axis oriented to the north  (in OIM A1 axis)
    y = -np.cos(np.deg2rad(theta)) # (in OIM A2 axis)This has a negative sign because our SEM is oriented left to right. I may need to change this if the Euler angles dont represent the orientation corrected 
    z = 0.0 #(in OIM A3 axis)
    return [x,y,z]
    
def Shortest_angle(uvw,axis):
    
    dot = np.dot(uvw,axis.T)
    vector_norm = np.linalg.norm(uvw)
    axis_norm = np.linalg.norm(axis)
    
    short_angle = np.arccos(dot/(axis_norm * vector_norm))
    
    short_angle =  np.rad2deg(short_angle)
    
    return short_angle
    
    
def vector_direction(theta,phi1,Phi,phi2): 
#Sample Rotated using Crystal 2 Sample (Works as Expected for AMNH EBSD)
# """Theta: Input the angle of your profile relative to the positive x axis (East) on your image
#Input Euler angles for your crystal"""
    
    uvw = np.matrix(profile_angle(theta))
    
    #uvw = Sample2Crystal(xyz,phi1,Phi,phi2)
    

    Ac = Crystal2Sample(np.matrix([1,0,0]),phi1,Phi,phi2)
    Bc = Crystal2Sample(np.matrix([0,1,0]),phi1,Phi,phi2)# B and A2 aligned to the East
    Cc = Crystal2Sample(np.matrix([0,0,1]),phi1,Phi,phi2)
    
    print "Apha (to 100) :" + str(min(Shortest_angle(uvw,Ac),Shortest_angle(uvw,-Ac)))
    print "Beta (to 010):" + str(min(Shortest_angle(uvw,Bc),Shortest_angle(uvw,-Bc)))
    print "Gamma (to 001):" + str(min(Shortest_angle(uvw,Cc),Shortest_angle(uvw,-Cc)))
    
    
    

    
    
def vector_direction2(theta,phi1,Phi,phi2): # Vector Rotated using Sample to Crystal (Works as Expected for AMNH EBSD)
    xyz = profile_angle(theta)
    uvw = Sample2Crystal(xyz,phi1,Phi,phi2)
    
    
    
    #Ac = Sample2Crystal(np.matrix([0,1,0]),phi1,Phi,phi2)
    #Bc = Sample2Crystal(np.matrix([1,0,0]),phi1,Phi,phi2)#B and A2 aligned to the East 
    #Cc = Sample2Crystal(np.matrix([0,0,1]),phi1,Phi,phi2)
    Ac = np.matrix([1,0,0])
    Bc = np.matrix([0,1,0])# B and Y aligned to the East
    Cc = np.matrix([0,0,1])
    print "Apha (to 100):" + str(min(Shortest_angle(uvw,Ac),Shortest_angle(uvw,-Ac)))
    print "Beta (to 010):" + str(min(Shortest_angle(uvw,Bc),Shortest_angle(uvw,-Bc)))
    print "Gamma (to 001):" + str(min(Shortest_angle(uvw,Cc),Shortest_angle(uvw,-Cc)))
    
    
Seperate module for plotting stereonets required for the following code. 
#def trend_plunge(vector):
#  xy_hyp = m.sqrt(vector[0,0]**2+vector[0,1]**2) # calculates the the hypotenuse of the XY vector
#  plunge = m.degrees(m.atan2(vector[0,2],xy_hyp)) # Positive dip is postive z meaning up from sample surface #stereonet should be plotted in only 1 hemisphere
#  trend = m.degrees(m.atan2(-vector[0,1],vector[0,0]))
#  return trend, plunge
#  
#
#def plot_olivine_stereonet(theta,phi1,Phi,phi2):
#    """Theta: Input the angle of your profile relative to the positive x axis (East) on your image
#    Input Euler angles for your crystal
#    stereonets may be 180˚ rotations from what you see in OIM data collection because it takes into account rotations of the sample stage"
#    """
#    #Plots a Lower Hemisphere Stereonet
#    
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='stereonet')
#    
#    uvw = np.matrix(profile_angle(theta))
#    trend, plunge = trend_plunge(uvw)
#    ax.line(plunge,trend,'go')
#    trend, plunge = trend_plunge(-uvw)
#    ax.line(plunge,trend,'go')
#    
#    Ac = Crystal2Sample(np.matrix([1,0,0]),phi1,Phi,phi2)
#    trend, plunge = trend_plunge(Ac)
#    ax.line(plunge,trend,'bo')
#    trend, plunge = trend_plunge(-Ac)
#    ax.line(plunge,trend,'bo')
#    
#    Bc = Crystal2Sample(np.matrix([0,1,0]),phi1,Phi,phi2)# B and A1 aligned to the East
#    trend, plunge = trend_plunge(Bc)
#    ax.line(plunge,trend,'ro')
#    trend, plunge = trend_plunge(-Bc)
#    ax.line(plunge,trend,'ro')
#
#    Cc = Crystal2Sample(np.matrix([0,0,1]),phi1,Phi,phi2)
#    trend, plunge = trend_plunge(Cc)
#    ax.line(plunge,trend,'ko')
#    trend, plunge = trend_plunge(-Cc)
#    ax.line(plunge,trend,'ko')
#    
#    ax.grid()
#
#    plt.show()