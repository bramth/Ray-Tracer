# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 10:25:57 2019

@author: R. Wildeboer, June 2014
@ported: B. ter Huurne, February 2019
"""

import numpy as np
import matplotlib.pyplot as plt
from optics_lib import snells,propagate2flat,propagate2sph,calc_n_K,calc_n_F,wavelength_to_rgb

def calc_dist(d,z,index):
    dist = sum(d[:index+1]) - sum(z[:index+1])
    return dist

def plot_optical_axis(ax,d):
    ax.plot([0,sum(d)],[0,0],':k')

def plot_lens(ax,D,R,d):
    Y = np.linspace(-D/2,D/2,30)                # spherical interface height
    X = sum(d)+R-np.sqrt(R**2-Y**2)*R/np.abs(R) # spherical interface location
    ax.plot(X,Y,':k')

def plot_flat(ax,D,d):
    sum_d = sum(d)
    ax.plot([sum_d,sum_d],[-D/2,D/2],':k')             # flat surface
    
def plot_ray(ax,Z,R,color):
    X = [0]
    for item in Z:
        X.append(item+X[-1])        
     
    ax.plot(X,R,color=color)            # draw rays
    
    

# Initialization
D = 60                          # component diameter [mm]
#L = [400,550,700]               # wavelength [nm]
L = [550]
f = 100                         # focal length of set [mm]
NA = 0.3                        # numerical aperature 
double = True

# Figure 
fig = plt.figure()              
ax = fig.add_subplot(111)
ax.set_xlabel('optical axis z [mm]')
ax.set_ylabel('height r [mm]')
ax.set_aspect('equal')


for l in L: 
    n0 = 1                          # vacuum refractive index [mm]
    n1 = calc_n_K(l)                # lens 1 refractive index [mm]
    n2 = calc_n_F(l)                # lens 2 refractive index [mm]
    
    if double:    
        R = [47,-47,0]
        d = [50,25,10,110]
        n = [(n0,n1),(n1,n2),(n2,n0)]
    else:
        R = [55,0]
        d = [50,25,110]
        n = [(n0,n1),(n1,n0)]

    color = wavelength_to_rgb(l)  
    
    RR = []
    THETA = [] 
    Z = []
    
    for jdx,r_0 in enumerate(np.linspace(-(D/2),(D/2-2*D/10),20)):    # starting positions from axis
        z = []
        r = [r_0]
        theta = [0]     # incident angle of beam
        
        # Interfaces
        for idx,R_lens in enumerate(R):
            if R_lens == 0: 
                dist = calc_dist(d,z,idx)
                z_temp,r_temp,theta_n,theta_z = propagate2flat(r[-1],theta[-1],dist,D)
                #if z_temp == -1: continue
                z.append(z_temp)
                r.append(r_temp)
                theta.append(snells(theta_n,theta_z,n[idx]))
            else:
                dist = calc_dist(d,z,idx)
                z_temp,r_temp,theta_n,theta_z = propagate2sph(r[-1],theta[-1],dist,D,R[idx])
                #if z_temp == -1: continue
                z.append(z_temp)
                r.append(r_temp)
                theta.append(snells(theta_n,theta_z,n[idx]))
        
        
        # EVALUATION
        RR.append(r[-1])
        THETA.append(-np.tan(theta[-1]))
        Z.append(sum(z))
            
        # Final 'end' interface
        dist = calc_dist(d,z,idx+1)
        z_temp,r_temp,theta_n,theta_z = propagate2flat(r[-1],theta[-1],dist,D)
        z.append(z_temp)
        r.append(r_temp)
    
        # VISUALIZATION: Plot rays
        plot_ray(ax,z,r,color)
        
    # EVALUATION C'NTD
    g = 1
    cross = []
    for j in range(jdx-1):
        for k in range(j+1,jdx):
            cross.append((RR[j]-RR[k]+THETA[k]*(Z[j]-Z[k]))/(THETA[j]-THETA[k]) + Z[k] - d[0] - np.sum(d[1:-1])/2)
            g += 1
    f = np.mean(cross)
    df = np.std(cross)
    
    # PRINT FOCAL LENGTH
    print(f'f = {f}, df = {df}')
       


# VISUALIZATION: Optical system
plot_optical_axis(ax,d)
for idx,R_lens in enumerate(R):
    if R_lens == 0:
        plot_flat(ax,D,d[:idx+1])                  # flat surface
    else:
        plot_lens(ax,D,R_lens,d[:idx+1])       # spherical surface