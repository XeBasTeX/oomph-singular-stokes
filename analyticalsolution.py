#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 10:13:11 2019

@author: bastien

This code computes a \lambda, solution of the algebric equation, and computes
the corresponding stream function near the corner. It finally plot it with its
contour and save the figures in the folder /figures/

"""

import numpy as np
import matplotlib.pyplot as plt

import scipy.optimize

#%% Find the lambda_1 value
# We find a root of f, which gives us a lambda verifying the algebric equation
# presented in the papers

# Constants definition
alpha =  (3*np.pi/2) /2
sec = lambda angle : 1/np.cos(angle)    # The secant function
K = 4   # Defined as A*sec(firstLambda - 2)*alpha

# Algebric equation we want to solve, by looking for the roots of this function
f = lambda x : x*np.tan(x*alpha) - (x - 2)*np.tan((x - 2)*alpha)

# Alternate definition of the algebric equation, this one is clearly smoother
# We stick with this form from now
# But you may easily check that they both gives the correct firstLambda
f2 = lambda x : np.sin(2*(x - 1)*alpha) + (x - 1)*np.sin(2*alpha)



# Plot the characteristic function
borneInf = 0
borneSup = 2
fVectorise = np.vectorize(f2)
xVecteur = np.arange(borneInf + 0.1, borneSup - 0.1, 1/5000)

fig = plt.figure()
plt.plot(xVecteur, fVectorise(xVecteur), label='$g \, (\lambda)$')
plt.axhline(y = 0.0, color='orange', linestyle='--', label='$x$-axis')
plt.title('Algebric equation $g$ and its roots', fontsize=18)
plt.xlabel('Exponent $\lambda$', fontsize=16)
plt.ylabel('$g \, (\lambda)$', fontsize=16)
plt.grid(True)
plt.legend()
plt.savefig('figures/characteristic_equation_plot.eps', format='eps',
            dpi=1000, bbox_inches='tight', pad_inches=0.03)
plt.show()

# What is the lambda that solves the equation?
firstLambda = scipy.optimize.newton(f2, 1.33) # 0.455516263218 # 1.54448373678
print("We find λ = ", firstLambda)


## We will know compute the streamfunction on a grid
# We define the Moffat's stream function
streamFunction = lambda r, theta : K*pow(r, firstLambda) *\
            (np.cos((firstLambda - 2)*alpha)*np.cos(firstLambda*theta/alpha)*\
             - np.cos(firstLambda*alpha)*np.cos((firstLambda - 2)*theta/alpha))
    
# For which value streamFunction(., theta) = 0 ? We plot it
fig = plt.figure()
thetaVecteur = np.arange(-np.pi, np.pi, 0.01)
plt.plot(thetaVecteur,
         streamFunction(np.ones(len(thetaVecteur)), thetaVecteur))
plt.axhline(y = 0.0, color='orange', linestyle='--', label='$x$-axis')
plt.title(r'Stream fonction along $\theta$ variation ($r$ constant)',
          fontsize=18)
plt.xlabel(r'Angle $\theta$', fontsize=16)
plt.ylabel(r'$f_\lambda (\theta)$', fontsize=16)
plt.grid(True)
plt.show()

#%% Compute and plot stream function (and velocities) on the corner

# Another and simpler definition of the stream function, in the free surface case
streamFunction2 = lambda r, theta : K*pow(r, 2 - np.pi/(2*alpha))*\
                    np.cos((theta - np.pi/4)*np.pi/(2*alpha))

# Definition of velocities
u_r = lambda r, phi : K*pow(r, firstLambda - 1)*\
    ( - firstLambda/alpha*np.cos((firstLambda - 2)*alpha)*np.sin(firstLambda*phi/alpha)*\
     + (firstLambda - 2)/alpha*np.cos(firstLambda*alpha)*np.sin((firstLambda - 2)*phi/alpha))

# u_theta = -d\psi/dr
u_theta = lambda r, phi : -firstLambda*K*pow(r, firstLambda - 1) *\
            (np.cos((firstLambda - 2)*alpha)*np.cos(firstLambda*phi/alpha)*\
             - np.cos(firstLambda*alpha)*np.cos((firstLambda - 2)*phi/alpha))
    


def cart2pol(x, y):
    '''This function takes (x,y) coordinates and return the
    corresponding polar coordinates (r, theta)'''
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(r, phi)

def pol2cart(rho, phi):
    '''This function takes (r,theta) coordinates and return the
    corresponding cartesian coordinates (x,y)'''
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

# Initialize the grids
N = 200
Nx = N
Ny = N
dim = (Ny, Nx)
gridStreamFunction = np.zeros(dim) # Grid for the correct solution
gridAlternateStreamFunction = np.zeros(dim) # Grid for an alternate solution
gridFormerSolution = np.zeros(dim) # Grid for the Laplace case

# Save the data in the RESLT folder
data = open("RESLT/python_moffat_singular_function.dat", "w" )
data.write('ZONE I=' + str(N) + ', J=' + str(N) +'\n') # Header of .dat

for i in range(Ny):
    for j in range(Nx):
        x = (j/N - 1/2)
        y = -(i/N - 1/2)
        if i < Ny/2 or j >= Nx/2:
            (rho, phi) = cart2pol(x, y)
            gridStreamFunction[i, j] =  streamFunction(rho, phi - np.pi/4)
            gridAlternateStreamFunction[i, j] = streamFunction2(rho, phi)
            gridFormerSolution[i, j] =  (rho**(2/3))*np.sin(2/3*(phi + np.pi/2))
            u_x = np.cos(phi)*u_r(rho, phi - np.pi/4) - np.sin(phi)*u_theta(rho, phi - np.pi/4)
            u_y = np.sin(phi)*u_r(rho, phi - np.pi/4) + np.cos(phi)*u_theta(rho, phi - np.pi/4)
        line = str(x) + ' ' + str(y) + ' ' + str(gridStreamFunction[i, j]) + ' ' + str(u_x) + ' ' + str(u_y)
        data.write(line + '\n')
 
data.close()
   
#fig = plt.figure()
#plt.imshow(grid)
##plt.imshow(np.log(grid))
#plt.colorbar()
#plt.title("Norm of velocity ($\lambda$ = {0:5.3f})".format(firstLambda))
#plt.xlabel('$x$ coordinate')
#plt.ylabel('$y$ coordinate')
#
## On sauvegarde nos données
##plt.savefig('singular_velocity_plot/v_plot_lambda_is_{0:5.3f}.png'.format(firstLambda), dpi=300)
#
#
### On calcule le gradient de la solution sur la grille
#
#fig = plt.figure()
#grid_grad = np.gradient(grid, axis=1)
#plt.imshow(grid_grad)
#plt.colorbar()
#plt.title("Gradient of vector $v$ ($\lambda$ = {0:5.3f})".format(firstLambda))
#plt.xlabel('$x$ coordinate')
#plt.ylabel('$y$ coordinate')


## Plot every grid
#fig = plt.figure()
#plt.imshow(gridStreamFunction, cmap = plt.cm.binary)
#plt.colorbar()
#plt.title("$\psi$ values ($\lambda$ = {0:5.3f})".format(firstLambda))
#plt.xlabel('$x$ coordinate')
#plt.ylabel('$y$ coordinate')
#
#
#fig = plt.figure()
#plt.contour(np.flipud(gridStreamFunction), 20)
#plt.colorbar()
#plt.title("Contour of singular $\psi$ values ($\lambda$ = {0:5.3f})".format(firstLambda))
#plt.xlabel('$x$ coordinate')
#plt.ylabel('$y$ coordinate')
#plt.gca().set_aspect('equal', adjustable='box')

#plt.savefig('stream_function_plot_lambda_is_{0:5.3f}.eps'.format(firstLambda), format='eps', dpi=1000)
## If you want to remove the frame of the plot:
#plt.savefig('stream_function_plot_lambda_is_{0:5.3f}.eps'.format(firstLambda), format='eps', dpi=1000, bbox_inches='tight', pad_inches=0)



# 2 subplots of the stream function and its contour
fig = plt.figure(figsize=(12,4))


plt.subplot(1, 2, 1)
plt.imshow(gridStreamFunction, cmap=plt.cm.binary,
           interpolation='nearest', extent=[0.5, -0.5, -0.5, 0.5])
plt.colorbar()
plt.title("$\psi$ values ($\lambda$ = {0:5.3f})".format(firstLambda),
          fontsize=18)
plt.xlabel('$x$ coordinate', fontsize=16)
plt.ylabel('$y$ coordinate', fontsize=16)

plt.subplot(1, 2, 2)
plt.contour(np.flipud(gridStreamFunction), 20,
            interpolation='nearest', extent=[-0.5, 0.5, -0.5, 0.5])
plt.colorbar()
plt.grid(True)
plt.title("Contour of $\psi$", fontsize=18)
plt.xlabel('$x$ coordinate', fontsize=16)
plt.ylabel('$y$ coordinate', fontsize=16)
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig('figures/stream_function_plot.eps', format='eps',
            dpi=1000, bbox_inches='tight', pad_inches=0.03)
plt.show()



# Plot of the alternate streamfunction in free surface case
fig = plt.figure()
plt.contour(np.flipud(gridAlternateStreamFunction), 20)
plt.colorbar()
plt.grid(True)
plt.title("Alternate singular function", fontsize=18)
plt.xlabel('$x$ coordinate', fontsize=16)
plt.ylabel('$y$ coordinate', fontsize=16)
plt.gca().set_aspect('equal', adjustable='box')


