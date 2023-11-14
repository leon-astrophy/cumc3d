##########################################################
#
# Solving for grid variables in spherical simulations
#
##########################################################

#import#
import math
import numpy as np
import pandas as pd

##########################################################

# inner radial boundary #
xl = 1.5

# outer radial boundary #
xr = 31.5

# dr/(r dphi), should not be too small #
ratio_x = 0.5

# (r sin(theta)dphi)/(dr) at pole #
ratio_y = 0.05

# number of z-grid #
nz = 128

##########################################################

# solve for alpha #
alpha_y = np.arcsin(ratio_y)*nz/math.pi

##########################################################

# define function #
def strech_y (alpha, ny, nz):
  n_grid = float(ny/2)
  r_ratio = alpha**(1/(n_grid-1))
  return (r_ratio**n_grid - 1)/(r_ratio - 1) - float(nz/4)

#assign#
ny_left = 1
ny_right = 256
f_left = strech_y (alpha_y, ny_left, nz)
f_right = strech_y (alpha_y, ny_right, nz)

# bisection #
for n in range (0, 100):
  ny_center = 0.5*(ny_left+ny_right)
  f_center = strech_y (alpha_y, ny_center, nz)
  if(f_center*f_left > 0):
    ny_left = ny_center
  elif(f_center*f_right > 0):
    ny_right = ny_center

# get value #
ny_center = int(ny_center)
if(ny_center % 2 != 0): ny_center = ny_center + 1

# define function #
def log_grid (xl, xr, nx, nz):
  left = nx*np.log((nz + ratio_x*math.pi)/(nz - ratio_x*math.pi))
  right = np.log(xr/xl)
  return left - right

#assign#
nx_left = 1
nx_right = 256
f_left = log_grid (xl, xr, nx_left, nz)
f_right = log_grid (xl, xr, nx_right, nz)

# bisection #
for n in range (0, 100):
  nx_center = 0.5*(nx_left+nx_right)
  f_center = log_grid (xl, xr, nx_center, nz)
  if(f_center*f_left > 0):
    nx_left = nx_center
  elif(f_center*f_right > 0):
    nx_right = nx_center

# get value #
nx_center = int(nx_center)
if(nx_center % 2 != 0): nx_center = nx_center + 1

print('nx = '+str(nx_center))
print('ny = '+str(ny_center))
print('nz = '+str(nz))

##########################################################
# get r grid #

# solve for epsilon #
deps = np.log(xr/xl)/nx_center

# array #
r_face = np.zeros(nx_center+7)

# assign ghost shell#
dr = xl/3
for i in range (0, 4):
  r_face[i] = i*dr

# assign active cells #
for i in range (4, len(r_face)):
  r_face[i] = r_face[i-1]*np.exp(deps)

##########################################################
# get theta grid #

# solve for the strecting ratio #
r_left = 1.01
r_right = 3.01

# define function #
def strech_solve (r_in, alpha, ny, nz):
  n_grid = float(ny/2)
  return (r_in**n_grid - 1)/(r_in - 1) - float(nz/4)

# assign left and right #
f_left = strech_solve(r_left, alpha_y, ny_center, nz)
f_right = strech_solve(r_right, alpha_y, ny_center, nz)

# bisection #
for n in range (0, 100):
  r_center = 0.5*(r_left+r_right)
  f_center = strech_solve(r_center, alpha_y, ny_center, nz)
  if(f_center*f_left > 0):
    r_left = r_center
  elif(f_center*f_right > 0):
    r_right = r_center

# array #
th_face = np.zeros(ny_center+7)

# dphi at equator #
dphi = 2*math.pi/nz

# loop over theta > pi/2 #
th_face[3 + ny_center//2] = math.pi/2
for i in range (3+ny_center//2 + 1, 3+ny_center+1):
  th_face[i] = th_face[i-1] + dphi*r_center**(i - 3 - ny_center//2 - 1)
th_face[3+ny_center] = math.pi 

# ghost shell 
for i in range (3+ny_center+1, 3+ny_center+4):
  th_face[i] = th_face[i-1] + (th_face[2*(3+ny_center) - i + 1] - th_face[2*(3+ny_center) - i])

# loop over the theta < pi/2 #
for i in range (3+ny_center//2 - 1, 2, -1):
  th_face[i] = math.pi - th_face[6+ny_center - i]

# ghost shell #
for i in range (2, -1, -1):
  th_face[i] = - th_face[6-i]

##########################################################
# z-grid, trivial #

# array #
phi_face = np.zeros(nz+7)

# loop over #
for i in range (0, len(phi_face)):
  phi_face[i] = dphi*(i-3)

##########################################################
# output #

r_face = pd.DataFrame(r_face)
r_face.to_csv('r_grid.dat', index=False, header=False, sep=' ')

th_face = pd.DataFrame(th_face)
th_face.to_csv('th_grid.dat', index=False, header=False, sep=' ')

phi_face = pd.DataFrame(phi_face)
phi_face.to_csv('phi_grid.dat', index=False, header=False, sep=' ')

##########################################################



