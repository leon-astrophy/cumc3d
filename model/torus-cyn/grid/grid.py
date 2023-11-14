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
xr = 21.5

# dr/(r dphi), should not be too small #
ratio_x = 0.4

# number of phi-grid #
nphi = 128

##########################################################

# define function #
def log_grid (xl, xr, nx, nphi):
  left = nx*np.log((nphi + ratio_x*math.pi)/(nphi - ratio_x*math.pi))
  right = np.log(xr/xl)
  return left - right

#assign#
nx_left = 1
nx_right = 1024
f_left = log_grid (xl, xr, nx_left, nphi)
f_right = log_grid (xl, xr, nx_right, nphi)

# bisection #
for n in range (0, 100):
  nx_center = 0.5*(nx_left+nx_right)
  f_center = log_grid (xl, xr, nx_center, nphi)
  if(f_center*f_left > 0):
    nx_left = nx_center
  elif(f_center*f_right > 0):
    nx_right = nx_center

# get value #
nx_center = int(nx_center)
if(nx_center % 2 != 0): nx_center = nx_center + 1

# z has same grid number as x #
nz_center = nx_center

print('nr = '+str(nx_center))
print('nphi = '+str(nphi))
print('nz = '+str(nx_center))

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
# get z grid #

# first, find dx #
dr = np.zeros(nx_center+6)

# assign #
for i in range (0, len(dr)):
  dr[i] = r_face[i+1] - r_face[i]

# now get z-grid #
z_face = np.zeros(nz_center+7)

# assign #
for i in range (3+nz_center//2+1, 3+nz_center+4):
  i_temp = i - (3+nz_center//2)
  z_face[i] = z_face[3+nz_center//2] + np.sum(dr[3:3+i_temp])

# reverse #
for i in range (3+nz_center//2-1, -1, -1):
  i_temp = (3+nz_center//2-1-i) + 3+nz_center//2 + 1
  z_face[i] = -z_face[i_temp]

##########################################################
# phi-grid, trivial #

# dphi #
dphi = 2*math.pi/nphi

# array #
phi_face = np.zeros(nphi+7)

# loop over #
for i in range (0, len(phi_face)):
  phi_face[i] = dphi*(i-3)

##########################################################
# output #

r_face = pd.DataFrame(r_face)
r_face.to_csv('r_grid.dat', index=False, header=False, sep=' ')

th_face = pd.DataFrame(z_face)
th_face.to_csv('z_grid.dat', index=False, header=False, sep=' ')

phi_face = pd.DataFrame(phi_face)
phi_face.to_csv('phi_grid.dat', index=False, header=False, sep=' ')

##########################################################



