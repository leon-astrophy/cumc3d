############################################################################################
#
# Plotting utilities for cartesian 3D output from CUMC3D
# Written by Ho Sang (Leon) Chan at 2023
# The Chinese University of Hong Kong and University of Colorado
#
############################################################################################

#import packages#
import sys
import h5py
import math
import function
import numpy as np

############################################################################################
# define function #

#plotting function #
def plot(z_in,var,vec='vel',stream=False):

  ##########################################################################################
  # set primitive variables #
	prim1 = z_in[:,:,0].T

	# set axis #
	x1 = X_xy
	y1 = Y_xy

	# set vector #
	if(stream):
		if(vec == 'vel'):
			vx1 = velx[:,:,0].T
			vy1 = vely[:,:,0].T
		elif(vec == 'bfield'):
			vx1 = bx[:,:,0].T
			vy1 = by[:,:,0].T
	else:
		vx1 = None
		vy1 = None
	
	fname = imgdir+var+'-xyfull-'+str(time) 
	function.contour(x1,y1,prim1,vx1,vy1,None,None,fname,time)

	##########################################################################################
	#plot#
	if(z_in.shape[2] > 1): 
    
		########################################################################################
		# set primitive variables #
		prim1 = z_in[:,0,:].T

		# set axis #
		x1 = X_xz
		y1 = Z_xz

		# set vector #
		if(stream):
			if(vec == 'vel'):
				vx1 = velx[:,0,:].T
				vy1 = velz[:,0,:].T
			elif(vec == 'bfield'):
				vx1 = bx[:,0,:].T
				vy1 = bz[:,0,:].T
		else:
			vx1 = None
			vy1 = None
		
		fname = imgdir+var+'-xzfull-'+str(time) 
		function.contour(x1,y1,prim1,vx1,vy1,None,None,fname,time)

		########################################################################################
		# set primitive variables #
		prim1 = z_in[0,:,:].T

		# set axis #
		x1 = Y_yz
		y1 = Z_yz

		# set vector #
		if(stream):
			if(vec == 'vel'):
				vx1 = vely[0,:,:].T
				vy1 = velz[0,:,:].T
			elif(vec == 'bfield'):
				vx1 = by[0,:,:].T
				vy1 = bz[0,:,:].T
		else:
			vx1 = None
			vy1 = None
		
		fname = imgdir+var+'-yzfull-'+str(time) 
		function.contour(x1,y1,prim1,vx1,vy1,None,None,fname,time)

############################################################################################
#load command line parameters #

# read path #
hdf5file=sys.argv[1]

# get path #
imgdir = './figure/'

# load hdf5 file #
f = h5py.File(hdf5file, 'r')

############################################################################################
# file input output #

#load grid#
dset = f['x-interface']
xface = dset[:]
dset = f['y-interface']
yface = dset[:]
dset = f['z-interface']
zface = dset[:]

#get cell centered coordiante#
xaxis = np.zeros(len(xface)-3)
yaxis = np.zeros(len(yface)-3)
zaxis = np.zeros(len(zface)-3)

# compute #
for i in range (0, len(xaxis)):
  xaxis[i] = 0.5*(xface[i+2] + xface[i+1])
for i in range (0, len(yaxis)):
  yaxis[i] = 0.5*(yface[i+2] + yface[i+1])
for i in range (0, len(zaxis)):
  zaxis[i] = 0.5*(zface[i+2] + zface[i+1])

############################################################################################
#gridding#

#mesh grid, half x-y plane#
X_xy, Y_xy = np.meshgrid(xaxis, yaxis)

#mesh grid, half x-z plane#
X_xz, Z_xz = np.meshgrid(xaxis, zaxis)

#mesh grid, half y-z plane#
Y_yz, Z_yz = np.meshgrid(yaxis, zaxis)

#check shape#
if(X_xy.shape[0] != xaxis.shape[0]): X_xy = X_xy.T
if(Y_xy.shape[0] != xaxis.shape[0]): Y_xy = Y_xy.T
if(X_xz.shape[0] != xaxis.shape[0]): X_xz = X_xz.T
if(Z_xz.shape[0] != xaxis.shape[0]): Z_xz = Z_xz.T
if(Y_yz.shape[0] != yaxis.shape[0]): Y_yz = Y_yz.T
if(Z_yz.shape[0] != yaxis.shape[0]): Z_yz = Z_yz.T

########################################################################################
# main plotting functions #

# load primitive variables #
dset = f['primitive']
primitive = dset[:]
primitive = primitive.T

# get primitive variables #
rho = primitive[0,1:-1,1:-1,1:-1]
velx = primitive[1,1:-1,1:-1,1:-1]
vely = primitive[2,1:-1,1:-1,1:-1]
velz = primitive[3,1:-1,1:-1,1:-1]
p = primitive[4,1:-1,1:-1,1:-1]

# load magnetic field #
dset = f['bfield']
bfield = dset[:]
bfield = bfield.T

#magnetic fields#
bfacex = bfield[0,:,:,:]
bfacey = bfield[1,:,:,:]
bfacez = bfield[2,:,:,:]

#get cell centered coordiante#
bx = rho.copy()
by = rho.copy()
bz = rho.copy()

# compute #
for i in range (0, len(xaxis)):
  bx[i,:,:] = 0.5*(bfacex[i+2,2:-1,2:-1] + bfacex[i+1,2:-1,2:-1])
for i in range (0, len(yaxis)):
  by[:,i,:] = 0.5*(bfacey[2:-1,i+2,2:-1] + bfacey[2:-1,i+1,2:-1])
for i in range (0, len(zaxis)):
  bz[:,:,i] = 0.5*(bfacez[2:-1,2:-1,i+1] + bfacez[2:-1,2:-1,i+1])

#load epsilon#
dset = f['epsilon']
epsilon = dset[:]
epsilon = epsilon.T
epsilon = epsilon[1:-1,1:-1,1:-1]

#time#
dset = f['time']
time = dset[:][0]
time = np.round(time, 3)

########################################################################################
#plot#

# density #
z = rho
plot(z,'rho',vec='vel',stream=True)

# plasma beta #
z = 2*p/(bx**2 + by**2 + bz**2 + 1e-10)
plot(z,'beta',vec='bfield',stream=True)

########################################################################################