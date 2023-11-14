########################################################################################
#
# Plotting utilities for spherical 3D output from CUMC3D
# Written by Ho Sang (Leon) Chan at 2023
# The Chinese University of Hong Kong and University of Colorado
#
########################################################################################

# os #
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

#import packages#
import sys
import h5py
import function
import numpy as np

########################################################################################
#plotting function #
def plot(z_in,var,vec='vel',stream=False,bh=False):

	######################################################################################
	# black hole radius #

  #black hole radius #
	rbh = np.min(xaxis)

  ######################################################################################
  # set primitive variables #
	prim1 = z_in[:,:,0]
	prim2 = z_in[:,:,mirror]

	# set axis #
	x1 = x_xz
	y1 = z_xz

	# set vector #
	if(stream):
		if(vec == 'vel'):
			vx1 = velx[:,:,0]
			vy1 = vely[:,:,0]
		elif(vec == 'bfield'):
			vx1 = bx[:,:,0]
			vy1 = by[:,:,0]
	else:
		vx1 = None
		vy1 = None

	# set filename and plot #
	fname = imgdir+var+'-xzfull-'+str(time)
	function.contour_interpolate_stream(x1,y1,prim1,xaxis,yaxis,vx1,vy1,prim2,rbh,fname,time)

	######################################################################
	#plot#
	if(z_in.shape[2] > 1): 

		####################################################################################
		# set primitive variables #
		prim1 = z_in[:,z_in.shape[1]//2,:]

		# set axis #
		x1 = x_xy
		y1 = y_xy

		# set vector #
		if(stream):
			if(vec == 'vel'):
				vx1 = velx[:,z_in.shape[1]//2,:]
				vy1 = velz[:,z_in.shape[1]//2,:]
			elif(vec == 'bfield'):
				vx1 = bx[:,z_in.shape[1]//2,:]
				vy1 = bz[:,z_in.shape[1]//2,:]
		else:
			vx1 = None
			vy1 = None

		# set filename and plot #
		fname = imgdir+var+'-xyfull-'+str(time)
		function.contour(x1,y1,prim1,None,None,None,rbh,fname,time)

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

#mesh grid#
R_xz, Theta_xz = np.meshgrid(xaxis, yaxis)

#check shape#
if(R_xz.shape[0] != xaxis.shape[0]): R_xz = R_xz.T
if(Theta_xz.shape[0] != xaxis.shape[0]): Theta_xz = Theta_xz.T

# get x and z #
x_xz = R_xz * np.sin(Theta_xz)
z_xz = R_xz * np.cos(Theta_xz)

#mesh grid#
R_xy, Theta_xy = np.meshgrid(xaxis, zaxis)

#check shape#
if(R_xy.shape[0] != xaxis.shape[0]): R_xy = R_xy.T
if(Theta_xy.shape[0] != xaxis.shape[0]): Theta_xy = Theta_xy.T

# get x and y #
x_xy = R_xy * np.cos(Theta_xy)
y_xy = R_xy * np.sin(Theta_xy)

########################################################################################
# file input output #

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
  bx[i,:,:] = 1.5*(xface[i+2]-xface[i+1])/(xface[i+2]**3-xface[i+1]**3)\
							* (xface[i+2]**2*bfacex[i+2,2:-1,2:-1] + xface[i+1]**2*bfacex[i+1,2:-1,2:-1])
for i in range (0, len(yaxis)):
  by[:,i,:] = 0.5*(yface[i+2]-yface[i+1])\
							* (np.sin(yface[i+2])*bfacey[2:-1,i+2,2:-1] + np.sin(yface[i+1])*bfacey[2:-1,i+1,2:-1])\
							/ np.abs(np.cos(yface[i+2])-np.cos(yface[i+1]))
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

# set mirror #
if(rho.shape[2] > 1):
	mirror = rho.shape[2]//2
else:
	mirror = 0

########################################################################################
#plot#

# want black hole?#
bh = True

# density #
z = np.log10(rho/rho.max())
plot(z,'rho',vec='bfield',stream=False,bh=bh)

# speed of sound #
z = np.sqrt(p/rho)
#plot(z,'cs',stream=False,bh=bh)

# internal energy #
z = epsilon
#plot(z,'epsilon',stream=False,bh=bh)

# pressure
z = p
#plot(z,'pressure',stream=False,bh=bh)

# velocity x 
z = velx
#plot(z,'velx',stream=False,bh=bh)

# velocity y 
z = vely
#plot(z,'vely',stream=False,bh=bh)

# velocity z
z = velz
#plot(z,'velz',stream=False,bh=bh)

# total velocity
z = np.sqrt(velx**2+vely**2+velz**2)
plot(z,'vel',stream=False,bh=bh)

# plasma beta #
z = 2*p/(bx**2 + by**2 + bz**2 + 1e-10)
z = np.log10(z)
plot(z,'beta',stream=False,bh=bh)

# bfield x #
z = bx
#plot(z,'bx',stream=False,bh=bh)

# bfield y #
z = by
#plot(z,'by',stream=False,bh=bh)

# bfield z #
z = bz
#plot(z,'bz',stream=False,bh=bh)

# total bfield #
z = np.sqrt(bx**2 + by**2 + bz**2)
plot(z,'bfield',stream=False,bh=bh)

# alven velocity #
z = np.sqrt(bx**2 + by**2 + bz**2)/np.sqrt(rho)
plot(z,'alven',stream=False,bh=bh)

# total energy density #
z = rho*epsilon + 0.5*rho*(velx**2+vely**2+velz**2) + 0.5*(bx**2 + by**2 + bz**2)
#plot(z,'energy',stream=False,bh=bh)

########################################################################################
