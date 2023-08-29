############################################################################################
#
# Plotting utilities for spherical 3D output from CUMC3D
# Written by Ho Sang (Leon) Chan at 2023
# The Chinese University of Hong Kong and University of Colorado
#
#############################################################################################

# os #
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

#import packages#
import sys
import h5py
import math
import numpy as np
from matplotlib import ticker
import matplotlib.pyplot as plt

############################################################################################
# define function #

# plot contour #
def plot_contour(prim,filename,domain='xzhalf',scale='log',stream=False,bh=False):

	######################################################################
	# set figure #
	fig, ax = plt.subplots(1, 1)

	######################################################################
	# set label #
	if(domain == 'xyplane'):
		plt.xlabel('x-direction', size=15)
		plt.ylabel('y-direction', size=15)
	else:
		plt.xlabel('x-direction', size=15)
		plt.ylabel('z-direction', size=15)

	######################################################################
	# set x,y axis #
	if(domain == 'xzhalf'):
		x_plot = X_half
		y_plot = Z_half
		z_plot = prim[:,:,0].T
		if(z_plot.shape != x_plot.shape): z_plot = z_plot.T
	elif(domain == 'xzfull'):
		x1_plot = X_half
		y1_plot = Z_half
		x2_plot = -X_half
		y2_plot = Z_half
		z1_plot = prim[:,:,0].T
		z2_plot = prim[:,:,int(prim.shape[2]/2)].T
		if(z1_plot.shape != x1_plot.shape): z1_plot = z1_plot.T
		if(z2_plot.shape != x2_plot.shape): z2_plot = z2_plot.T
	elif(domain == 'xyplane'):
		x_plot = X_plane
		y_plot = Y_plane
		z_plot = prim[:,0,:].T
		if(z_plot.shape != x_plot.shape): z_plot = z_plot.T

	######################################################################
	# log plot #
	if(scale == 'log'):
		zmax = int(math.ceil(np.log10(np.max(prim))))
		zmin = int(math.floor(np.log10(np.min(prim))))
		if(domain == 'xzfull'):
			cp = ax.contourf(x1_plot, y1_plot, z1_plot, np.logspace(zmin, zmax, 100), locator=ticker.LogLocator(), cmap='plasma',extend='both')
			cp = ax.contourf(x2_plot, y2_plot, z2_plot, np.logspace(zmin, zmax, 100), locator=ticker.LogLocator(), cmap='plasma',extend='both')
		else:
			cp = ax.contourf(x_plot, y_plot, z_plot, np.logspace(zmin, zmax, 100), locator=ticker.LogLocator(), cmap='plasma',extend='both')
		cbar = fig.colorbar(cp)
		rang = np.arange(int(np.log10(np.min(prim))), int(np.log10(np.max(prim))) + 1, 1)
		loca = 10 ** (np.array(rang).astype(float))
		cbar.set_ticks(loca)
		cbar.minorticks_off()
		labels = ['10$^{%.0f}$' % x for x in rang]
		cbar.ax.set_yticklabels(labels, fontsize=15)

	######################################################################
	# linear plot #
	elif(scale == 'linear'):
		if(domain == 'xzfull'):
			cp = ax.contourf(x1_plot, y1_plot, z1_plot, 100, cmap='plasma', extend='both')
			cp = ax.contourf(x2_plot, y2_plot, z2_plot, 100, cmap='plasma', extend='both')
		else:
			cp = ax.contourf(x_plot, y_plot, z_plot, 100, cmap='plasma', extend='both')
		cbar = fig.colorbar(cp)
		cbar.ax.tick_params(labelsize=15)

	######################################################################
	# want stream lines? #
	if(stream):
		if(domain == 'xzhalf'):
			skip = (slice(None, None, 8), slice(None, None, 4))
			b_scale = np.sqrt(bx_plot**2 + by_plot**2)
			U = bx_plot/b_scale
			V = by_plot/b_scale
			ax.quiver(x_plot[skip], y_plot[skip], U[skip], V[skip], scale=10, width=.01)
		#elif(domain == 'xzfull'): 
		#elif(domain == 'xyplane'): 

	######################################################################
	# draw black hole #
	if(bh):
		circle1 = plt.Circle((0, 0), np.min(R), color='black')
		ax.add_patch(circle1)

	######################################################################
	# set domain #
	if(domain == 'xzhalf'):
		plt.xlim(0,1000)#np.max(x_plot))
		plt.ylim(-1000,1000)#np.min(y_plot), np.max(y_plot)) 
	else:
		plt.xlim(-1000,1000)#np.min(x_plot), np.max(x_plot))
		plt.ylim(-1000,1000)#np.min(y_plot), np.max(y_plot))       

	######################################################################
	# further plot parameters # 
	ax.tick_params(axis='both', labelsize=15)
	plt.title('Time = ' + '%.3f' % time, size=15)
	plt.gca().set_aspect('equal')
	plt.tight_layout()
	plt.savefig(imgdir+str(filename)+'-'+str(domain)+'-'+'%.3f' % time +'.png')
	plt.clf()
	plt.close()

############################################################################################
#plotting function #
def plot(z_in,fname,stream=False):

	######################################################################
	# choose plotting scale #
	if((np.min(z_in) == np.max(z_in)) or (np.min(z_in) <= 0) \
      or (np.max(z_in) <= 0) or (np.isnan(z_in).any())):
		scale = 'linear'
	else:
		scale = 'log'
	
	######################################################################
	# plot #
	plot_contour(z,fname,domain='xzhalf',scale=scale,stream=stream,bh=True)

	######################################################################
	#plot#
	plot_contour(z,fname,domain='xzfull',scale=scale,stream=False,bh=True)

	######################################################################
	#plot#
	if(z.shape[2] > 1): plot_contour(z,fname,domain='xyplane',scale=scale,stream=False,bh=True)

############################################################################################
#load command line parameters #

# read path #
gridfile=sys.argv[1]
hdf5file=sys.argv[2]

# get path #
imgdir = './figure/'

############################################################################################
# file input output #

#load grid#
f = h5py.File(gridfile, 'r')
dset = f['x-direction']
xaxis = dset[:]
dset = f['y-direction']
yaxis = dset[:]
dset = f['z-direction']
zaxis = dset[:]

########################################################################################
# loading variables #

# load hdf5 file #
f = h5py.File(hdf5file, 'r')

# load primitive variables #
dset = f['primitive']
primitive = dset[:]
primitive = primitive.T

# get primitive variables #
rho = primitive[0,:,:,:]
velx = primitive[1,:,:,:]
vely = primitive[2,:,:,:]
velz = primitive[3,:,:,:]
p = primitive[4,:,:,:]

# load magnetic field #
dset = f['bfield']
bfield = dset[:]
bfield = bfield.T

#magnetic fields#
bx = bfield[0,:,:,:]
by = bfield[1,:,:,:]
bz = bfield[2,:,:,:]

#load epsilon#
dset = f['epsilon']
epsilon = dset[:]
epsilon = epsilon.T
epsilon = epsilon[:,:,:]

#time#
dset = f['time']
time = dset[:][0]

############################################################################################

#mesh grid#
R, Theta = np.meshgrid(xaxis, yaxis)
X_half = R * np.sin(Theta)
Z_half = R * np.cos(Theta)

# vector field #
b_total = np.sqrt(bx**2+by**2+bz**2)
bx_plot = bx[:,:,0].T*np.sin(Theta) + by[:,:,0].T*np.cos(Theta)
by_plot = bx[:,:,0].T*np.cos(Theta) - by[:,:,0].T*np.sin(Theta)

#mesh grid#
R, Theta = np.meshgrid(xaxis, zaxis)
X_plane = R * np.cos(Theta)
Y_plane = R * np.sin(Theta)

########################################################################################
#plot#

# density #
z = rho
plot(z,'rho',stream=True)

# speed of sound #
z = np.sqrt(p/rho)
plot(z,'cs',stream=True)

# density #
z = epsilon
#plot(z,'epsilon',stream=False)

#pressure
z = p
#plot(z,'pressure',stream=False)

# velocity x 
z = velx
plot(z,'velx',stream=False)

# velocity y 
z = vely
#plot(z,'vely',stream=False)

# velocity z
z = velz
#plot(z,'velz',stream=False)

# velocity
z = np.sqrt(velx**2+vely**2+velz**2)
#plot(z,'vel',stream=False)

# inverse plasma beta #
z = 2*p/(bx**2 + by**2 + bz**2)
plot(z,'beta',stream=False)

# inverse plasma beta #
z = bx
#plot(z,'bx',stream=False)

# inverse plasma beta #
z = by
#plot(z,'by',stream=False)

# inverse plasma beta #
z = bz
#plot(z,'bz',stream=False)

# inverse plasma beta #
z = np.sqrt(bx**2 + by**2 + bz**2)
#plot(z,'bfield',stream=False)

# alven velocity #
z = np.sqrt(bx**2 + by**2 + bz**2)/np.sqrt(rho)
plot(z,'alven',stream=False)

# density #
z = rho*epsilon + 0.5*rho*(velx**2+vely**2+velz**2) + 0.5*(bx**2 + by**2 + bz**2)
#plot(z,'energy',stream=False)

########################################################################################
