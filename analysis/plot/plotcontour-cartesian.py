############################################################################################
# Plotting utilities for cartesian 3D output from CUMC3D
############################################################################################

#import packages#
import os
import sys
import h5py
import math
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker, cm, colors
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import FormatStrFormatter

############################################################################################
# define function #

# plot contour #
def plot_contour(X1_in,X2_in,prim,filename,domain='xy',scale='log',stream=False,bh=False):

    # set figure #
    fig, ax = plt.subplots(1, 1)

    # set axis ticks #
    if(domain == 'xy'):
        plt.xlabel('x-direction', size=15)
        plt.ylabel('y-direction', size=15)
    elif(domain == 'xz'):
        plt.xlabel('x-direction', size=15)
        plt.ylabel('z-direction', size=15)
    elif(domain == 'yz'):
        plt.xlabel('y-direction', size=15)
        plt.ylabel('z-direction', size=15)

    # log scale plot #
    if(scale == 'log'):
        zmax = int(math.ceil(np.log10(np.max(prim))))
        zmin = int(math.floor(np.log10(np.min(prim))))
        cp = ax.contourf(X1_in, X2_in, prim, np.logspace(zmin, zmax, 100), locator=ticker.LogLocator(), cmap='RdGy',extend='both')
        cbar = fig.colorbar(cp)
        rang = np.arange(int(np.log10(np.min(prim))), int(np.log10(np.max(prim))) + 1, 1)
        loca = 10 ** (np.array(rang).astype(float))
        cbar.set_ticks(loca)
        cbar.minorticks_off()
        labels = ['10$^{%.0f}$' % x for x in rang]
        cbar.ax.set_yticklabels(labels, fontsize=15)

    # linear scale plot #
    elif(scale == 'linear'):
        cp = ax.contourf(X1_in, X2_in, prim, 100, cmap='RdGy', extend='both')
        cbar = fig.colorbar(cp)
        cbar.ax.tick_params(labelsize=15)
    
    # plot bfield stream line ?#
    if(stream):
        if(domain == 'xy'):
            b1 = xy_full(bx)
            b2 = xy_full(by)
        elif(domain == 'xz'):
            b1 = xz_full(bx)
            b2 = xz_full(bz)
        elif(domain == 'yz'):
            b1 = yz_full(by)
            b2 = yz_full(bz)
        plt.streamplot(X1_in, X2_in, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)

    #plot a black hole at the center #
    if(bh):
        circle1 = plt.Circle((0, 0), np.min(R), color='black')
        ax.add_patch(circle1) 

    # axis properties # 
    plt.xlim(np.min(X1_in), np.max(X1_in))
    plt.ylim(np.min(X2_in), np.max(X2_in))
    ax.tick_params(axis='both', labelsize=15)
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(imgdir+str(filename)+'-'+str(domain)+'-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

############################################################################################
# define function #

# array for planar view #
def xy_full(prim):
    p_temp = np.ndarray(shape=(prim.shape[0],prim.shape[1]), dtype=float)
    p_temp[:,:] = prim[:,:,0]
    p_temp = p_temp.T
    if(p_temp.shape != X_xy.shape):
      p_temp = p_temp.T
    return p_temp

# array for birdeye view #
def xz_full(prim):
    p_temp = np.ndarray(shape=(prim.shape[0],prim.shape[2]), dtype=float)
    p_temp[:,:] = prim[:,0,:]
    p_temp = p_temp.T
    if(p_temp.shape != X_xz.shape):
      p_temp = p_temp.T
    return p_temp

# array for angular average #
def yz_full(prim):
    p_temp = np.ndarray(shape=(prim.shape[1],prim.shape[2]), dtype=float)
    p_temp[:,:] = prim[0,:,:]
    p_temp = p_temp.T
    if(p_temp.shape != Y_yz.shape):
      p_temp = p_temp.T
    return p_temp

############################################################################################
# define function #

#plotting function #
def plot(z_in,fname,stream=False):

    # choose plotting scale #
    if((np.min(z_in) == np.max(z_in)) or (np.min(z_in) <= 0) or (np.max(z_in) <= 0) or np.isnan(z_in).any()):
        scale = 'linear'
    else:
        scale = 'log'
    
    ######################################################################################
    # x-y projection variables #
    z = xy_full(z_in)

    #plot#
    plot_contour(X_xy,Y_xy,z,fname,domain='xy',scale=scale,stream=stream,bh=False)

    ######################################################################################
    # x-z projection variables #
    z = xz_full(z_in)

    #plot#
    try:
        plot_contour(X_xz,Z_xz,z,fname,domain='xz',scale=scale,stream=stream,bh=False)
    except:
        pass

    ######################################################################################
    # y-z projection variables #
    z = yz_full(z_in)

    #plot#
    try:
        plot_contour(Y_yz,Z_yz,z,fname,domain='yz',scale=scale,stream=stream,bh=False)
    except:
        pass

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

############################################################################################
#gridding#

#mesh grid, half x-y plane#
X_xy, Y_xy = np.meshgrid(xaxis, yaxis)

#mesh grid, half x-z plane#
X_xz, Z_xz = np.meshgrid(xaxis, zaxis)

#mesh grid, half y-z plane#
Y_yz, Z_yz = np.meshgrid(yaxis, zaxis)

########################################################################################
# main plotting functions #

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

#magnetic fields#
neq = primitive.shape[0]
bx = primitive[neq-3,:,:,:]
by = primitive[neq-2,:,:,:]
bz = primitive[neq-1,:,:,:]

#load epsilon#
dset = f['epsilon']
epsilon = dset[:]
epsilon = epsilon.T
epsilon = epsilon[:,:,:]

#time#
dset = f['time']
time = dset[:][0]

########################################################################################
#plot#

# density #
z = rho
plot(z,'rho',stream=True)

# plasma beta #
z = 2*p/(bx**2 + by**2 + bz**2 + 1e-10) + 1e-10
plot(z,'beta',stream=False)

########################################################################################