############################################################################################
#
# Plotting utilities for cylindrical 3D output from CUMC3D
# Written by Ho Sang (Leon) Chan at 2023
# The Chinese University of Hong Kong and University of Colorado
#
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
def plot_contour(prim,filename,domain='xzhalf',scale='log',stream=False,bh=False):

    # set figure #
    fig, ax = plt.subplots(1, 1)

    # set axis label #
    if(domain == 'xyplane'):
        plt.xlabel('x-direction', size=15)
        plt.ylabel('y-direction', size=15)
    else:
        plt.xlabel('x-direction', size=15)
        plt.ylabel('z-direction', size=15)

    # set x,y axis #
    if(domain == 'xzhalf'):
        x_plot = X_half
        y_plot = Z_half
    elif(domain == 'xzfull'):
        x_plot = X_full
        y_plot = Z_full
    elif(domain == 'xyplane'):
        x_plot = X_plane
        y_plot = Y_plane

    #log scale plot#
    if(scale == 'log'):
        zmax = int(math.ceil(np.log10(np.max(prim))))
        zmin = int(math.floor(np.log10(np.min(prim))))
        cp = ax.contourf(x_plot, y_plot, prim, np.logspace(zmin, zmax, 100), locator=ticker.LogLocator(), cmap='plasma',extend='both')
        cbar = fig.colorbar(cp)
        rang = np.arange(int(np.log10(np.min(prim))), int(np.log10(np.max(prim))) + 1, 1)
        loca = 10 ** (np.array(rang).astype(float))
        cbar.set_ticks(loca)
        cbar.minorticks_off()
        labels = ['10$^{%.0f}$' % x for x in rang]
        cbar.ax.set_yticklabels(labels, fontsize=15)

    # linear scale plot #
    elif(scale == 'linear'):
        cp = ax.contourf(x_plot, y_plot, prim, 100, cmap='plasma', extend='both')
        cbar = fig.colorbar(cp)
        cbar.ax.tick_params(labelsize=15)

    # want stream lines? #
    if(stream):
        if(domain == 'xzhalf'):
            b1 = angular_average(bx)
            b2 = angular_average(bz)
            plt.streamplot(x_plot, y_plot, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)

    # draw a black hole #
    if(bh):
        circle1 = plt.Circle((0, 0), np.min(R), color='black')
        ax.add_patch(circle1)

    # other plot settings #
    plt.xlim(np.min(x_plot), np.max(x_plot))
    plt.ylim(np.min(y_plot), np.max(y_plot))
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
def planar_view(prim):
    p_temp = np.ndarray(shape=(2*prim.shape[0],prim.shape[2]), dtype=float)
    p_temp[0:int(p_temp.shape[0]/2),:] = prim[::-1,0,:]
    p_temp[int(p_temp.shape[0]/2):p_temp.shape[0],:] = prim[:,int(prim.shape[1]/2),:]
    p_temp = p_temp.T
    if(p_temp.shape != X_full.shape):
        p_temp = p_temp.T
    return p_temp

# array for birdeye view #
def birdeye_view(prim):
    p_temp = np.ndarray(shape=(prim.shape[0],prim.shape[1]+1), dtype=float)
    p_temp[:,0:prim.shape[1]] = prim[:,:,int(prim.shape[2]/2)]
    p_temp[:,prim.shape[1]] = prim[:,0,int(prim.shape[2]/2)]
    p_temp = p_temp.T
    if(p_temp.shape != X_plane.shape):
        p_temp = p_temp.T
    return p_temp

# array for angular average #
def angular_average(prim):
    p_temp = np.average(prim, axis=1)
    p_temp = p_temp.T
    if(p_temp.shape != X_half.shape):
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
    # get angular averaged variables #
    z = angular_average(z_in)

    # plot #
    plot_contour(z,fname,domain='xzhalf',scale=scale,stream=stream,bh=False)

    ######################################################################################
    # get planar variables #
    z = planar_view(z_in)

    # plot #
    plot_contour(z,fname,domain='xzfull',scale=scale,stream=False,bh=True)

    ######################################################################################
    # get birdeye variables #
    z = birdeye_view(z_in)

    # plot #
    plot_contour(z,fname,domain='xyplane',scale=scale,stream=False,bh=True)

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

#mesh grid, half x-z plane#
X_half, Z_half = np.meshgrid(xaxis, zaxis)

#mesh grid, x-y plane#
ytemp = np.zeros(yaxis.shape[0] + 1)
ytemp[0:yaxis.shape[0]] = yaxis[:]
ytemp[yaxis.shape[0]] = 2*math.pi + yaxis[0]
R, phi = np.meshgrid(xaxis, ytemp)
X_plane, Y_plane = R*np.cos(phi), R*np.sin(phi)
ytemp = []

#mesh grid, full x-z plane#
xtemp = np.zeros(2*xaxis.shape[0])
xtemp[0:int(xtemp.shape[0]/2)] = -xaxis[::-1]
xtemp[int(xtemp.shape[0]/2):xtemp.shape[0]] = xaxis[:]
X_full, Z_full = np.meshgrid(xtemp, zaxis)
xtemp = []

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
# main plotting functions #

# density #
z = rho
plot(z,'rho',stream=True)

# inverse plasma beta #
z = 2*p/(bx**2 + by**2 + bz**2 + 1e-10)
plot(z,'beta',stream=False)

########################################################################################