############################################################################################
# Plotting utilities for 1D output from CUMC3D
############################################################################################

#import#
import os
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

############################################################################################
#load command line parameters #

path=sys.argv[1]

############################################################################################
# file input output #

# get path #
imgdir = './figure/'

# initialise filename #
filename = []

# get all filename #
for root, dirs, files in os.walk(path):
    for file in files:
        if file.endswith("nm.hdf5"):
            filename.append(os.path.join(file))

# also get grid parameters #
f = h5py.File(path+'./outfile/grid_param.hdf5', 'r')
dset = f['x-direction']
xaxis = dset[:]

############################################################################################
# plotting here #

# loop over all files #
for i in range (0, len(filename)):

    # load hdf5 #
    f = h5py.File(path+'outfile/'+filename[i], 'r')

    # load primitive variables #
    dset = f['primitive']
    primitive = dset[:]
    primitive = primitive.T

    # get primitive variables #
    rho = primitive[0,:,0,0]
    velx = primitive[1,:,0,0]
    vely = primitive[2,:,0,0]
    velz = primitive[3,:,0,0]
    p = primitive[4,:,0,0]

    #magnetic fields#
    neq = primitive.shape[0]
    bx = primitive[neq-3,:,0,0]
    by = primitive[neq-2,:,0,0]
    bz = primitive[neq-1,:,0,0]

    #load epsilon#
    dset = f['epsilon']
    epsilon = dset[:]
    epsilon = epsilon.T
    epsilon = epsilon[:,0,0]

    #time#
    dset = f['time']
    time = dset[:][0]

    ########################################################################################
    # plot all profiles #

    # x, y, z velocity #
    fig, axs = plt.subplots(3, 1, sharex='all', gridspec_kw={'hspace': 0, 'wspace': 0})
    ((ax1, ax2, ax3)) = axs
    ax1.plot(xaxis, velx, label='x-velocity')
    ax2.plot(xaxis, vely, label='y-velocity', color='tab:red')
    ax3.plot(xaxis, velz, label='z-velocity', color='tab:purple')
    ax1.set_title('Time = '+'%.3f (Code Unit)' % time, size=15)
    ax3.set_xlabel('x-direction (Code Unit)', size=15)
    ax2.set_ylabel('Velocity (Code Unit)', size=15)
    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax3.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.savefig(imgdir+'vel-profile-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

    ########################################################################################
    
    # x, y, z magnetic fields #
    fig, axs = plt.subplots(3, 1, sharex='all', gridspec_kw={'hspace': 0, 'wspace': 0})
    ((ax1, ax2, ax3)) = axs
    ax1.plot(xaxis, bx, label='x-bfield')
    ax2.plot(xaxis, by, label='y-bfield', color='tab:red')
    ax3.plot(xaxis, bz, label='z-bfield', color='tab:purple')
    ax1.set_title('Time = '+'%.3f (Code Unit)' % time, size=15)
    ax3.set_xlabel('x-direction (Code Unit)', size=15)
    ax2.set_ylabel('Magnetic Fields (Code Unit)', size=15)
    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax3.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.savefig(imgdir+'bfield-profile-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

    ########################################################################################
    
    # density, epsilon, pressure #
    fig, axs = plt.subplots(3, 1, sharex='all', gridspec_kw={'hspace': 0, 'wspace': 0})
    ((ax1, ax2, ax3)) = axs
    ax1.plot(xaxis, rho, label='Density')
    ax2.plot(xaxis, p, label='Pressure', color='tab:red')
    ax3.plot(xaxis, epsilon, label='Epsilon', color='tab:purple')
    ax1.set_title('Time = '+'%.3f (Code Unit)' % time, size=15)
    ax3.set_xlabel('x-direction (Code Unit)', size=15)
    ax2.set_ylabel('(Code Units)', size=15)
    ax1.grid()
    ax2.grid()
    ax3.grid()
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax1.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    ax3.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.savefig(imgdir+'rhopeps-profile-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

    ########################################################################################

