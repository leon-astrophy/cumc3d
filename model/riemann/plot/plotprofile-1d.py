#import#
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt

############################################################################################
# filename#
filename = []

#load#
for root, dirs, files in os.walk('../outfile/'):
    for file in files:
        if file.endswith("nm.hdf5"):
            filename.append(os.path.join(file))

#load#
f = h5py.File('../outfile/grid_param.hdf5', 'r')
dset = f['x-direction']
xaxis = dset[:]

############################################################################################
#loop#
for i in range (0, len(filename)):

    #load#
    f = h5py.File('../outfile/'+filename[i], 'r')

    #data set#
    dset = f['primitive']
    primitive=dset[:]
    primitive = primitive.T

    #assign#
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

    #time#
    dset = f['time']
    time = dset[:][0]

    ########################################################################################
    # plot#
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
    plt.savefig('velprofile-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

    ########################################################################################
    # plot#
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
    plt.savefig('bfieldprofile-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

    ########################################################################################
    # plot#
    fig, axs = plt.subplots(2, 1, sharex='all', gridspec_kw={'hspace': 0, 'wspace': 0})
    ((ax1, ax2)) = axs
    ax1.plot(xaxis, rho)
    ax2.plot(xaxis, p, color='tab:red')
    ax1.set_title('Time = '+'%.3f (Code Unit)' % time, size=15)
    ax1.set_ylabel('Density (Code Unit)', size=15)
    ax2.set_ylabel('Pressure', size=15)
    ax2.set_xlabel('x-direction (Code Unit)', size=15)
    ax1.grid()
    ax2.grid()
    plt.savefig('prhoprofile-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

