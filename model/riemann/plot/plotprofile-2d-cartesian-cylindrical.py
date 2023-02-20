#import packages#
import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors

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
dset = f['y-direction']
yaxis = dset[:]

############################################################################################
#loop#
for i in range (0, len(filename)):

    #load#
    f = h5py.File('../outfile/'+filename[i], 'r')
    dset = f['primitive']
    primitive=dset[:]
    primitive = primitive.T
    neq = primitive.shape[0]

    # assign #
    rho = primitive[0,:,:,0].T
    velx = primitive[1,:,:,0].T
    vely = primitive[2, :, :, 0].T
    p = primitive[4,:,:,0].T
    bx = primitive[neq-3,:,:,0].T
    by = primitive[neq-2,:,:,0].T

    #assign#
    dset = f['time']
    time = dset[:][0]

    ########################################################################################
    #plot#
    X, Y = np.meshgrid(xaxis, yaxis)
    fig, ax = plt.subplots(1, 1)
    cp = ax.contourf(X, Y, rho, 1000, cmap='coolwarm')
    cbar = fig.colorbar(cp)
    plt.streamplot(X,Y,velx,vely, density=1.4, linewidth=0.5, color='white')
    plt.title('Time = ' + '%.3f' % time)
    plt.xlabel('x-direction')
    plt.ylabel('y-direction')
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(xaxis), np.max(xaxis))
    plt.ylim(np.min(yaxis), np.max(yaxis))
    plt.savefig('rho-vstream-contour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

    ########################################################################################
    #plot#
    X, Y = np.meshgrid(xaxis, yaxis)
    fig, ax = plt.subplots(1, 1)
    cp = ax.contourf(X, Y, p, 1000, cmap='coolwarm')
    cbar = fig.colorbar(cp)
    plt.streamplot(X,Y,velx,vely, density=1.4, linewidth=0.5, color='white')
    plt.title('Time = ' + '%.3f' % time)
    plt.xlabel('x-direction')
    plt.ylabel('y-direction')
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(xaxis), np.max(xaxis))
    plt.ylim(np.min(yaxis), np.max(yaxis))
    plt.savefig('p-vstream-contour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

    ########################################################################################
    #plot#
    X, Y = np.meshgrid(xaxis, yaxis)
    fig, ax = plt.subplots(1, 1)
    cp = ax.contourf(X, Y, rho, 1000, cmap='coolwarm')
    cbar = fig.colorbar(cp)
    plt.streamplot(X,Y,bx,by, density=1.4, linewidth=0.5, color='white')
    plt.title('Time = ' + '%.3f' % time)
    plt.xlabel('x-direction')
    plt.ylabel('y-direction')
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(xaxis), np.max(xaxis))
    plt.ylim(np.min(yaxis), np.max(yaxis))
    plt.savefig('rho-bstream-contour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

    ########################################################################################
    #plot#
    X, Y = np.meshgrid(xaxis, yaxis)
    fig, ax = plt.subplots(1, 1)
    cp = ax.contourf(X, Y, p, 1000, cmap='coolwarm')
    cbar = fig.colorbar(cp)
    plt.streamplot(X,Y,bx,by, density=1.4, linewidth=0.5, color='white')
    plt.title('Time = ' + '%.3f' % time)
    plt.xlabel('x-direction')
    plt.ylabel('y-direction')
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(xaxis), np.max(xaxis))
    plt.ylim(np.min(yaxis), np.max(yaxis))
    plt.savefig('p-bstream-contour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()