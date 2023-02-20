#import packages#
import os
import h5py
import math
import matplotlib
import numpy as np
import matplotlib.ticker as tick
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker, cm, colors
from matplotlib.ticker import FormatStrFormatter

############################################################################################
# define function #
def linear_plot(prim,filename):
    X, Y = np.meshgrid(xaxis, zaxis)
    fig, ax = plt.subplots(1, 1)
    z = prim
    cp = ax.contourf(X, Y, z, 100, cmap='plasma')
    cbar = fig.colorbar(cp)
    for art in ax.get_children():
        if not isinstance(art, matplotlib.patches.FancyArrowPatch):
            continue
        art.remove()
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.xlabel('x-direction', size=15)
    plt.ylabel('z-direction', size=15)
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(np.min(Y), np.max(Y))
    plt.tight_layout()
    plt.savefig(str(filename)+'-contour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

# define function #
def log_plot(prim,filename,stream=False):
    X, Y = np.meshgrid(xaxis, zaxis)
    fig, ax = plt.subplots(1, 1)
    z = prim
    zmax = int(np.log10(np.max(z)))
    zmin = int(np.log10(np.min(z)))
    cp = ax.contourf(X, Y, z, np.logspace(zmin, zmax, 100), locator=ticker.LogLocator(), cmap='plasma',extend='both')
    cbar = fig.colorbar(cp)
    rang = np.arange(int(np.log10(np.min(z))), int(np.log10(np.max(z))) + 1, 1)
    loca = 10 ** (np.array(rang).astype(float))
    cbar.set_ticks(loca)
    cbar.minorticks_off()
    labels = ['10$^{%.0f}$' % x for x in rang]
    cbar.ax.set_yticklabels(labels, fontsize=15)
    if(stream):
        plt.streamplot(X,Y,bx,bz, density=1.5, linewidth=None, color='black', broken_streamlines=True)
    for art in ax.get_children():
        if not isinstance(art, matplotlib.patches.FancyArrowPatch):
            continue
        art.remove()
    ax.tick_params(axis='both', labelsize=15)
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.xlabel('x-direction', size=15)
    plt.ylabel('z-direction', size=15)
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(np.min(Y), np.max(Y))
    plt.tight_layout()
    plt.savefig(str(filename)+'-contour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

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
dset = f['z-direction']
zaxis = dset[:]

############################################################################################
#loop over all files#

for i in range (0, len(filename)):

    ########################################################################################
    #load#
    f = h5py.File('../outfile/'+filename[i], 'r')
    dset = f['primitive']
    primitive=dset[:]
    primitive = primitive.T
    neq = primitive.shape[0]

    # assign #
    rho = primitive[0,:,0,:].T
    velx = primitive[1,:,0,:].T
    vely = primitive[2,:,0,:].T
    velz = primitive[3,:,0,:].T
    p = primitive[4,:,0,:].T
    bx = primitive[neq-3,:,0,:].T
    by = primitive[neq-2,:,0,:].T
    bz = primitive[neq-1,:,0,:].T
  
    #load#
    dset = f['epsilon']
    epsilon = dset[:]
    epsilon = epsilon.T
    epsilon = epsilon[:,0,:].T

    #assign#
    dset = f['time']
    time = dset[:][0]

    #magnetic fields#
    threshold = 1e-10
    for i in range (0, bx.shape[0]):
        for j in range(0, bx.shape[1]):
            if(abs(bx[i,j]) < threshold):
                bx[i,j] = 0
            if(abs(by[i,j]) < threshold):
                by[i,j] = 0
            if(abs(bz[i, j]) < threshold):
                bz[i, j] = 0

    ########################################################################################
    #plot#
    z = rho
    if(np.min(z) == np.max(z)):
        linear_plot(z,'rho')
    else:
        log_plot(z,'rho',stream=False)


