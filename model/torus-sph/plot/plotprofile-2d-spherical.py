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
#function for plot#
def plotxz_linear(prim,fname,b1,b2,stream=False):
    fig, ax = plt.subplots(1, 1)
    cp = ax.contourf(X, Z, prim, 100, cmap='plasma')
    cbar = fig.colorbar(cp,format=tick.FormatStrFormatter('%.2f'))
    cbar.ax.tick_params(labelsize=15)
    ax.tick_params(axis='both', labelsize=15)
    if(stream):
        plt.streamplot(X, Z, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)
    for art in ax.get_children():
        if not isinstance(art, matplotlib.patches.FancyArrowPatch):
            continue
        art.remove()
    circle1 = plt.Circle((0, 0), R.min(), color='black')
    ax.add_patch(circle1)
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.xlabel('x-direction', size=15)
    plt.ylabel('y-direction', size=15)
    plt.gca().set_aspect('equal')
    plt.xlim(0, np.max(X))
    plt.ylim(np.min(Z), np.max(Z))
    plt.tight_layout()
    plt.savefig(str(fname)+'-xzcontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

#function for plot#
def plotxz_log(prim,fname,b1,b2,stream=False):
    fig, ax = plt.subplots(1, 1)
    zmax = int(np.log10(np.max(prim)))
    zmin = int(np.log10(np.min(prim)))
    cp = ax.contourf(X, Z, prim, np.logspace(zmin, zmax, 100), locator=ticker.LogLocator(), cmap='plasma',extend='both')
    cbar = fig.colorbar(cp)
    rang = np.arange(int(np.log10(np.min(prim))), int(np.log10(np.max(prim))) + 1, 1)
    loca = 10 ** (np.array(rang).astype(float))
    cbar.set_ticks(loca)
    cbar.minorticks_off()
    labels = ['10$^{%.0f}$' % x for x in rang]
    cbar.ax.set_yticklabels(labels, fontsize=15)
    ax.tick_params(axis='both', labelsize=15)
    if(stream):
        plt.streamplot(X, Z, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)
    for art in ax.get_children():
        if not isinstance(art, matplotlib.patches.FancyArrowPatch):
            continue
        art.remove()
    circle1 = plt.Circle((0, 0), R.min(), color='black')
    ax.add_patch(circle1)
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.xlabel('x-direction', size=15)
    plt.ylabel('z-direction', size=15)
    plt.gca().set_aspect('equal')
    plt.xlim(0, np.max(X))
    plt.ylim(np.min(Z), np.max(Z))
    plt.tight_layout()
    plt.savefig(str(fname)+'-xzcontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

############################################################################################
# function for patching poles #
def patch_pole(prim):
    p_temp = np.ndarray(shape=(prim.shape[0],prim.shape[1]+2), dtype=float)
    p_temp[:,1:prim.shape[1] + 1] = prim[:,:]
    p_temp[:,p_temp.shape[1]-1] = p_temp[:,p_temp.shape[1]-2]
    p_temp[:,0] = p_temp[:,1]
    return p_temp

############################################################################################
#plotting function #
def plot(z_in,fname,stream=False):
    z = patch_pole(z_in)
    b1 = patch_pole(bx)
    b2 = patch_pole(bz)
    if(np.min(z_in) == np.max(z_in)):
        plotxz_linear(z,fname,b1,b2,stream=stream)
    else:
        plotxz_log(z,fname,b1,b2,stream=stream)

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

#patch data at poles#
ytemp = np.zeros(yaxis.shape[0] + 2)
ytemp[1:yaxis.shape[0] + 1] = yaxis[:]
ytemp[ytemp.shape[0]-1] = ytemp[1] + math.pi
ytemp[0] = - ytemp[1]
yaxis = ytemp
ytemp = []

#mesh grid#
R, Theta = np.meshgrid(xaxis, yaxis)
X = R * np.sin(Theta)
Z = R * np.cos(Theta)
X = X.T
Z = Z.T

############################################################################################
#loop#
for i in range (0, len(filename)):

    #load#
    f = h5py.File('../outfile/'+filename[i], 'r')
    dset = f['primitive']
    primitive=dset[:]
    primitive = primitive.T
    neq = primitive.shape[0]

    #assign#
    dset = f['time']
    time = dset[:][0]

    # assign #
    rho = primitive[0,:,:,0]
    velx = primitive[1,:,:,0]
    vely = primitive[2, :, :, 0]
    velz = primitive[3, :, :, 0]
    p = primitive[4,:,:,0]
    bx = primitive[neq-3,:,:,0]
    by = primitive[neq-2,:,:,0]
    bz = primitive[neq-1, :, :, 0]

    #load#
    dset = f['epsilon']
    epsilon=dset[:]
    epsilon = epsilon.T
    epsilon = epsilon[:, :, 0]

    ########################################################################################
    #plot#

    #density#
    plot(rho,'rho',stream=False)