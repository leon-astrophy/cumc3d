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
def plotxz_half_linear(prim,fname,b1,b2,stream=False):
    fig, ax = plt.subplots(1, 1)
    cp = ax.contourf(X_half, Z_half, prim, 100, cmap='plasma')
    cbar = fig.colorbar(cp,format=tick.FormatStrFormatter('%.2f'))
    cbar.ax.tick_params(labelsize=15)
    ax.tick_params(axis='both', labelsize=15)
    if(stream):
        plt.streamplot(X_half, Z_half, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)
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
    plt.xlim(0, np.max(X_half))
    plt.ylim(np.min(Z_half), np.max(Z_half))
    plt.tight_layout()
    plt.savefig(str(fname)+'-xzhalfcontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

#function for plot#
def plotxz_half_log(prim,fname,b1,b2,stream=False):
    fig, ax = plt.subplots(1, 1)
    zmax = int(np.log10(np.max(prim)))
    zmin = int(np.log10(np.min(prim)))
    cp = ax.contourf(X_half, Z_half, prim, np.logspace(zmin, zmax, 100), locator=ticker.LogLocator(), cmap='plasma',extend='both')
    cbar = fig.colorbar(cp)
    rang = np.arange(int(np.log10(np.min(prim))), int(np.log10(np.max(prim))) + 1, 1)
    loca = 10 ** (np.array(rang).astype(float))
    cbar.set_ticks(loca)
    cbar.minorticks_off()
    labels = ['10$^{%.0f}$' % x for x in rang]
    cbar.ax.set_yticklabels(labels, fontsize=15)
    ax.tick_params(axis='both', labelsize=15)
    if(stream):
        plt.streamplot(X_half, Z_half, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)
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
    plt.xlim(0, np.max(X_half))
    plt.ylim(np.min(Z_half), np.max(Z_half))
    plt.tight_layout()
    plt.savefig(str(fname)+'-xzhalfcontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

############################################################################################
#function for plot#
def plotxz_full_linear(prim,fname,b1,b2,stream=False):
    fig, ax = plt.subplots(1, 1)
    cp = ax.contourf(X_full, Z_full, prim, 100, cmap='plasma')
    cbar = fig.colorbar(cp,format=tick.FormatStrFormatter('%.2f'))
    cbar.ax.tick_params(labelsize=15)
    ax.tick_params(axis='both', labelsize=15)
    if(stream):
        plt.streamplot(X_full, Z_full, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)
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
    plt.xlim(np.min(X_full), np.max(X_full))
    plt.ylim(np.min(Z_full), np.max(Z_full))
    plt.tight_layout()
    plt.savefig(str(fname)+'-xzfullcontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

#function for plot#
def plotxz_full_log(prim,fname,b1,b2,stream=False):
    fig, ax = plt.subplots(1, 1)
    zmax = int(np.log10(np.max(prim)))
    zmin = int(np.log10(np.min(prim)))
    cp = ax.contourf(X_full, Z_full, prim, np.logspace(zmin, zmax, 100), locator=ticker.LogLocator(), cmap='plasma',extend='both')
    cbar = fig.colorbar(cp)
    rang = np.arange(int(np.log10(np.min(prim))), int(np.log10(np.max(prim))) + 1, 1)
    loca = 10 ** (np.array(rang).astype(float))
    cbar.set_ticks(loca)
    cbar.minorticks_off()
    labels = ['10$^{%.0f}$' % x for x in rang]
    cbar.ax.set_yticklabels(labels, fontsize=15)
    ax.tick_params(axis='both', labelsize=15)
    if(stream):
        plt.streamplot(X_full, Z_full, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)
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
    plt.xlim(np.min(X_full), np.max(X_full))
    plt.ylim(np.min(Z_full), np.max(Z_full))
    plt.tight_layout()
    plt.savefig(str(fname)+'-xzfullcontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

############################################################################################
#function for plot#
def plotxy_plane_linear(prim,fname,b1,b2,stream=False):
    fig, ax = plt.subplots(1, 1)
    cp = ax.contourf(X_plane, Y_plane, prim, 100, cmap='plasma')
    cbar = fig.colorbar(cp,format=tick.FormatStrFormatter('%.2f'))
    cbar.ax.tick_params(labelsize=15)
    ax.tick_params(axis='both', labelsize=15)
    if(stream):
        plt.streamplot(X_plane, Y_plane, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)
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
    plt.xlim(np.min(X_plane), np.max(X_plane))
    plt.ylim(np.min(Y_plane), np.max(Y_plane))
    plt.tight_layout()
    plt.savefig(str(fname)+'-xyplanecontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

#function for plot#
def plotxy_plane_log(prim,fname,b1,b2,stream=False):
    fig, ax = plt.subplots(1, 1)
    zmax = int(np.log10(np.max(prim)))
    zmin = int(np.log10(np.min(prim)))
    cp = ax.contourf(X_plane, Y_plane, prim, np.logspace(zmin, zmax, 100), locator=ticker.LogLocator(), cmap='plasma',extend='both')
    cbar = fig.colorbar(cp)
    rang = np.arange(int(np.log10(np.min(prim))), int(np.log10(np.max(prim))) + 1, 1)
    loca = 10 ** (np.array(rang).astype(float))
    cbar.set_ticks(loca)
    cbar.minorticks_off()
    labels = ['10$^{%.0f}$' % x for x in rang]
    cbar.ax.set_yticklabels(labels, fontsize=15)
    ax.tick_params(axis='both', labelsize=15)
    if(stream):
        plt.streamplot(X_plane, Y_plane, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)
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
    plt.xlim(np.min(X_plane), np.max(X_plane))
    plt.ylim(np.min(Y_plane), np.max(Y_plane))
    plt.tight_layout()
    plt.savefig(str(fname)+'-xyplanecontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

############################################################################################
# function for patching poles #
def patch_pole_half(prim):
    p_temp = np.ndarray(shape=(prim.shape[0],prim.shape[1]+2), dtype=float)
    p_temp[:,1:prim.shape[1] + 1] = prim[:,:,0]
    p_temp[:,p_temp.shape[1]-1] = p_temp[:,p_temp.shape[1]-2]
    p_temp[:,0] = p_temp[:,1]
    return p_temp

# function for patching poles #
def patch_full(prim):
    p_temp = np.ndarray(shape=(prim.shape[0],prim.shape[1]*2+1), dtype=float)
    p_temp[:,0:prim.shape[1]] = prim[:,:,prim.shape[2]-1]
    p_temp[:,prim.shape[1]:prim.shape[1]*2] = prim[:,:,0]
    p_temp[:,prim.shape[1]*2] = prim[:,prim.shape[1]-1,0]
    return p_temp

# function for patching poles #
def patch_plane(prim):
    p_temp = np.ndarray(shape=(prim.shape[0],prim.shape[2]+1), dtype=float)
    p_temp[:,0:p_temp.shape[1]-1] = prim[:,int(prim.shape[1]/2),:]
    p_temp[:,p_temp.shape[1]-1] = prim[:,int(prim.shape[1]/2),0]
    return p_temp

############################################################################################
#plotting function #
def plot(z_in,fname,stream=False):
    if(np.min(z_in) == np.max(z_in)):
        z = patch_pole_half(z_in)
        b1 = patch_pole_half(bx)
        b2 = patch_pole_half(by)
        plotxz_half_linear(z,fname,b1,b2,stream=stream)
        z = patch_full(z_in)
        b1 = patch_full(bx)
        b2 = patch_full(by)
        plotxz_full_linear(z,fname,b1,b2,stream=stream)
        z = patch_plane(z_in)
        b1 = patch_plane(bx)
        b2 = patch_plane(bz)
        plotxy_plane_linear(z,fname,b1,b2,stream=stream)
    else:
        z = patch_pole_half(z_in)
        b1 = patch_pole_half(bx)
        b2 = patch_pole_half(by)
        plotxz_half_log(z,fname,b1,b2,stream=stream)
        z = patch_full(z_in)
        b1 = patch_full(bx)
        b2 = patch_full(by)
        plotxz_full_log(z,fname,b1,b2,stream=stream)
        z = patch_plane(z_in)
        b1 = patch_plane(bx)
        b2 = patch_plane(bz)
        plotxy_plane_log(z,fname,b1,b2,stream=stream)

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

#create full angle#
ytemp = np.zeros(2*yaxis.shape[0]+1)
ytemp[yaxis.shape[0]:2*yaxis.shape[0]] = yaxis[:]
ytemp[0:yaxis.shape[0]] = -yaxis[::-1]
ytemp[2*yaxis.shape[0]] = ytemp[0] + 2*math.pi

#mesh grid#
R, Theta = np.meshgrid(xaxis, ytemp)
X_full = R * np.sin(Theta)
Z_full = R * np.cos(Theta)
X_full = X_full.T
Z_full = Z_full.T

#patch data at poles#
ytemp = np.zeros(yaxis.shape[0] + 2)
ytemp[1:yaxis.shape[0] + 1] = yaxis[:]
ytemp[ytemp.shape[0]-1] = ytemp[1] + math.pi
ytemp[0] = - ytemp[1]

#mesh grid#
R, Theta = np.meshgrid(xaxis, ytemp)
X_half = R * np.sin(Theta)
Z_half = R * np.cos(Theta)
X_half = X_half.T
Z_half = Z_half.T

#patch data at azimutal#
ztemp = np.zeros(zaxis.shape[0] + 1)
ztemp[0:zaxis.shape[0]] = zaxis[:]
ztemp[ztemp.shape[0]-1] = ztemp[0] + 2*math.pi

#mesh grid#
R, Theta = np.meshgrid(xaxis, ztemp)
X_plane = R * np.cos(Theta)
Y_plane = R * np.sin(Theta)
X_plane = X_plane.T
Y_plane = Y_plane.T

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
    rho = primitive[0,:,:,:]
    velx = primitive[1,:,:,:]
    vely = primitive[2,:,:,:]
    velz = primitive[3,:,:,:]
    p = primitive[4,:,:,:]
    bx = primitive[neq-3,:,:,:]
    by = primitive[neq-2,:,:,:]
    bz = primitive[neq-1,:,:,:]

    #load#
    dset = f['epsilon']
    epsilon=dset[:]
    epsilon = epsilon.T
    epsilon = epsilon[:,:,:]

    ########################################################################################
    #plot#

    #density#
    plot(rho,'rho',stream=False)