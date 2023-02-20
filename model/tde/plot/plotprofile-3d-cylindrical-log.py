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
def linear_xzhalf(prim,filename,b1,b2,stream=False):
    fig, ax = plt.subplots(1, 1)
    zmax = np.max(prim)
    zmin = np.min(prim)
    cp = ax.contourf(X_half, Z_half, prim, 100, cmap='plasma', extend='both')
    cbar = fig.colorbar(cp,format=tick.FormatStrFormatter('%.2f'))
    cbar.ax.tick_params(labelsize=15)
    ax.tick_params(axis='both', labelsize=15)
    if(stream):
        plt.streamplot(X_half, Z_half, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)
    for art in ax.get_children():
        if not isinstance(art, matplotlib.patches.FancyArrowPatch):
            continue
        art.remove()
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.xlabel('x-direction', size=15)
    plt.ylabel('z-direction', size=15)
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(X_half), np.max(X_half))
    plt.ylim(np.min(Z_half), np.max(Z_half))
    plt.tight_layout()
    plt.savefig(str(filename)+'-halfxzcontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

# define function #
def linear_xzplot(prim,filename):
    fig, ax = plt.subplots(1, 1)
    zmax = np.max(prim)
    zmin = np.min(prim)
    cp = ax.contourf(X, Z, prim, 100, cmap='plasma', extend='both')
    cbar = fig.colorbar(cp,format=tick.FormatStrFormatter('%.2f'))
    cbar.ax.tick_params(labelsize=15)
    ax.tick_params(axis='both', labelsize=15)
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.xlabel('x-direction', size=15)
    plt.ylabel('z-direction', size=15)
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(np.min(Z), np.max(Z))
    plt.tight_layout()
    plt.savefig(str(filename)+'-xzcontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

# define function #
def linear_xyplot(prim,filename):
    fig, ax = plt.subplots(1, 1)
    zmax = np.max(prim)
    zmin = np.min(prim)
    circle1 = plt.Circle((0, 0), R.min(), color='black')
    ax.add_patch(circle1)
    cp = ax.contourf(X_plane, Y_plane, prim, 100, cmap='plasma',extend='both')
    cbar = fig.colorbar(cp,format=tick.FormatStrFormatter('%.2f'))
    cbar.ax.tick_params(labelsize=15)
    ax.tick_params(axis='both', labelsize=15)
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.xlabel('x-direction', size=15)
    plt.ylabel('y-direction', size=15)
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(X_plane), np.max(X_plane))
    plt.ylim(np.min(Y_plane), np.max(Y_plane))
    plt.savefig(str(filename)+'-xycontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

############################################################################################
# define function #
def log_xzhalf(prim,filename,b1,b2,stream=False):
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
    if(stream):
        plt.streamplot(X_half, Z_half, b1, b2, density=1.5, linewidth=None, color='black', broken_streamlines=True)
    for art in ax.get_children():
        if not isinstance(art, matplotlib.patches.FancyArrowPatch):
            continue
        art.remove()
    ax.tick_params(axis='both', labelsize=15)
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.xlabel('x-direction', size=15)
    plt.ylabel('z-direction', size=15)
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(X_half), np.max(X_half))
    plt.ylim(np.min(Z_half), np.max(Z_half))
    plt.tight_layout()
    plt.savefig(str(filename)+'-xzcontourhalf-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

# define function #
def log_xzplot(prim,filename):
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
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.xlabel('x-direction', size=15)
    plt.ylabel('z-direction', size=15)
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(X), np.max(X))
    plt.ylim(np.min(Z), np.max(Z))
    plt.tight_layout()
    plt.savefig(str(filename)+'-xzcontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

# define function #
def log_xyplot(prim,filename):
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
    circle1 = plt.Circle((0, 0), R.min(), color='black')
    ax.add_patch(circle1)
    ax.tick_params(axis='both', labelsize=15)
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.xlabel('x-direction', size=15)
    plt.ylabel('y-direction', size=15)
    plt.gca().set_aspect('equal')
    plt.xlim(np.min(X_plane), np.max(X_plane))
    plt.ylim(np.min(Y_plane), np.max(Y_plane))
    plt.savefig(str(filename)+'-xycontour-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

############################################################################################
# define function #
def axis_duplicate(prim):
    p_temp = np.ndarray(shape=(2*prim.shape[0],prim.shape[1],prim.shape[2]), dtype=float)
    p_temp[0:int(p_temp.shape[0]/2),:,:] = prim[::-1,:,:]
    p_temp[int(p_temp.shape[0]/2):p_temp.shape[0],:,:] = prim[:,:,:]
    p_temp = np.average(p_temp, axis=1)
    return p_temp

# define function #
def polar_duplicate(prim):
    p_temp = np.ndarray(shape=(prim.shape[0],prim.shape[1]+1), dtype=float)
    p_temp[:,0:prim.shape[1]] = prim[:,:]
    p_temp[:,prim.shape[1]] = prim[:,0]
    return p_temp

# define function #
def axis_average(prim):
    p_temp = np.average(prim, axis=1)
    return p_temp

# define function #
def axis_phi0(prim,i_in):
    p_temp = prim[:,i_in,:]
    return p_temp

############################################################################################
#plotting function #
def plot(z_in,i_in,fname):
    if(np.min(z_in) == np.max(z_in)):
        z = axis_phi0(z_in,i_in).T
        b1 = axis_phi0(bx,i_in).T
        b2 = axis_phi0(bz,i_in).T
        fname = str(fname) + '-' + str(i_in)
        linear_xzhalf(z,str(fname), b1, b2, stream=True)
        z = axis_duplicate(z_in).T
        linear_xzplot(z,str(fname))
        z = polar_duplicate(z_in[:,:,int(z_in.shape[2]/2)]).T
        linear_xyplot(z,str(fname))
    else:
        z = axis_average(z_in).T
        b1 = axis_average(bx).T
        b2 = axis_average(bz).T
        log_xzhalf(z,str(fname), b1, b2, stream=True)
        z = axis_duplicate(z_in).T
        log_xzplot(z,str(fname))
        z = polar_duplicate(z_in[:,:,int(z_in.shape[2]/2)]).T
        log_xyplot(z,str(fname))

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
#gridding#

#mesh grid, x-y plane#
ytemp = np.zeros(yaxis.shape[0] + 1)
ytemp[0:yaxis.shape[0]] = yaxis[:]
ytemp[yaxis.shape[0]] = 2*math.pi + yaxis[0]
yaxis = ytemp
ytemp = []
R, phi = np.meshgrid(xaxis, yaxis)
X_plane, Y_plane = R*np.cos(phi), R*np.sin(phi)
phi, R = np.meshgrid(yaxis, xaxis)

#mesh grid, half x-z plane#
X_half, Z_half = np.meshgrid(xaxis, zaxis)

#mesh grid, full x-z plane#
xtemp = np.zeros(2*xaxis.shape[0])
xtemp[0:int(xtemp.shape[0]/2)] = -xaxis[::-1]
xtemp[int(xtemp.shape[0]/2):xtemp.shape[0]] = xaxis[:]
xaxis = xtemp
xtemp = []
X, Z = np.meshgrid(xaxis, zaxis)

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
    epsilon = dset[:]
    epsilon = epsilon
    epsilon = epsilon[:,:,:].T

    ########################################################################################
    #magnetic fields#

    threshold = 1e-10
    for i in range (0, bx.shape[0]):
        for j in range(0, bx.shape[1]):
            for k in range(0, bx.shape[2]):
                if(abs(bx[i,j,k]) < threshold):
                    bx[i,j,k] = 0
                if(abs(by[i,j,k]) < threshold):
                    by[i,j,k] = 0
                if(abs(bz[i,j,k]) < threshold):
                    bz[i,j,k] = 0

    ########################################################################################
    #plot#

    #bx#
    plot(rho,0,'rho')