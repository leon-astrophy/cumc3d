############################################################################################
# Plotting utilities for spherical 3D output from CUMC3D
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

# get phi component vector potential #
def get_aphi_half(vr_in, vth_in):
    a_phi = np.ndarray(shape=(xaxis.shape[0],yaxis.shape[0]), dtype=float)
    for i in range (0, a_phi.shape[0]):
        for j in range (0, a_phi.shape[1]):

            #integrate the Bz component #
            if(i == 0):
                int1 = 0
            else:
                dummy_var = xaxis[0:i]
                integrand = xaxis[0:i]*np.reshape(vth_in[0:i,0], (-1,))
                int1 = np.trapz(integrand, x=dummy_var)
                int1 = -int1*np.sin(yaxis[0])/xaxis[i]/np.sin(yaxis[j])

            # integrate the Br component #
            if(j == 0):
                int2 = 0
            else:
                dummy_var = yaxis[0:j]
                integrand = np.sin(yaxis[0:j])*np.reshape(vr_in[i,0:j], (-1,))
                int2 = np.trapz(integrand, x=dummy_var)
                int2 = int2*xaxis[i]/np.sin(yaxis[j])

            # sum them up #
            a_phi[i,j] = int1 + int2

    a_phi = a_phi - np.min(a_phi)
    return a_phi

# get phi component vector potential #
def get_aphi_full(aphi1, aphi2):
    a_phi = np.ndarray(shape=(xaxis.shape[0],2*yaxis.shape[0]+1), dtype=float)
    a_phi[:,0:yaxis.shape[0]] = aphi1[:,::-1]
    a_phi[:,yaxis.shape[0]:2*yaxis.shape[0]] = aphi2[:,:]
    a_phi[:,2*yaxis.shape[0]] = a_phi[:,0]
    a_phi = a_phi - np.min(a_phi)
    return a_phi

############################################################################################
# define function #

# get stream function #
def stream_function(v1_in, v2_in, domain='xzhalf'):
    if(domain=='xzhalf'):
        v1 = angular_average(v1_in)
        v2 = angular_average(v2_in)
        avec = get_aphi_half(v1,v2)
        avec = np.reshape(avec, (avec.shape[0],avec.shape[1], 1))
        avec = patch_pole_xzhalf(avec)

    elif(domain=='xzfull'):
        v1 = v1_in[:,:,int(v1_in.shape[2]/2)]
        v2 = v2_in[:,:,int(v2_in.shape[2]/2)]
        aphi1 = get_aphi_half(v1,v2)
        v1 = v1_in[:,:,0]
        v2 = v2_in[:,:,0]
        aphi2 = get_aphi_half(v1,v2)
        avec = get_aphi_full(aphi1, aphi2).T
        if(avec.shape != X_full.shape):
            avec = avec.T

    return avec

############################################################################################
# define function #

# plot contour #
def plot_contour(X1_in,X2_in,prim,filename,domain='xzhalf',scale='log',stream=False,bh=False):

    # set figure #
    fig, ax = plt.subplots(1, 1)

    # set label #
    if(domain == 'xyplane'):
        plt.xlabel('x-direction', size=15)
        plt.ylabel('y-direction', size=15)
    else:
        plt.xlabel('x-direction', size=15)
        plt.ylabel('z-direction', size=15)

    # log plot #
    if(scale == 'log'):
        zmax = int(np.log10(np.max(prim)))
        zmin = int(np.log10(np.min(prim)))
        cp = ax.contourf(X1_in, X2_in, prim, np.logspace(zmin, zmax, 100), locator=ticker.LogLocator(), cmap='RdGy',extend='both')
        cbar = fig.colorbar(cp)
        rang = np.arange(int(np.log10(np.min(prim))), int(np.log10(np.max(prim))) + 1, 1)
        loca = 10 ** (np.array(rang).astype(float))
        cbar.set_ticks(loca)
        cbar.minorticks_off()
        labels = ['10$^{%.0f}$' % x for x in rang]
        cbar.ax.set_yticklabels(labels, fontsize=15)

    # linear plot #
    elif(scale == 'linear'):
        cp = ax.contourf(X1_in, X2_in, prim, 100, cmap='RdGy', extend='both')
        cbar = fig.colorbar(cp)
        cbar.ax.tick_params(labelsize=15)

    # want stream lines? #
    if(stream):
        avec = stream_function(bx,by,domain=domain)
        plt.contour(X1_in,X2_in,avec,levels=10,colors='black')

    # draw black hole #
    if(bh):
        circle1 = plt.Circle((0, 0), np.min(R), color='black')
        ax.add_patch(circle1)

    # set domain #
    if(domain == 'xzhalf'):
        plt.xlim(0, np.max(X1_in))
        plt.ylim(np.min(X2_in), np.max(X2_in))
    else:
        plt.xlim(np.min(X1_in), np.max(X1_in))
        plt.ylim(np.min(X2_in), np.max(X2_in))       

    # further plot parameters # 
    ax.tick_params(axis='both', labelsize=15)
    plt.title('Time = ' + '%.3f' % time, size=15)
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(imgdir+str(filename)+'-'+str(domain)+'-'+'%.3f' % time +'.png')
    plt.clf()
    plt.close()

############################################################################################
#define function#

# array for angular average #
def angular_average(prim):
    p_temp = np.average(prim, axis=2)
    return p_temp

# function for patching poles #
def patch_pole_xzhalf(prim):
    p_temp = np.ndarray(shape=(prim.shape[0],prim.shape[1]+2), dtype=float)
    prim = angular_average(prim)
    p_temp[:,1:prim.shape[1] + 1] = prim[:,:]
    p_temp[:,p_temp.shape[1]-1] = p_temp[:,p_temp.shape[1]-2]
    p_temp[:,0] = p_temp[:,1]
    p_temp = p_temp.T
    if(p_temp.shape != X_half.shape):
        p_temp = p_temp.T
    return p_temp

# function for patching poles #
def patch_pole_xzfull(prim):
    p_temp = np.ndarray(shape=(prim.shape[0],prim.shape[1]*2 + 1), dtype=float)
    p_temp[:,0:prim.shape[1]] = prim[:,:,int(prim.shape[2]/2)]
    p_temp[:,prim.shape[1]:prim.shape[1]*2] = prim[:,:,0]
    p_temp[:,2*prim.shape[1]] = prim[:,0,0]
    p_temp = p_temp.T
    if(p_temp.shape != X_full.shape):
        p_temp = p_temp.T
    return p_temp

# function for patching poles #
def patch_pole_xyplane(prim):
    p_temp = np.ndarray(shape=(prim.shape[0],prim.shape[2]+1), dtype=float)
    p_temp[:,0:p_temp.shape[1]-1] = prim[:,int(prim.shape[1]/2),:]
    p_temp[:,p_temp.shape[1]-1] = prim[:,int(prim.shape[1]/2),0]
    p_temp = p_temp.T
    if(p_temp.shape != X_plane.shape):
        p_temp = p_temp.T
    return p_temp

############################################################################################
#plotting function #
def plot(z_in,fname,stream=False):

    # choose plotting scale #
    if((np.min(z_in) == np.max(z_in)) or (np.min(z_in) <= 0) or (np.max(z_in) <= 0) or np.isnan(z_in).any()):
        scale = 'linear'
    else:
        scale = 'log'
    
    ##########################################################################################
    # get angular averaged quantity
    z = patch_pole_xzhalf(z_in)

    # plot #
    plot_contour(X_half,Z_half,z,fname,domain='xzhalf',scale=scale,stream=stream,bh=True)

    ##########################################################################################
    # get x-z plane full quantity
    z = patch_pole_xzfull(z_in)

    #plot#
    plot_contour(X_full,Z_full,z,fname,domain='xzfull',scale=scale,stream=stream,bh=True)

    ##########################################################################################
    # get x-y plane full quantity
    z = patch_pole_xyplane(z_in)

    #plot#
    plot_contour(X_plane,Y_plane,z,fname,domain='xyplane',scale=scale,stream=False,bh=True)

############################################################################################
#load command line parameters #

path=sys.argv[1]

############################################################################################
# file input output #

# get path #
imgdir = './figure/'

# filename#
filename = []

#load file names#
for root, dirs, files in os.walk(path+'outfile/'):
    for file in files:
        if file.endswith("nm.hdf5"):
            filename.append(os.path.join(file))

#grid setting#
f = h5py.File(path+'outfile/grid_param.hdf5', 'r')
dset = f['x-direction']
xaxis = dset[:]
dset = f['y-direction']
yaxis = dset[:]
dset = f['z-direction']
zaxis = dset[:]

############################################################################################

#create full angle#
ytemp = np.zeros(2*yaxis.shape[0] + 1)
ytemp[2*yaxis.shape[0]] = yaxis[0] + math.pi
ytemp[yaxis.shape[0]:2*yaxis.shape[0]] = yaxis[:]
ytemp[0:yaxis.shape[0]] =  2*math.pi - yaxis[::-1]

#mesh grid#
R, Theta = np.meshgrid(xaxis, ytemp)
X_full = R * np.sin(Theta)
Z_full = R * np.cos(Theta)
X_full = X_full
Z_full = Z_full

#patch data at poles#
ytemp = np.zeros(yaxis.shape[0] + 2)
ytemp[1:yaxis.shape[0] + 1] = yaxis[:]
ytemp[ytemp.shape[0]-1] = ytemp[1] + math.pi
ytemp[0] = - ytemp[1]

#mesh grid#
R, Theta = np.meshgrid(xaxis, ytemp)
X_half = R * np.sin(Theta)
Z_half = R * np.cos(Theta)
X_half = X_half
Z_half = Z_half

#patch data at azimutal#
ztemp = np.zeros(zaxis.shape[0] + 1)
ztemp[0:zaxis.shape[0]] = zaxis[:]
ztemp[ztemp.shape[0]-1] = ztemp[0] + 2*math.pi

#mesh grid#
R, Theta = np.meshgrid(xaxis, ztemp)
X_plane = R * np.cos(Theta)
Y_plane = R * np.sin(Theta)
X_plane = X_plane
Y_plane = Y_plane

############################################################################################
# main trunck #

#loop over all files #
for i in range (0, len(filename)):

    #load#
    f = h5py.File(path+'outfile/'+filename[i], 'r')
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

    #load epsilon #
    dset = f['epsilon']
    epsilon=dset[:]
    epsilon = epsilon.T
    epsilon = epsilon[:,:,:]

    ########################################################################################
    #plot#

    # density #
    z = rho
    plot(z,'rho',stream=True)

    # toroidal magnetic field #
    z = bz
    plot(z,'bphi',stream=False)

    # plasma beta #
    z = 2*p/(bx**2 + by**2 + bz**2 + 1e-10)
    plot(z,'beta',stream=False)