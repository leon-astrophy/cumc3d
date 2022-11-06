import os
import h5py
import numpy as np
import matplotlib.pyplot as plt

# filename#
filename = []

for root, dirs, files in os.walk('../outfile/'):
    for file in files:
        if file.endswith(".hdf5"):
            filename.append(os.path.join(file))

#loop#
for i in range (0, len(filename)):
    f = h5py.File('../outfile/'+filename[i], 'r')
    dset = f['primitive']
    primitive=dset[:]
    primitive = primitive.T

    rho = primitive[:,:,0,0]
    velx = primitive[:,:,0,1]
    p = primitive[:,:,0,4]
    dset = f['x-direction']
    xaxis = dset[:]
    dset = f['y-direction']
    yaxis = dset[:]
    dset = f['time']
    time = dset[:][0]

    X, Y = np.meshgrid(xaxis, yaxis)
    fig, ax = plt.subplots(1, 1)
    cp = ax.contourf(X, Y, rho, 1000)
    fig.colorbar(cp)
    plt.title('Time = ' + str(time))
    plt.xlabel('x-direction')
    plt.ylabel('y-direction')
    plt.savefig('rhocontour-'+str(time)+'.png')
    plt.clf()
    plt.close()