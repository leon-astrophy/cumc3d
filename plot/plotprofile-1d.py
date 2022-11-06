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

    rho = primitive[:,0,0,0]
    velx = primitive[:,0,0,1]
    p = primitive[:,0,0,4]
    dset = f['x-direction']
    xaxis = dset[:]
    dset = f['time']
    time = dset[:][0]

    #plot#
    plt.plot(xaxis, rho)
    plt.grid()
    plt.title('Time = '+str(time))
    plt.xlabel('x-direction')
    plt.ylabel('Density')
    plt.savefig('rhoprofile-'+str(time)+'.png')
    plt.clf()

    #plot#
    plt.plot(xaxis, velx)
    plt.grid()
    plt.title('Time = '+str(time))
    plt.xlabel('x-direction')
    plt.ylabel('Velocity')
    plt.savefig('velprofile-'+str(time)+'.png')
    plt.clf()

    #plot#
    plt.plot(xaxis, p)
    plt.grid()
    plt.title('Time = '+str(time))
    plt.xlabel('x-direction')
    plt.ylabel('Pressure')
    plt.savefig('preprofile-'+str(time)+'.png')
    plt.clf()