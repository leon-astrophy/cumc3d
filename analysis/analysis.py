############################################
# Plotting utilties for CUMC3D simulations
############################################

#import#
import os
import sys
import h5py
import numpy as np
from subprocess import call, Popen

################################################################################

print('\t')
print('--------------------------------------------')
print('- Welcome to plotting utilities for CUMC3D -')
print('--------------------------------------------')
print('\t')

################################################################################

# definine function #
def getfname(x):
  x = x.split('/')
  x = x[2]
  return x

################################################################################

#define path#
path = '../model/'

#get all model directories#
models = [f.path for f in os.scandir(path) if f.is_dir()]

#get model names#
models = [getfname(x) for x in models]

#print out#
print('Here are the models available:'+'\n')
for i in range (0, len(models)):
  print(str(i+1)+'.'+models[i])
print('\t')
print('Choose the model for plotting analysis:')
print('\t')

# take input #
target = int(input())

#check error input#
if(target < 1 or target > len(models)):
  print('\t')
  print('Error, wrong indicies')
  exit()

################################################################################

#set target models#
models = models[target-1]
print('\t')
print('The selected model is '+str(models))

#set path#
path = str(path) + str(models) + '/'

################################################################################

#arrays for coordinate and dimensions#
n_dim = ['1d', '2d', '3d']
coords = ['cartesian', 'cylindrical', 'spherical']

# read grid parameter files, get model dimension and coordinate #
f = h5py.File(path+'outfile/grid_param.hdf5', 'r')
dset = f['dimension']
dim = int(dset[:][0])
dset = f['coordinate']
coord = int(dset[:][0])

################################################################################

# filename#
filename = []

#load#
for root, dirs, files in os.walk(path+'outfile/'):
    for file in files:
        if file.endswith("nm.hdf5"):
            filename.append(os.path.join(file))

################################################################################

# assign #
n_dim = n_dim[dim - 1]
coords = coords[coord - 1]

# run process according to user input #
print('\t')
print('Now, run the plotting script ...')
print('\t')
if(n_dim == '1d'):
    script = './plot/plotprofile-1d.py'
else:
    script = './plot/plotcontour-'+str(coords)+'.py'

#run python script, loop over filename #
for i in range (0, len(filename)):
  hdf5_file = path+'outfile/'+filename[i]
  grid_file = path+'outfile/grid_param.hdf5'
  Popen(['python3', script, grid_file, hdf5_file], close_fds=True)

################################################################################

print('-------------------------')
print('- End plotting analysis -')
print('-------------------------')
print('\t')

################################################################################
