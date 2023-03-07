############################################
# Generate animations from figures
############################################

#import#
import os
import sys
import h5py
import numpy as np
from subprocess import call
import moviepy.video.io.ImageSequenceClip

################################################################################

print('\t')
print('---------------------------------------------')
print('- Welcome to animation utilities for CUMC3D -')
print('---------------------------------------------')
print('\t')

################################################################################

# definine function #
def getfname(x):
  x = x.split('-')
  x = x[0] + '-' + x[1]
  return x

#sorting function #
def keyfunc(x):
  test = x.split('-')
  test = test[2].split('.png')
  test = float(test[0])
  return test

################################################################################

#define path#
path = './figure/'

#load unique filename in figure folder#
filename = []
for root, dirs, files in os.walk(path):
    for file in files:
        if file.endswith(".png"):
            filename.append(os.path.join(file))

#get model names#
filename = [getfname(x) for x in filename]
filename = np.unique(filename)

################################################################################
# generate animation

#frame per second#
fps=10

#now loop over all unique png file#
for i in range (0, len(filename)):
  image_files = []
  image_files = [os.path.join(path,img)
                for img in os.listdir(path)
                if (img.startswith(filename[i]) and img.endswith(".png"))]
  image_files = sorted(image_files, key=keyfunc)
  clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
  clip.write_videofile('./movie/'+filename[i]+'.mp4')

################################################################################

print('---------------------------')
print('- End animation generator -')
print('---------------------------')
print('\t')

################################################################################