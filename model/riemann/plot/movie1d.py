#import packages#
import os
import moviepy.video.io.ImageSequenceClip

#sorting function #
def keyfunc(x):
  test = x.split('-')
  test = test[1].split('.png')
  test = test[0]
  return test

#folder#
image_folder='./'

#frame per second#
fps=10

#####################################################################################
# magnetic fields #
image_files = []
image_files = [os.path.join(image_folder,img)
               for img in os.listdir(image_folder)
               if img.startswith("bfieldprofile-")]
image_files = sorted(image_files, key=keyfunc)
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('bfield-1d.mp4')

#####################################################################################
# pressure density fields #
image_files = []
image_files = [os.path.join(image_folder,img)
               for img in os.listdir(image_folder)
               if img.startswith("prhoprofile-")]
image_files = sorted(image_files, key=keyfunc)
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('prho-1d.mp4')

#####################################################################################
# velocity fields #
image_files = []
image_files = [os.path.join(image_folder,img)
               for img in os.listdir(image_folder)
               if img.startswith("velprofile-")]
image_files = sorted(image_files, key=keyfunc)
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('vel-1d.mp4')