#import packages#
import os
import moviepy.video.io.ImageSequenceClip

#sorting function #
def keyfunc(x):
  test = x.split('-')
  test = test[2].split('.png')
  test = float(test[0])
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
               if img.startswith("rho-xy")]
image_files = sorted(image_files, key=keyfunc)
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('rho-xy-3d.mp4')