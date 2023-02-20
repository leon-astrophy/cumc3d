#import packages#
import os
import moviepy.video.io.ImageSequenceClip

#sorting function #
def keyfunc(x):
  test = x.split('-')
  test = test[3].split('.png')
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
               if img.startswith("rho-bstream-")]
image_files = sorted(image_files, key=keyfunc)
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('rhobstream-2d.mp4')

#####################################################################################
# magnetic fields #
image_files = []
image_files = [os.path.join(image_folder,img)
               for img in os.listdir(image_folder)
               if img.startswith("rho-vstream-")]
image_files = sorted(image_files, key=keyfunc)
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('rhovstream-2d.mp4')

#####################################################################################
# magnetic fields #
image_files = []
image_files = [os.path.join(image_folder,img)
               for img in os.listdir(image_folder)
               if img.startswith("p-bstream-")]
image_files = sorted(image_files, key=keyfunc)
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('pbstream-2d.mp4')

#####################################################################################
# magnetic fields #
image_files = []
image_files = [os.path.join(image_folder,img)
               for img in os.listdir(image_folder)
               if img.startswith("p-vstream-")]
image_files = sorted(image_files, key=keyfunc)
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('pvstream-2d.mp4')
