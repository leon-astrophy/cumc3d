import os
import moviepy.video.io.ImageSequenceClip
image_folder='./'
fps=10

image_files = [os.path.join(image_folder,img)
               for img in os.listdir(image_folder)
               if img.startswith("rho")]
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps=fps)
clip.write_videofile('rhomovie.mp4')