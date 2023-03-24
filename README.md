# Image-Processing
3 MATLAB scripts that filter noise from an image, piece back and align a cut image and carry out edge detection on an image, respectively.

The image filtration folder contains the matlab script as well as the input photo and the output filtrations. Filtration was obtained by using fourier transforms to remove sinusoidal noise. A PDF document also shows the step by step process as well as the theory behind the code.

The image stitching folder contains a folder with 4 input image-pieces of the same map, taken from different perspectives as well as the matlab script and helper function files along with some output images. The task was to stitch together the 4 images relative to a main piece and have a fully alligned photo of the map using the Computer Vision toolbox in Matlab. The final output result is shown as image "croppedmap.jpg". A PDF document also shows the step by step process as well as the theory behind the code.

The edge detection folder contains a matlab script and an input photo as well as some function helper files and output images. The task was to take the input image of blood vessels and perform edge detection using Marr Hildreth edge detection and canny edge detection after applying a gaussian blur and getting a binary image. The final result can be seen in the images "ex3_7_canny.jpg" and "ex3_7_marr.jpg". A PDF document also shows the step by step process as well as the theory behind the code.
