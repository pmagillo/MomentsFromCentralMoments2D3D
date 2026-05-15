# MomentsFromCentralMoments2D3D
Python code for computing geometric moments of 2D and 3D images from precomputed central moments. The code uses quadtree/octree decomposition or block decomposition. It includes implementation of corresponding state-of-the-art algorithms for comparison purpose. Test data included.

MORE INFORMATION:

Computation of discrete geometric moments on 2D and 3D binary images, black=1 is the foreground (object) and white=0 is the background.

Approach based on decomposing the foreground into basic shapes:
squares (cubes in 3D) of a quadtree (octree) or
rectangles (cuboids in 3D) of a block decomposition.
The moments of each basic shape are computed and the overall moments are the sum of them.

The software is written in Python with no additional packages with the exception of reading2D.py which uses pillow and numpy to read 2D image.
However 2D images may also be read from raw binary format, and the part of code using the libraries can be commented away. This part is marked in the source file reading2D.py

It implements state-of-the-art methods based on quadtree/octree and on 2D/3D block decomposition, and
a new version which precomputes central moments and uses them for computing the ordinary moments.

Details on accepted data formats are in: DATA_FORMAT.TXT
Details of the implemented programs are in: PROGRAMS.TXT
Details for performing the experiments are in: EXPERIMENTS.TXT
Test data are in: DATA_2D, DATA_3D
