# MomentsFromCentralMoments2D3D
Python code for computing geometric moments of 2D and 3D images from precomputed central moments. The code uses quadtree/octree decomposition or block decomposition. It includes implementation of corresponding state-of-the-art algorithms for comparison purpose. Test data included.

MORE INFORMATION:

Computation of discrete geometric moments on 2D and 3D binary images, 
black is the foreground (object) and white is the background.

Approach based on decomposing the foreground into basic shapes:
squares (cubes in 3D) of a quadtree (octree) or
rectangles (cuboids in 3D) of a block decomposition.
The moments of each basic shape are computed and the overall
moments are the sum of them.

The software is written in Python with no additional packages.

It implements state-of-the-art methods based on quadtree/octree
and on 2D/3D block decomposition, and
a new version which precomputes central moments and uses them
for computing the ordinary moments.

Details of the implemented programs are in:
PROGRAMS.TXT

Details of the experiments (including used test data) are in:
EXPERIMENTS.TXT
