"""
Common functions to all 3D methods.
"""

def readCubes(file_name):
  """
  Read from file the coordinates of the full cubes
  and return them in list of ternes.
  """
  ternes = []
  f = open(file_name,"r")
  for L in f:
    coord = L.split()
    if len(coord)==0: continue
    x,y,z = [int(c) for c in coord]
    ternes.append((x,y,z))
  f.close()
  return ternes


def printMoments(MOME):
  #print("Momenti su quadtree con algoritmo di Wu et al.")
  #print(type(black_cubes))
  for p in range(4):
    for q in range(4):
      for r in range(4):
        if p+q+r>3: continue
        print("Moment ",p,q,r, " = ", int(MOME[(p,q,r)]))

#Global variable defining the moments to be computed
orders = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]
orders = orders + [(1,1,0),(1,0,1),(0,1,1),(2,0,0),(0,2,0),(0,0,2)]
orders = orders + [(0,1,2),(0,2,1),(1,2,0),(1,0,2),(2,0,1),(2,1,0)]
orders = orders + [(3,0,0),(0,3,0),(0,0,3),(1,1,1)]

def main(arg, decomposF, preprocF, momentsF, title):
   """
   Apply the function decompF to build a decomposition of the 3D
   image, and apply the function momentsF to compute the moments
   based on the decomposition.
   Print the string titleF, which describes the method.
   """
   
   print(title)
   if len(arg)<=1:
        print("Need 3D image name")

   else:
        print("---Read cubes from file "+ arg[1])
        input_cubes = readCubes(arg[1])
        print("Number of black voxels:", len(input_cubes))
        elements = decomposF(input_cubes)
        print("Number of black elements:", elements.num_elem())
        if hasattr(elements,"max_triplet"):
           print("Max sides:", elements.max_triplet())
        
        times_to_repeat = 1
        try:
           t = int(arg[2])
           if t>1: times_to_repeat = t
        except:
           pass
        for t in range(times_to_repeat):
           if preprocF!=None: preprocF(elements)
           moments = momentsF(elements)
        printMoments(moments)
