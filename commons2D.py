"""
Common functions to all 2D methods.
"""

def printMoments(MOME):
  for p in range(4):
    for q in range(4):
        if p+q>3: continue
        print("Moment ",p,q, " = ", int(MOME[(p,q)]))
        #print(type(MOME[(p,q)]))#***************
        #OK #assert type(MOME[(p,q)]) is int
        
# Global variable defining the moments to be computed
orders = [(0,0),(1,0),(0,1)]
orders = orders + [(1,1),(2,0),(0,2)]
orders = orders + [(3,0),(0,3),(2,1),(1,2)]

from reading2D import adaptiveRead as readInput

def main(arg, decomposF, preprocF, momentsF, title):
   """
   Apply the function decompF to build a decomposition of the 2D
   image, and apply the function momentsF to compute the moments
   based on the decomposition.
   Print the string titleF, which describes the method.
   """

   print(title)
   if len(arg)<=1:
        print("Need 2D image name")

   else:
        print("---Read pixels from file "+ arg[1])
        input_pixels = readInput(arg[1])
        print("Number of black pixels:", input_pixels.number(1))
        elements = decomposF(input_pixels)
        print("Number of black elements:", elements.num_elem())
        if hasattr(elements,"max_pair"):
           print("Max sides:", elements.max_pair())

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
