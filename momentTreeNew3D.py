from commons3D import orders
from octree import OCT_Tree

def setCentralMoments(max_side):
  """
  Compute and store all central moments of cubes
  with edge equal to a power of two, from 1 to max_side.
  - for edge==1 the only non-zero moment is mu_000 
  - for edge>1 the only non-zero moments are mu_000
    and mu_200=mu_020=mu_002
  """
  CentrMom0 = []
  CentrMom2 = []
  # edge 1: mu_000=1 and mu_200=0
  CentrMom0.append(1)
  CentrMom2.append(0)
  # edge 2: mu_000=8 and mu_200=2
  CentrMom0.append(8)
  CentrMom2.append(2)
  # other edges
  prev_edge = 2
  edge = 4
  sommaQ = 0.25 # (0.5)^2 
  while edge<=max_side:
      for t in range(prev_edge+1,edge,2):
           sommaQ += ((0.5*t)**2)
      sqr_edge = edge*edge;
      CentrMom0.append(edge*sqr_edge) # mu_000=edge^3
      CentrMom2.append(int(2*sqr_edge*sommaQ))
      #print("Per ",edge," viene ",CentrMom2[-1])
      prev_edge = edge
      edge *= 2
  return (CentrMom0, CentrMom2)

def preprocessing(QT):
  global stored0, stored2
  stored0, stored2 = setCentralMoments(QT.side)
      
def octreeMoments(OT):
  """
  Compute moments of order up to 3 from the 3D image,
  that has been encoded in the octree OT, exploiting
  precomputed central moments.
  Return a dictionary where key is the triplet
  (p,q,r) and value is the moment m_{p,q,r}
  """
  stored0, stored2 = setCentralMoments(OT.side)
  MM = {key:0 for key in orders}
  for node in OT.leaves:
    node = OT.leaves[node] 
    if node.color==0: continue
    x,y,z = node.xcen, node.ycen, node.zcen
    #s = node.side()
    #print("Nodo con baricentro ",x,y, " indice",node.exponent)
    m000 = stored0[node.exponent]
    #
    m100 = x*m000
    m010 = y*m000
    m001 = z*m000
    #
    m011 = z*m010
    m101 = x*m001
    m110 = y*m100

    m200 = stored2[node.exponent] + x*m100
    m020 = stored2[node.exponent] + y*m010
    m002 = stored2[node.exponent] + z*m001
    
    m021 = z*m020
    m210 = y*m200
    m120 = x*m020
    m102 = x*m002
    m012 = y*m002
    m201 = z*m200

    m300 = 3*x*m200 -2*x*x*m100
    m030 = 3*y*m020 -2*y*y*m010
    m003 = 3*z*m002 -2*z*z*m001

    m111 = x*m011

    MM[(0,0,0)] += m000
    MM[(1,0,0)] += m100
    MM[(0,1,0)] += m010
    MM[(1,1,0)] += m110
    MM[(2,0,0)] += m200
    MM[(0,2,0)] += m020
    MM[(2,1,0)] += m210
    MM[(1,2,0)] += m120
    MM[(3,0,0)] += m300
    MM[(0,3,0)] += m030
    MM[(0,0,3)] += m003
    
    MM[(0,0,1)] += m001
    MM[(1,0,1)] += m101
    MM[(0,1,1)] += m011
    MM[(1,1,1)] += m111
    MM[(2,0,1)] += m201
    MM[(0,2,1)] += m021
    
    MM[(0,0,2)] += m002
    MM[(1,0,2)] += m102
    MM[(0,1,2)] += m012

  return MM

#-------------------MAIN-------------------

from octree import OCT_Tree, buildOctree
from commons3D import main
import sys

if __name__ == "__main__":
   main(sys.argv, buildOctree, preprocessing, octreeMoments, "====3D Octree, new method.")
