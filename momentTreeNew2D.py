from commons2D import orders
from quadtree import QTR_Tree

def setCentralMoments(max_side):
  """
  Compute and store all central moments of squares
  with edge equal to a power of two, from 1 to max_side.
  - for edge==1 the only non-zero moment is mu_00 
  - for edge>1 the only non-zero moments are mu_00
    and mu_20=mu_02
  """
  CentrMom0 = []
  CentrMom2 = []
  # edge 1: mu_00=1 and mu_20=0
  CentrMom0.append(1)
  CentrMom2.append(0)
  # edge 2: mu_00=4 and mu_20=1
  CentrMom0.append(4)
  CentrMom2.append(1)
  # other edges
  prev_edge = 2
  edge = 4
  sommaQ = 0.25 # (0.5)^2 
  #term = 0.5
  while edge<=max_side:
    for t in range(prev_edge+1,edge,2):
         sommaQ += ((0.5*t)**2)
    sqr_edge = edge*edge;
    CentrMom0.append(sqr_edge) # mu_00=edge^2
    CentrMom2.append(int(2*edge*sommaQ)) #CAMBIATO 14 marzo
    #print("Per ",edge," viene ",CentrMom2[-1])
    prev_edge = edge
    edge *= 2
  return (CentrMom0, CentrMom2)

def preprocessing(QT):
  global stored0, stored2 
  stored0, stored2 = setCentralMoments(QT.side)

    
def quadtreeMoments(QT):
  """
  [NEW] Compute moments of order up to 3 from the 2D image,
  that has been encoded in the quadtree QT, exploiting 
  precomputed central moments.
  Return a dictionary where key is the pair
  (p,q) and value is the moment m_{p,q}
  """
  MM = {key:0 for key in orders}
  for node in QT.leaves:
    node = QT.leaves[node] 
    if node.color==0: continue
    x,y = node.xcen, node.ycen
    #s = node.side()
    #print("Nodo con baricentro ",x,y, " indice",node.exponent)
    m00 = stored0[node.exponent]
    #assert type(m00) is int #**************
    m10 = x*m00
    #assert int(m10)==10
    m01 = y*m00
    #assert int(m01)==m01
    #
    m11 = y*m10 #x*m01
    #assert int(m11)==m11
    m20 = stored2[node.exponent] + x*m10
    #assert int(m20)==m20
    m02 = stored2[node.exponent] + y*m01
    #assert int(m02)==m02
    #    
    m12 = x*m02
    #assert int(m12)==12
    m21 = y*m20
    #assert int(m21)==m21
    #
    m30 = 3*x*m20 -2*x*x*m10
    #assert int(m30)==m30
    m03 = 3*y*m02 -2*y*y*m01
    #assert int(m03)==m03

    MM[(0,0)] += m00
    MM[(1,0)] += int(m10)
    MM[(0,1)] += int(m01)
    MM[(1,1)] += int(m11)
    MM[(2,0)] += int(m20)
    MM[(0,2)] += int(m02)
    MM[(1,2)] += int(m12)
    MM[(2,1)] += int(m21)
    MM[(3,0)] += int(m30)
    MM[(0,3)] += int(m03)
      
  return MM

#-------------------MAIN-------------------

from quadtree import QTR_Tree, buildQuadtree
from commons2D import main
import sys

if __name__ == "__main__":
   main(sys.argv, buildQuadtree, preprocessing, quadtreeMoments, "====2D Quadtree, new method.")
