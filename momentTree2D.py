from commons2D import orders
from quadtree import QTR_Tree

def factorG(order, coord, side):
  if order==0:
    return side
  side2=side*side
  if order==1:
    return side*coord+0.5*(side2-side)
  side3=side2*side
  coord2=coord*coord
  if order==2:
    return side*coord2 + coord*(side2-side)+(2*side3-3*side2+side)/6
  # here order==3
  side4=side3*side
  coord3=coord2*coord
  return coord3*side +1.5*coord2*(side2-side) +coord*side3 -0.5*coord*(3*side2-side) +0.25*(side4-2*side3+side2)


def prova(xx,yy,lato):
  print("Monenti di ("+str(xx)+","+str(yy)+") lato "+str(lato))
  #(0,0) lato 1 ok, lato 2 ok
  #(1,2) lato 1 ok, lato 2 ok
  for p in range(4):
    for q in range(4):
      if p+q>3: continue
      mom = factorG(p, xx, lato)*factorG(q, yy, lato)
      print("Momento ",p,q,"   = ", mom)
  

def quadtreeMoments(QT):
  """
  Compute moments of order up to 3 from the image,
  that is encoded in the quadtree QT.
  Return a dictionary where key is the pair
  (p,q) and value is the moment m_{p,q}
  """
  # orders is the global variable imported from commons2D
  MM = {key:0 for key in orders}
  for node in QT.leaves:
    node = QT.leaves[node] 
    if node.color==0: continue
    x,y = node.xmin, node.ymin
    L = node.side()
    for p,q in orders:
       MM[(p,q)] += int( factorG(p, x, L)*factorG(q, y, L) )
  return MM


#-------------------MAIN-------------------

from quadtree import QTR_Tree, buildQuadtree
from commons2D import main
import sys

if __name__ == "__main__":
   #provaPrecalcoli()
   main(sys.argv, buildQuadtree, None, quadtreeMoments, "====2D Quadtree, old method.")
