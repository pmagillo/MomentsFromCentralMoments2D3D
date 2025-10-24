from commons3D import orders
from octree import OCT_Tree

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


def prova(xx,yy,zz,lato):
  print("Monenti di "+str((xx,yy,zz))+" lato "+str(lato))
  #(0,0) lato 1 ok, lato 2 ok
  #(1,2) lato 1 ok, lato 2 ok
  for p in range(4):
    for q in range(4):
      for r in range(4):
        if p+q>3: continue
        mom = factorG(p, xx, lato)*factorG(q, yy, lato)*factorG(r, zz, lato)
        print("Momento ",p,q,"   = ", mom)
  
  
def octreeMoments(OT):
  """
  Compute moments of order up to 3 from the image,
  that is encoded in the octree OT.
  Return a dictionary where key is the triplet
  (p,q,r) and value is the moment m_{p,q,r}
  """
  # orders is the global variable imported from commons3D
  MM = {key:0 for key in orders}
  for node in OT.leaves:
    node = OT.leaves[node] 
    if node.color==0: continue
    x,y,z = node.xmin, node.ymin, node.zmin
    L = node.side()
    for p,q,r in orders:
       MM[(p,q,r)] += ( factorG(p, x, L)*factorG(q, y, L)*factorG(r, z, L) )
  return MM


#-------------------MAIN-------------------

from octree import OCT_Tree, buildOctree
from commons3D import main
import sys

if __name__ == "__main__":
   main(sys.argv, buildOctree, None, octreeMoments, "====3D Octree, old method.")
