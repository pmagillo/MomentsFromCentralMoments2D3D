
"""
The quadtree encodes a 2D image of size LxL where
L = 2^E for some natural E.
The pixel of minimum coordinates in the image is (0,0) 
and the pixel with largest coordinates is (L-1,L-1).
The root node of the quadtree covers the entire image.
Let C=(cx,xy) be the center of the current node.
The area of a parent node is subdivided among its four
children in this way:
- child 0: x<cx, y<cy
- child 1: x>cx, y<cy
- child 2: x<cx, y>cy
- child 3: x>cx, y>cy
Each node is identified by a location code, a string on the
alphabet 0,1,2,3.
The location code of the root is the empty sequence.
The location code of another node is equal to the location code
of its parent plus a last element equal to the position of
this node as a child if its parent.
The quadtree stores:
- the lengths of two sides of the original image along x,y
- side of the square domain of the quadtree (power of two)
- exponent of the power of two
- set of leaves implemented as a dictionary with key the location code
Each node (representing a square) stores:
- the coordinates of an anchor point (xmin,ymin)
- the coordinates of the center (xcen,ycen)
- the exponent whose power of two gives the side of the square
- the color (0 or 1)
"""

def parent(loc_code):
  """
  Return the location code of the parent of the node with the 
  given location code. The location code of the parent node
  is obtained by deleting the last element.
  """
  L = len(loc_code)
  if L==0:  return None #the root has no parent
  return tuple(loc_code[0:L-1])

def children(loc_code):
  """
  Return a list of the location codes of the eight children
  of the node with the given location code. The location codes
  of the children are obtained by appending one of the numbers 
  0,1,2,3 at the end of the given code.
  """
  return [tuple(loc_code+[i]) for i in range(4)]

class QTR_Node:
  """
  An object of this class represents a node within a quadtree.  
  """
  def __init__(self, x0,y0, expon=0, color=0):
    """
    The constructor creates a node with pixel of minimum 
    coordinates (x0,y0), consisting of K x K squares,
    where K=2^expon, and having the given color.
    """
    self.xmin = x0
    self.ymin = y0
    self.exponent = expon
    self.color = color
    if expon==0:
      self.xcen = x0
      self.ycen = y0
    else:
      delta = 2**(expon-1)-0.5
      self.xcen = x0+delta
      self.ycen = y0+delta
      
  def side(self):
    return 2**self.exponent

  def __str__(self):
      S = str()
      S = S+"Node "+str((self.xmin,self.ymin))+" -- "
      S = S+str((self.xmin+self.side()-1,self.ymin+self.side()-1))
      #S = S +"; colore "+str(self.color)
      return S
 
def power_of_two(L1, L2):
  """
  Return the exponent E of the minimum power of two
  which is greater then or equal to L1 and L2.
  L1, L2 must be positive integers (they are the
  lengths of the three sides of the iumage).
  """
  L = 1
  E = 0
  while L<L1 or L<L2:
    L *= 2
    E += 1
  return E
  
class QTR_Tree:
  """
  An object of this class is a quadtree representing an image 
  """
  def __init__(self, SX = 1, SY = 1):
    """
    The constructor creates a quadtree representing a 2D image
    of SX x SY pixels, and all pixels are white.
    The side of the area covered by the quadtree is the
    smallest power of two which is greater of equal to the two
    sides of the given image.
    The quadtree stores:
    -size_x, size_y: dimensions of the image
    -exponent,side: dimension of the quadtree domain 
      (larger or equal w.r.t. to the previous ones),
      and as value
    - leaves: dictionary of quadtree leaves (key=location
      code and associated value=node).
    """
    #dimensions of the image
    self.x_side = SX
    self.y_side = SY
    #exponent E such that the side of the square covered
    #by the quadtree is 2^E
    self.exponent = power_of_two(self.x_side, self.y_side)
    self.side = 2**self.exponent
    #empty location code denotes the root, which is the
    #only node of this quadtree
    root = ()
    #dictionary containing all the leaves of the quadtree,
    #the key is the location code and the value is the node
    #with color, 0=while, 1=black
    #center = (self.side-1)*0.5
    #self.leaves = {root: QT_Node(center,center, self.exponent, color=0)}
    self.leaves = {root: QTR_Node(0,0, self.exponent, color=0)}

  def node_side(self,loc_code):
    """
    Return the length of side of the square covered by 
    the node with given location code.
    """
    return 2**(self.exponent-len(loc_code))
      
  def min_pixel(self,loc_code):
    """
    Return the (x,y,z) coordinates of the pixel with
    minimum coordinates inside the node with given
    location code.
    """
    s = 2**self.exponent
    x,y = 0,0
    for e in loc_code:
      s = s//2
      if e in (1,3): x += s
      if e in (2,3): y += s
    return (x,y)
  
  def code_for_pixel(self,x,y):
    """
    Return the location code of a node containing
    just one pixel with the given coordinates.
    """
    code = []
    x0,y0 = 0,0
    middle = 2**(self.exponent-1) # NON VA self.side//2
    while True:
      digit = 0
      if x>=x0+middle:
         digit += 2
         x0+=middle
      if y>=y0+middle:
         digit += 1
         y0+=middle
      code.append(digit)
      if middle==1: break
      else: middle = middle//2
    return tuple(code)

  def num_elem(self):
    """
    Return the number of black leaves of this quadtree.
    """
    black_leaves = [C for C in self.leaves if self.leaves[C].color==1]
    return len(black_leaves)

  def print_me(self, filename):
    F = open(filename,'w')
    for C in self.leaves:
       if self.leaves[C].color==1:
          F.write(str(self.leaves[C])+'\n')
    F.close()
    
def stampa(Q):
  print("Dimensioni ",Q.x_side, Q.y_side)
  print("Esponente ",Q.exponent)

  print("Foglie:")
  for aa,bb in zip(Q.leaves.keys(),Q.leaves.values()):
      print("  ",aa," --> ",bb)   

  NERI = [C for C in sorted(Q.leaves) if Q.leaves[C].color==1]
  BIAN = [C for C in sorted(Q.leaves) if Q.leaves[C].color==0]
  print("Foglie bianche")
  for C in BIAN:
    node = Q.leaves[C]
    #print("Nodo ",C," min coord ",Q.min_pixel(C), " lato ",Q.node_side(C))
    print("Nodo ",C," min coord ",node.xmin,node.ymin," lato ",node.side())
  print("Foglie nere")
  for C in NERI:
    node = Q.leaves[C]
    #print("Nodo ",C," min coord ",Q.min_pixel(C), " lato ",Q.node_side(C))
    print("Nodo ",C," min coord ",node.xmin,node.ymin," lato ",node.side())

def quadtreeBuild(SX,SY, black_pixels=[]):
  """
  Build the quadtree for the given 2D image
  which covers the square [0,SX-1] x [0,SY-1]
  and whose black pixels are contained in the list
  black_pixels, the other pixels are white.
  """
  Q = QTR_Tree(SX, SY)
  Qside = 2**Q.exponent
  #print("Costruisco quadtree di lato ",Qside)
  #keep a queue of location codes that are candidate to
  #be merged into a parent node
  candid = []
  #remove the root and add all individual black pixels as leaves
  del Q.leaves[()]
  #print("Appena cancellata root, foglie=",Q.leaves)
  for x,y in black_pixels:
    C = Q.code_for_pixel(x,y)
    Q.leaves[C] = QTR_Node(x,y)
    Q.leaves[C].color = 1
    if C[-1]==0: candid.append(parent(C))
    #if (len(Q.leaves)%100)==0: print(" Ora ci sono ",len(Q.leaves)," nodi")

  #print("Foglie messe (solo le nere), adesso provo a fondere")
  
  #scan the queue and try to merge nodes
  I = 0
  while I<len(candid):
    C = candid[I]
    I += 1
    #print("Provo a fondere i figli di ",C) ###C tupla loc code
    present = [(C+(i,) in Q.leaves) for i in range(4)]
    #print("   Tutti presenti",present)
    #qqq = input("Un tasto")
    if False in present: continue
    #se ci sono tutti, sono tutti neri
    children = [Q.leaves[C+(i,)] for i in range(4)]
    #print("   Allora fondo")
    #the eight children have equal color, so merge them
    for i in range(4): del Q.leaves[C+(i,)]
    #determina punto di aggancio come baricentro
    xx = children[0].xmin #0.25*sum([children.x for i in range(4)])
    yy = children[0].ymin #0.25*sum([children.y for i in range(4)])
    Q.leaves.update({C: QTR_Node(xx,yy, children[0].exponent+1, 1 )})
    #print("Ho messo nodo con location code", C)
    #print("Appena aggiornato, foglie=")
    #for aa,bb in zip(Q.leaves.keys(),Q.leaves.values()): print(aa,bb)   
    if C[-1]==0: candid.append(parent(C))
  return Q

def Funziona_ma_lento_quadtreeBuild(SX,SY, black_pixels=[]):
  """
  Build the quadtree for the given 2D image
  which covers the square [0,SX-1] x [0,SY-1]
  and whose black pixels are contained in the list
  black_pixels, the other pixels are white.
  """
  Q = QTR_Tree(SX, SY)
  Qside = 2**Q.exponent
  #print("Costruisco quadtree di lato ",Qside)
  #keep a queue of location codes that are candidate to
  #be merged into a parent node
  candid = []
  #remove the root and add all individual pixels, white
  del Q.leaves[()]
  #print("Appena cancellata root, foglie=",Q.leaves)
  for x in range(0,Qside):
     for y in range(0,Qside): 
          C = Q.code_for_pixel(x,y)
          #print(" metto ",x,y, " codice ",C)
          Q.leaves.update({C: QTR_Node(x,y)})
          if (x,y) in black_pixels: Q.leaves[C].color = 1 #****************
          if C[-1]==0: candid.append(parent(C))
          #if (len(Q.leaves)%100)==0: print(" Ora ci sono ",len(Q.leaves)," nodi")
  #print("Appena messi tutti bianchi, foglie=")
  #for aa,bb in zip(Q.leaves.keys(),Q.leaves.values()): print(aa,bb)
  #print(" e candidati ",candid)
  #print("Num foglie:",  len(Q.leaves))
  #print("Qtree edge:",Qside)
  #assert len(Q.leaves)==(Qside*Qside)
  #replace the color of all pixels in the list
  """
  **************** ALTERNATIVO A SOPRA
  for x,y in black_pixels:
    C = Q.code_for_pixel(x,y)
    Q.leaves[C].color = 1
  """
  #print("Appena messi anche i neri, foglie=")
  #for aa,bb in zip(Q.leaves.keys(),Q.leaves.values()): print(aa,bb)
  #assert len(Q.leaves)==(Qside*Qside)

  #print("Foglie messe, adesso provo a fondere")
  
  #scan the queue and try to merge nodes
  I = 0
  while I<len(candid):
    C = candid[I]
    I += 1
    #print("Provo a fondere i figli di ",C) ###C tupla loc code
    present = [(C+(i,) in Q.leaves) for i in range(4)]
    #qqq = input("Un tasto")
    if False in present: continue
    children = [Q.leaves[C+(i,)] for i in range(4)]
    child_col = [children[i].color for i in range(4)]
    if len(set(child_col))==1:
      #print("   Allora fondo")
      #the eight children have equal color, so merge them
      for i in range(4): del Q.leaves[C+(i,)]
      #determina punto di aggancio come baricentro
      xx = children[0].xmin #0.25*sum([children.x for i in range(4)])
      yy = children[0].ymin #0.25*sum([children.y for i in range(4)])
      Q.leaves.update({C: QTR_Node(xx,yy, children[0].exponent+1,child_col[0])})
      #print("Ho messo nodo con location code", C)
      #print("Appena aggiornato, foglie=")
      #for aa,bb in zip(Q.leaves.keys(),Q.leaves.values()): print(aa,bb)   
      if C[-1]==0: candid.append(parent(C))
  return Q

def buildQuadtree(black_pixels):
  """
  Build and return the quadtree for the 2D image 
  given as list of black squares.
  """
  maxX = 1+max([c[0] for c in black_pixels])
  maxY = 1+max([c[1] for c in black_pixels])
  return quadtreeBuild(maxX, maxY, black_pixels)

#---------------------MAIN-----------------------

from commons2D import readPixels
import sys

def main(arg):
    if len(arg)<=1:
         print("Need 2D image name")

    else:
         print("---Read pixels from file "+ arg[1])
         input_pixels = readPixels(arg[1])
         QT = buildQuadtree(input_pixels)
         print("Number of black nodes: ",QT.num_elem())
         QT.print_me("builtqtree_out.txt")
         print("Nodi stampati su builtqtree_out.txt")

if __name__ == "__main__":
   main(sys.argv)
