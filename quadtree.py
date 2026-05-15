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
      return str(self.xmin)+' '+str(self.ymin)+' '+str(self.exponent)
 
def power_of_two(L1, L2):
  """
  Return the exponent E of the minimum power of two
  which is greater then or equal to L1 and L2.
  L1, L2 must be positive integers (they are the
  lengths of the three sides of the image).
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
      if e in (1,3): y += s
      if e in (2,3): x += s
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
    F.write(str(self.x_side)+' '+str(self.y_side)+'\n')
    F.write(str(self.exponent)+'\n')
    F.write(str(self.side)+'\n')
    for C in self.leaves:
       if self.leaves[C].color==1:
          s = ""
          for digit in C: s = s+str(digit)
          F.write(str(s)+'\n')
          #F.write(str(s)+' '+str(self.leaves[C])+'\n')
    F.close()

def readCode(thestring):
  return tuple([int(c) for c in thestring])
  

def readQuadtree(filename):
    F = open(filename,'r')
    SX, SY = [int(p) for p in F.readline().split()]
    newtree = QTR_Tree(SX, SY)
    newtree.exponent = [int(p) for p in F.readline().split()][0]
    newtree.side = [int(p) for p in F.readline().split()][0]
    newtree.leaves = dict()
    for L in F:
      parts = L.split()
      C = readCode(parts[0]) #location code
      x0, y0 = newtree.min_pixel(C)
      side = newtree.node_side(C)
      expon = newtree.exponent - len(C)
      newtree.leaves[C] = QTR_Node(x0,y0, expon, 1)
    F.close()
    return newtree


def buildQuadtree(IMG):
  """
  Build the quadtree for the given 2D image.
  """
  Q = QTR_Tree(IMG.dimX, IMG.dimY)
  Qside = 2**Q.exponent
  #keep a queue of location codes that are candidate to
  #be merged into a parent node
  candid = []
  #remove the root and add all individual black pixels as leaves
  del Q.leaves[()]
  for x in range(IMG.dimX):
    for y in range(IMG.dimY):
      if IMG.get(x,y)==1:
        C = Q.code_for_pixel(x,y)
        Q.leaves[C] = QTR_Node(x,y)
        Q.leaves[C].color = 1
        if C[-1]==0: candid.append(parent(C))
  
  #scan the queue and try to merge nodes
  I = 0
  while I<len(candid):
    C = candid[I]
    I += 1
    present = [(C+(i,) in Q.leaves) for i in range(4)]
    if False in present: continue
    #all present implies all black
    children = [Q.leaves[C+(i,)] for i in range(4)]
    #the eight children have equal color, so merge them
    for i in range(4): del Q.leaves[C+(i,)]
    #anchor point
    xx = children[0].xmin
    yy = children[0].ymin
    Q.leaves.update({C: QTR_Node(xx,yy, children[0].exponent+1, 1 )})
    if C[-1]==0: candid.append(parent(C))
  return Q


#---------------------MAIN-----------------------

from reading2D import adaptiveRead as readInput
import sys

def main(arg):
    if len(arg)<=1:
         print("Need 2D image name")

    else:
         print("---Read pixels from file "+ arg[1])
         input_image = readInput(arg[1])
         QT = buildQuadtree(input_image)
         print("Number of black leaves: ",QT.num_elem())
         QT.print_me("builtqtree_out.txt")
         print("Leaves saved to builtqtree_out.txt")

if __name__ == "__main__":
   main(sys.argv)
