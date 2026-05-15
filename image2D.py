class Image2D:
  """
  Class for a 2D binary image. There are only two colors:
  black=1 and white=0.
  """
  def __init__(self, a,b):
    """
    Create a 3D image with dimensions a,b,c on the three sides.
    """
    self.create(a,b)
  
  def create(self, dimX, dimY):
      """
      Create a 2D image with dimX, dimY pixels
      on the two sides, where all pixels are zero (black).
      """
      assert dimX>0 and dimY>0
      self.dimX = dimX
      self.dimY = dimY
      self.cells = [0 for i in range(dimX)]
      for x in range(dimX):
        self.cells[x] = [0 for i in range(dimY)]
      #BREVE MA FUNZIONA? self.cells = [[0]*dimY] * dimX
      
  def put(self, x,y, col):
    """
    Put the given color col (an integer between 0 and 255)
    into pixel of coordinates x,y,z.
    """
    assert col==0 or col==1
    assert x>=0 and x<self.dimX
    assert y>=0 and y<self.dimY
    self.cells[x][y] = col

  def get(self, x,y):
    """
    Get and return the color of the pixel of coordinates x,y,z.
    Pixels outside image are considered as background (=0).
    """
    #assert x>=0 and x<self.dimX
    #assert y>=0 and y<self.dimY
    if x<0 or x>=self.dimX or y<0 or y>=self.dimY:
      return 0
    return self.cells[x][y]

  def pixels(self, color):
    """
    Return the pixels of the image with given color
    as a list of triplets (x,y,z).
    """
    L = []
    for x in range(self.dimX):
      for y in range(self.dimY):
         if self.cells[x][y]==color:
             L.append( (x,y) )
    return L

  def number(self, color):
    """
    Return the number of the image with given color
    """
    N = 0
    for x in range(self.dimX):
      for y in range(self.dimY):
         if self.cells[x][y]==color: N += 1
    return N

  def num_elem(self):
    return self.number(1)

def readBinary(filename):
  """
  Load the image from the binary file named filename.
  The file contains a header followed by the data.
  The header is 4 byte encoding the image dimensions 
  (one unsigned 2 byte each (little-endian byte order).
  Voxel data is 1 byte per voxel.
  """
  f = open(filename,"rb")
  content = f.read()
  f.close()
  #decode dimensions (on 2 bytes each)
  nx = int.from_bytes(content[0:2], byteorder='little', signed=False)
  ny = int.from_bytes(content[2:4], byteorder='little', signed=False)
  print("Dimensions",nx,ny)
  #decode the data
  img = Image2D(nx,ny)
  i = 4
  for x in range(nx):
    for y in range(ny):
        val = int.from_bytes(content[i:i+1], byteorder='little', signed=False)
        img.put(x,y, val)
        i += 1
  return img

 
def writeBinary(img, filename):
  """
  Write the image to the binary file named filename.
  The file contains a header followed by the data.
  The header is 4 byte encoding the image dimensions 
  (one unsigned 2 byte each (little-endian byte order).
  Voxel data is 1 byte per voxel.
  """
  f = open(filename,"wb")
  #encode dimensions (on 2 bytes each)
  print("Dimensions", img.dimX, img.dimY)
  for D in (img.dimX, img.dimY):
    f.write( D.to_bytes(2, byteorder='little', signed=False) )
  #encode the data
  for x in range(img.dimX):
    for y in range(img.dimY):
       b = img.get(x,y).to_bytes(1, byteorder='little', signed=False)
       f.write(b)
  f.close()

