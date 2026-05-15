class Image3D:
  """
  Class for a 3D binary image. There are only two colors:
  black=1 and white=0.
  """
  def __init__(self, a,b,c):
    """
    Create a 3D image with dimensions a,b,c on the three sides.
    """
    self.create(a,b,c)
  
  def create(self, dimX, dimY, dimZ):
      """
      Create a 3D image with dimX, dimY, dimZ pixels
      on the three sides, where all pixels are zero (black).
      """
      assert dimX>0 and dimY>0 and dimZ>0
      self.dimX = dimX
      self.dimY = dimY
      self.dimZ = dimZ
      self.cells = [0 for i in range(dimX)]
      for x in range(dimX):
        self.cells[x] = [0 for i in range(dimY)]
        for y in range(dimY):
          self.cells[x][y] = [0 for i in range(dimZ)]
      
  def put(self, x,y,z, col):
    """
    Put the given color col (an integer between 0 and 255)
    into pixel of coordinates x,y,z.
    """
    assert col==0 or col==1
    assert x>=0 and x<self.dimX
    assert y>=0 and y<self.dimY
    assert z>=0 and z<self.dimZ
    self.cells[x][y][z] = col

  def get(self, x,y,z):
    """
    Get and return the color of the voxel of coordinates x,y,z.
    Pixels outside image are considered as background (=0).
    """
    if x<0 or x>=self.dimX or y<0 or y>=self.dimY or z<0 or z>=self.dimZ:
      return 0
    #assert x>=0 and x<self.dimX
    #assert y>=0 and y<self.dimY
    #assert z>=0 and z<self.dimZ
    return self.cells[x][y][z]

  def voxels(self, color):
    """
    Return the voxels of the image with given color
    as a list of triplets (x,y,z).
    """
    L = []
    for z in range(self.dimZ):  
      for x in range(self.dimX):
        for y in range(self.dimY):
           if self.get(x,y,z)==color:
             L.append( (x,y,z) )
    return L

  def number(self, color):
    """
    Return the number of the image with given color
    """
    N = 0
    for x in range(self.dimX):
      for y in range(self.dimY):
        for z in range(self.dimZ):
           if self.cells[x][y][z]==color: N += 1
    return N


def readBinary(filename):
  """
  Load the image from the binary file named filename.
  The file contains a header followed by the data.
  The header is 6 byte encoding the image dimensions 
  (one unsigned 2 byte each (little-endian byte order).
  Voxel data is 1 byte per voxel.
  """
  f = open(filename,"rb")
  content = f.read()
  f.close()
  #decode dimensions (on 2 bytes each)
  nx = int.from_bytes(content[0:2], byteorder='little', signed=False)
  ny = int.from_bytes(content[2:4], byteorder='little', signed=False)
  nz = int.from_bytes(content[4:6], byteorder='little', signed=False)
  print("Dimensions",nx,ny,nz)
  #decode the data
  img = Image3D(nx,ny,nz) # all zeros
  i = 6
  for z in range(nz):
    for y in range(ny):
      for x in range(nx):
        val = int.from_bytes(content[i:i+1], byteorder='little', signed=False) 
        if val!=0: img.put(x,y,z, val)
        i += 1
  return img

 
def writeBinary(img, filename):
  """
  Write the image to the binary file named filename.
  The file contains a header followed by the data.
  The header is 6 byte encoding the image dimensions 
  (one unsigned 2 byte each (little-endian byte order).
  Voxel data is 1 byte per voxel.
  """
  f = open(filename,"wb")
  #encode dimensions (on 2 bytes each)
  for D in (img.dimX, img.dimY, img.dimZ):
    f.write( D.to_bytes(2, byteorder='little', signed=False) )
  #encode the data
  for z in range(img.dimZ):
    for y in range(img.dimY):
      for x in range(img.dimX):
         b = img.get(x,y,z).to_bytes(1, byteorder='little', signed=False)
         f.write(b)
  f.close()

