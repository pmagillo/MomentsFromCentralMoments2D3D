"""
Alternative to image3D.py to be used in case of big 3D images.
Compact representation packing the color of 8 voxels in one byte.
The access to voxels is a bit slower.
The name of the class is the same as the name of the class 
defined in image3D.py, so they can be used in exchanged way.
"""

# Auxiliary to pack and unpack the bits
Ternes = [(0,0,0),(0,1,0),(1,0,0),(1,1,0),(0,0,1),(0,1,1),(1,0,1),(1,1,1)]
PosMasks = [1,   2,   4,   8,   16,   32,  64, 128]

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

  POS = dict(zip(Ternes,PosMasks))

  def create(self, dimX, dimY, dimZ):
      """
      Create a 3D image with dimX, dimY, dimZ voxels
      on the three sides, where all pixels are zero (black).
      The image is stored as an array of bytes, where byte
      [i,j,k] stores the color (one bit) for the eight image
      voxels [i,i+1] x [j,j+1] x [k,k+1]
      """
      assert dimX>0 and dimY>0 and dimZ>0
      # dimensions of the encoded image
      self.dimX = dimX
      self.dimY = dimY
      self.dimZ = dimZ
      # dimensions of the byte array
      self.arrayDX = (dimX+1)//2
      self.arrayDY = (dimY+1)//2
      self.arrayDZ = (dimZ+1)//2
      # data matrix encoded as a linearized array (all zero)
      self.data = [0] * (self.arrayDX*self.arrayDY*self.arrayDZ)
  
  def put(self, x,y,z, col):
    """
    Put the given color col (an integer between 0 and 255)
    into pixel of coordinates x,y,z.
    """
    assert col==0 or col==1
    assert x>=0 and x<self.dimX
    assert y>=0 and y<self.dimY
    assert z>=0 and z<self.dimZ
    IND = (x//2) + self.arrayDX*(y//2) + self.arrayDX*self.arrayDY*(z//2)
    pos_mask = self.POS[(x%2,y%2,z%2)]
    if col==0:
      self.data[IND] = self.data[IND] & (~pos_mask)
    else:
      self.data[IND] = self.data[IND] | pos_mask

  def get(self, x,y,z):
    """
    Get and return the color of the voxel of coordinates x,y,z.
    Pixels outside image are considered as background (=0).
    """
    if x<0 or x>=self.dimX or y<0 or y>=self.dimY or z<0 or z>=self.dimZ:
      return 0
    IND = (x//2) + self.arrayDX*(y//2) + self.arrayDX*self.arrayDY*(z//2)
    pos_mask = self.POS[(x%2,y%2,z%2)]
    if (self.data[IND] & pos_mask): return 1
    else: return 0

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
           if self.get(x,y,z)==color: N += 1
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
  img = Image3D(nx,ny,nz) #the image is all zeros
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


def readDimensions(filename):
  """
  Get the dimensions of the image from the binary file named filename.
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
  f.close()
  return (nx,ny,nz)

