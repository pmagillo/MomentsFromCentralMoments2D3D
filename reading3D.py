"""
Functions to read a 3D image from binary or ascii VTK image format.
Otherwise a 3D image can be read from binary format with the
function readBinary defined in image3D.py
A unified function is provided, named adaptiveRead, which 
determines the way to read a file based on the file extension.
"""

from image3D import Image3D

def searchKeyword(f, key, elements=1):
  """
  Read lines from file f until the line starts with key
  (case insensitive). Return the line starting with key,
  or None if not found.
  The parameter key may be a list of more keyword, in that
  case we search for one of them.
  In addition, check that the number of elements in the
  split line is the same as element.
  """
  if type(key)==str: key = [key]
  while True:
    line = f.readline()
    if line=="": break #end of file
    good = [line.upper().startswith(k) for k in key]
    if True in good: break
  parts = line.split()
  if len(parts)!=elements: return None
  return line

def readVTKheader(filename):
  f = open(filename,"r")
  # first line must be "# vtk DataFile Version" followed by a number
  line = searchKeyword(f, "# VTK DATAFILE VERSION",5)
  if line==None: return None
  # next line is comment
  line = f.readline()
  # next line must be ASCII or BINARY
  line = searchKeyword(f, ["ASCII","BINARY"],1)
  if line==None: return None
  encoding = line.upper().split()[0]
  # next line must be DATASET STRUCTURED_POINTS
  line = searchKeyword(f, "DATASET STRUCTURED_POINTS",2)
  if line==None: return None
  # next line must be DIMENSIONS followed by three integers
  line = searchKeyword(f, "DIMENSIONS",4)
  if line==None: return None
  line = line.split()
  dimX, dimY, dimZ = int(line[1])-1, int(line[2])-1, int(line[3])-1
  if dimX<1 or dimY<1 or dimZ<1: return None
  # next line must be ASPECT_RATIO followed by three real numbers
  line = searchKeyword(f, "ASPECT_RATIO",4)
  if line==None: return None
  line = line.split()
  a,b,c = float(line[1]), float(line[2]), float(line[3])
  # next line must be ORIGIN followed by three real numbers
  line = searchKeyword(f, "ORIGIN",4)
  if line==None: return None
  line = line.split()
  a,b,c = float(line[1]), float(line[2]), float(line[3])
  # next line must be CELL_DATA followed by dimX*dimY*dimZ
  line = searchKeyword(f, "CELL_DATA",2)
  if line==None: return None
  line = line.split()
  voxels = int(line[1])
  assert voxels == dimX*dimY*dimZ
  # next line must be SCALARS followed by a string and "unsigned_char 1"
  line = searchKeyword(f, "SCALARS",4)
  if line==None: return None
  line = line.split()
  if line[2]!="unsigned_char" or line[3]!="1": return None
  # next line must be LOOKUP_TABLE followed by a string 
  line = searchKeyword(f, "LOOKUP_TABLE",2)
  if line==None: return None
  line = line.split()
  # return encoding and dimensions
  f.close()
  return encoding, dimX, dimY, dimZ
  
def readVTKformat(filename):
  header = readVTKheader(filename)
  if header==None: return None
  encoding, dimX, dimY, dimZ = header
  if encoding=="ASCII": return ASCIItoImage(filename, dimX, dimY, dimZ)
  else: return BINARYtoImage(filename, dimX, dimY, dimZ)

def ASCIItoImage(filename, nx,ny,nz):
  f = open(filename,"r")
  # skip header
  line="dummy"
  while line!=None and not line.startswith("LOOKUP_TABLE"):
    line=f.readline()
  if line==None: return None
  # the rest of the file is data
  content = f.read().split()
  f.close()
  img = Image3D(nx,ny,nz) # all zeros
  i = 0
  for z in range(nz):
    for y in range(ny):
      for x in range(nx):
        val = int(content[i])
        if val!=0: img.put(x,y,z, val)
        i += 1
  return img

def BINARYtoImage(filename, nx,ny,nz):
  f = open(filename,"rb")
  content = f.read()
  f.close()
  headerlength = len(content) - (nx*ny*nz)
  # skip header
  i = headerlength
  # decode the data
  img = Image3D(nx,ny,nz) # all zeros
  for z in range(nz):
    for y in range(ny):
      for x in range(nx):
        val = content[i] #int.from_bytes(content[i], byteorder='little', signed=False) 
        if val!=0: img.put(x,y,z, val)
        i += 1
  return img

def adaptiveRead(filename):
  dot = filename.rfind(".")
  if dot==-1: return None
  extension = (filename[dot:]).lower()
  if extension==".vtk":
    return readVTKformat(filename)
  elif extension ==".binary":
    from image3D import readBinary
    return readBinary(filename)
  return None

if __name__=="__main__":
  import sys
  if len(sys.argv)<2:
    print("Need input file name")
    sys.exit(1)
  image = adaptiveRead(sys.argv[1])
  if image==None: 
    print("Unknown file format")
    sys.exit(1)
  print("Dimensions of read image:",image.dimX,image.dimY,image.dimZ)
