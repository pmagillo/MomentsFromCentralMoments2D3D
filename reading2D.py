"""
Functions to read a 2D image from standard image formats.
Otherwise a 3D image can be read from binary format with the
function readBinary defined in image2D.py
A unified function is provided, named adaptiveRead, which 
determines the way to read a file based on the file extension.
"""

from image2D import Image2D

#If PIL and numpy not installed, comment from here
from PIL import Image
import numpy

def isBinary(img_array):
  if len(img_array.shape)!=2: return False
  width,height = img_array.shape
  values = set([img_array[x][y] for x in range(width) for y in range(height)])
  return len(values)==2, min(values), max(values)

def toBinary(img_array):
  if len(img_array.shape)==2: 
    two, mini, maxi = isBinary(img_array)
    if two: return img_array
    med = max(1,(mini+maxi)//2)
    print("Greyscale image, converting to binary with threshold",med)
    return fromGraytoBinary(img_array, threshold=med)
  else:
    print("Color image, converting to binary")
    return fromRGBtoBinary(img_array)

def fromRGBtoBinary(img_array, threshold=128):
  width,height,colors = img_array.shape
  assert colors==3 or colors==4 # RGB or RGBA
  output = numpy.zeros((width,height), dtype=numpy.uint8)
  for x in range(width):
    for y in range(height):
      val = sum(img_array[x][y][0:3])/3
      if val>=threshold: output[x][y] = 1
  return output

def fromGraytoBinary(img_array, threshold=1):
  width,height = img_array.shape
  print("Image size",width,height)
  output = numpy.zeros((width,height), dtype=numpy.uint8)
  for x in range(width):
    for y in range(height):
      val = img_array[x][y]
      if val>=threshold: output[x][y] = 1
  return output

  
def readImageFormat(filename, verbose=False):
  # Read png image into a pillow image
  img = Image.open(filename)
  # Convert it to numpy array
  img = numpy.asarray(img)
  # Make it binary
  img = toBinary(img)
  # Create image in our format
  # Need to rotate because our images are the other way
  # and to swap 0/1 because in our images
  # 1=black(foreground) and 0=white(background)
  RUOTA = True
  width,height = img.shape
  total_pixels = width*height
  foreground_pixels = 0
  if RUOTA: IMG = Image2D(height,width)
  else: IMG = Image2D(width,height)
  for x in range(width):
    for y in range(height):
      if img[x][y]==0: 
         if RUOTA:  IMG.put(y,x, 1)
         else:      IMG.put(x,y, 1)
         foreground_pixels += 1
  if verbose:
    if RUOTA: print("Total pixels:",height,"x",width,"=",total_pixels)
    else:     print("Total pixels:",width,"x",height,"=",total_pixels)
    print("Foreground pixels:",foreground_pixels)
    print("Percent of foreground pixels:",foreground_pixels*100/total_pixels)
  return IMG

#If PIL and numpy not installed, comment up to here

def adaptiveRead(filename):
  dot = filename.rfind(".")
  if dot==-1: return None
  extension = (filename[dot:]).lower()
  if extension==".png" or extension==".jpg" or extension==".jpeg":
    return readImageFormat(filename)
  elif extension ==".binary":
    from image2D import readBinary
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
  print("Dimensions of read image:",image.dimX,image.dimY)
