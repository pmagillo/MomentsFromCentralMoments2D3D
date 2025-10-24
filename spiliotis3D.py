"""
Decomposition of a 3D image into cuboids
by Spiliotis and Mertzios.
The procedure works slice by slice.
On each 2D slice it considers the horizontal runs of
the image and tries to merge runs of equal length
existing in consecutive rows.
Then, it tries to merge equal rectangles existing in
consecutive slices.
"""
#---------------------CLASSES-----------------------

class BW_Block:
  """
  Type for a block of voxels covering the cuboid [x0,x1] x [y0,y1] x [z0,y0],
  including endpoints.
  """
  def __init__(self, b=None):
    """
    If a block b is given, copy b into this new block.
    Otherwise initialize this new block with dummy values.
    """
    if b==None:
      self.x0, self.y0 = 0,0
      self.x1, self.y1 = 0,0
      self.z1, self.z1 = 0,0
    else:
      assert isinstance(b,BW_Block)
      self.x0 = b.x0
      self.y0 = b.y0
      self.z0 = b.z0
      self.x1 = b.x1
      self.y1 = b.y1
      self.z1 = b.z1
      
  def set(self, x0,y0,z0, x1,y1,z1):
      """
      Set the cuboid of this block as [x0,x1] x [y0,y1] x [z0,z1].
      """
      self.x0 = x0
      self.y0 = y0
      self.z0 = z0
      self.x1 = x1
      self.y1 = y1
      self.z1 = z1

  def has_equal_xy(self, b):
      """
      Return true iff this block and the block b are the same cuboid.
      """
      return self.x0==b.x0 and self.x1==b.x1 and self.y0==b.y0 and self.y1==b.y1

  def __str__(self):
      return "Block "+str((self.x0,self.y0,self.z0))+" -- "+str((self.x1,self.y1,self.z1))

  def pixel_num(self):
      """
      Return the number of voxels of this block.
      """
      return (self.x1-self.x0+1)*(self.y1-self.y0+1)*(self.z1-self.z0+1)

class BW_BlockImage3D:
  """
  Image Block Representation by Spiliotis & Mertzios.
  The image is represented as a list of blocks.
  """

  def __init__(self):
    self.block = []
    self.origsize = 0

  def add_block(self, b):
    assert isinstance(b,BW_Block)
    self.block.append(b)

  def size(self):
    return len(self.block)

  def __str__(self):
     s = "Block 3D image\n"
     for b in self.block:
       s = s + str(b) + "\n"
     return s

  def num_elem(self):
     return len(self.block)
     
  def all_blocks(self):
     return self.block

  def max_triplet(self):
     """
     Return a triplet by the max height of a block
     the max width of a block, and the max depth of a block
     """
     max_deltaX, max_deltaY, max_deltaZ = 0,0,0
     for b in self.block:
       max_deltaX = max(max_deltaX,b.x1-b.x0+1)
       max_deltaY = max(max_deltaY,b.y1-b.y0+1)
       max_deltaZ = max(max_deltaZ,b.z1-b.z0+1)
     return (max_deltaX,max_deltaY,max_deltaZ)

  def statistiche(self):
     dim_uno = 0
     dim_piccole = 0 
     soglia = 7
     for b in self.block:
       if b.x1-b.x0<soglia and b.y1-b.y0<soglia and b.z1-b.z0<soglia:
         dim_piccole += 1
       if b.x1==b.x0 or b.y1==b.y0 or b.z1==b.z0:
         dim_uno += 1
     print("Di",len(self.block),"blocchi: hanno spessore 1 in",dim_uno," e dimensioni<8 in",dim_piccole)
     print("Max dimens: ",self.max_triplet())
     

#---------------------DECOMPOSITION-----------------------


def extractSliceBlocks(IMG, SX, SY, z):
  """
  Create the blocks for the 2D slice with given z.
  """
  #print("==================FORMO BLOCCHI DELLA FETTA z=",z)
  num = 0
  # Array used to store in-progress blocks. 
  # The array is indexed on x1.
  # For all used x, either block[x]=None or block[x].x1 = x.
  temp_block = [None for i in range(SX+1)]
  # Result to be returned
  slice = BW_BlockImage3D()
  slice.origsize = max([SX,SY,z])

  for y in range(SY+1):
     start_x = 0 # start of run
     inside_run = False # are we in a run? initially no 
     for x in range(SX+1):
         #print("===Processo pixel ",(x,y,z))
         if inside_run and not ((x,y,z) in IMG):
           #print("Fine run al pixel nero ",(x-1,y,z))
           #the run ended at x-1, and certainly x>0 because initially inside_run is off
           inside_run = False
           #if there is an in-progress block ending at x-1
           if temp_block[x-1]!=None:
           #if this block starts at start_x and extends up to previous y, then extend the block
              if (temp_block[x-1].x0==start_x) and (temp_block[x-1].y1==y-1):
                #print("   aggiorno blocco ",temp_block[x-1])
                temp_block[x-1].y1 = y
              else: #write the block and overwrite it 
                #print("   salvo blocco ",temp_block[x-1])
                slice.add_block(BW_Block(temp_block[x-1]))
                temp_block[x-1].set(start_x,y,z, x-1,y,z)
                num += 1
           else: # there is not an in-progress block ending at x-1, set a new block
             temp_block[x-1] = BW_Block()
             temp_block[x-1].set(start_x,y,z, x-1,y,z)
             #print("   creo blocco ", temp_block[x-1]);
         elif (not inside_run) and ((x,y,z) in IMG):
            #print("Inizio run al pixel nero ",(x,y,z))
            # we start being inside a run 
            start_x = x
            inside_run = True
      # the row ended, if we are stil inside a run, then a block ends at x=SX 
     if inside_run:
        #print("Fine riga col nero a",(SX,y,z))
        if temp_block[SX]!=None:
          #if this block starts at start_x and extends up to previous y, then extend the block
          if (temp_block[SX].x0==start_x) and (temp_block[SX].y1==y-1):
            temp_block[SX].y1 = y
            temp_block[SX].z1 = z
            #print("   aggiorno blocco ", temp_block[SX])
          else: # write the block and overwrite it
            #print("   salvo blocco ", temp_block[SX])
            slice.add_block(BW_Block(temp_block[SX]))
            temp_block[SX].set(start_x,y,z, SX,y,z)
            num += 1
            #print("   e aggiorno ", temp_block[SX])
        else: # there is not an in-progress block ending at SX, set a new one
          temp_block[SX] = BW_Block()
          temp_block[SX].set(start_x,y,z, SX,y,z)
          #print("   creo blocco ", temp_block[SX])
  # end for y
  #now write all remaining in-progress blocks
  for x in range(SX+1):
     if temp_block[x]!=None:
        #print("salvo blocco pendente ",x, temp_block[x])
        slice.add_block(BW_Block(temp_block[x]))
        num += 1
  #print("Number of blocks in slice z=",z, ": ",slice.size(),"\n")
  return slice
  
def extractBlocks(black_cubes):
  """
  Build the decomposition into blocks from a 3D image given
  as a list of black cubes, and return it.
  """
  # Convert the image from list of black voxels to 
  # a dictionary with key=(x,y,z) and value =1
  IMG = dict([(c,1) for c in black_cubes])
  SX = max([c[0] for c in black_cubes])
  SY = max([c[1] for c in black_cubes])
  SZ = max([c[2] for c in black_cubes])
  num = 0
  # Array used to store in-progress blocks of the previous slice
  # The array is indexed on x1.
  # For all used x, either block[x]=None or block[x].x1 = x.
  temp_block = [None for i in range(SX+1)]
  # Result to be returned
  final_blocks = BW_BlockImage3D()
  final_blocks.origsize = max([SX,SY,SZ])

  for z in range(SZ+1):
      slice_blocks = extractSliceBlocks(IMG, SX, SY, z)
      #CHECK SLICE
      #checkBlocks(slice_blocks,[c for c in black_cubes if c[2]==z])
      #print("======Blocchi fetta z=",z," vanno bene")
      #print("======FUSIONE")
      for BL in slice_blocks.all_blocks():
           #print("===Guardo se posso fondere il blocco ",BL)
           x = BL.x1+1
           if temp_block[x-1]!=None:
               if temp_block[x-1].has_equal_xy(BL):
                  #print("   si: fondo con blocco ",temp_block[x-1])
                  temp_block[x-1].z1 = z
               else: #write the block and overwrite it 
                  for xx in range(BL.x0,BL.x1+1):
                     if temp_block[xx]!=None:
                        #print("   no: salvo blocco ",temp_block[xx])
                        final_blocks.add_block(BW_Block(temp_block[xx]))
                        temp_block[xx] = None
                  temp_block[x-1] = BL
                  num += 1
                  #print("   e aggiorno ", temp_block[x-1])
           else: # there is not an in-progress block ending at x-1, set a new block
               for xx in range(BL.x0,BL.x1+1):
                  if temp_block[xx]!=None:
                     #print("   no: salvo blocco ",temp_block[xx])
                     final_blocks.add_block(BW_Block(temp_block[xx]))
                     temp_block[xx] = None
               temp_block[x-1] = BL
               #print("   creo blocco ", temp_block[x-1]);
      # end for BL
      for x in range(SX+1):
        if temp_block[x]!=None and temp_block[x].z1==z-1:
          # this block finishes at z-1 and has no chance to be merged
          # with the new block that will be built on next slice
          #print("salvo blocco pendente z-1",temp_block[x])
          final_blocks.add_block(BW_Block(temp_block[x]))
          temp_block[x] = None
          num += 1   
  # end for z
  #now write all remaining in-progress blocks
  for x in range(SX+1):
     if temp_block[x]!=None:
        #print("salvo blocco pendente finale",temp_block[x])
        final_blocks.add_block(BW_Block(temp_block[x]))
        num += 1
  #print("Number of blocks = ",final_blocks.size())
  return final_blocks

def checkBlocks(ibr, img):
   """
   Take an image block representation and a 3D image (given as a
   list of black cubes), and check that the block
   representation is equal to the image.
   In order to do this, for each block, change the color in the image:
   if it was black, put it white and vice versa.
   At the end, the image must be completely white.
   Return true iff it is all white.
   """

   assert isinstance(ibr, BW_BlockImage3D)
   assert isinstance(img,list)

   print("Are ",len(img)," pixels represented by ",ibr.size()," blocks?")
   OK = True
   # read ibr and change the color in current image
   for b in ibr.block:
     for x in range(b.x0,b.x1+1):
        for y in range(b.y0,b.y1+1):
           for z in range(b.z0,b.z1+1):
             if (x,y,z) in img:
                img.remove((x,y,z))
             else:
               print("Error, block pixel ",(x,y,z), " was white")
               OK = False
               #img[(x,y,z)]=1

   # the image now must be white
   if len(img)>0:
      for c in img:
           print("Error, black pixel ",c," not in block")
           OK = False
   if OK: print("TUTTO VA BENE")
   else:
      print("ERRORE")
      sys.exit(1)


#---------------------MAIN-----------------------

from commons3D import readCubes
import sys

def main(arg):
    if len(arg)<=1:
         print("Need 3D image name")

    else:
         print("---Read cubes from file "+ arg[1])
         input_cubes = readCubes(arg[1])
         BB = extractBlocks(input_cubes)
         #checkBlocks(BB, input_cubes)
         print("Number of blocks: ",BB.num_elem())
         BB.statistiche()

if __name__ == "__main__":
   main(sys.argv)
