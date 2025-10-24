#python3 -m cProfile main.py FileMpeg/dog01.txt 10

# ------- GLOBAL VARIABLES

#OPT =  optimization level: 0,1,2
#DIM =  input image dimension: 2,3
#UNA = True execute on each image with separated precomputation
#UNA = False execute on many images with one precomputation

# ------ AUX FUNCTIONS

def ripeti(funzione, argom, times = 1):
    for t in range(times):
        risultato = funzione(argom)
    return risultato  
        
def main_one_image(image_file, times_to_repeat=1, always=False):
   """
   image_file = name of image file (list of x y pairs)
   times_to_repeat = how many times the same computation
     must be repeated (to take average of times)
   always = if true, block_mom_new is computed in any case.
     if false only if the decomposition has many blocks
   """
   print("---Read pixels from file "+ image_file)
   input_pixels = readInput(image_file)
   print("Number of black pixels:", len(input_pixels))
   if DIM==2:
     max_x = max([x for x,y in input_pixels])
     max_y = max([y for x,y in input_pixels])
   elif DIM==3:
     max_x = max([x for x,y,z in input_pixels])
     max_y = max([y for x,y,z in input_pixels])
     max_z = max([z for x,y,z in input_pixels])
   print("Max x coordinate: ",max_x)
   print("Max y coordinate: ",max_y)
   if DIM==3: print("Max z coordinate: ",max_z)
   print("")
   
   print("---Tree")
   elements = ripeti(buildTree,input_pixels,times_to_repeat)
   print("Number of black leaves:", elements.num_elem())
   
   print("---Tree moments (old)")
   ripeti(tree_mom_old,elements,times_to_repeat)
   
   print("---Tree moments (new)")
   ripeti(tree_pre_new,elements,times_to_repeat)
   ripeti(tree_mom_new,elements,times_to_repeat)
   print("")
   
   print("---Block decomposition")
   elements = ripeti(extractBlocks,input_pixels,times_to_repeat)
   block_num = elements.num_elem()
   print("Number of blocks:", block_num)
   if DIM==2: print("Max sides:", elements.max_pair()) 
   elif DIM==3: print("Max sides:", elements.max_triplet())
   
   print("---Block moments (old)")
   ripeti(block_pre_old,elements,times_to_repeat)
   ripeti(block_mom_old,elements,times_to_repeat)

   print("---Block moments (new)")
   setOptimizationLevel(OPT)
   print("Optimization level: ",OPT);
   if (always or (block_num>max_x)):
     ripeti(block_pre_new,elements,times_to_repeat)
     ripeti(block_mom_new,elements,times_to_repeat)
     stampaGestione()
   else:
     print("Not computed")
   print()
   
def main_many_images(image_list, max_side, times_to_repeat=1, always=False):
   print("---Block preprocessing (new)")
   setOptimizationLevel(OPT)
   print("Optimization level: ",OPT);
   print("Only once with max_side ",max_side)
   ripeti(block_pre_once,max_side,times_to_repeat)
   
   for image_file in image_list:
     print("---Read pixels from file "+ image_file)
     input_pixels = readInput(image_file)
     print("Number of black pixels:", len(input_pixels))
     max_x = max([c[0] for c in input_pixels])
     max_y = max([c[1] for c in input_pixels])
     assert max_x<=max_side and max_y<=max_side
     print("Max x coordinate: ",max_x)
     print("Max y coordinate: ",max_y)
     if DIM==3:
       max_z = max([z for x,y,z in input_pixels])
       assert max_z<=max_side
       print("Max z coordinate: ",max_z)
     print("")
     
     print("---Tree")
     elements = ripeti(buildTree,input_pixels,times_to_repeat)
     print("Number of black leaves:", elements.num_elem())

     print("---Tree moments (old)")
     ripeti(tree_mom_old,elements,times_to_repeat)

     print("---Tree moments (new)")
     ripeti(tree_pre_new,elements,times_to_repeat)
     ripeti(tree_mom_new,elements,times_to_repeat)
     print("")

     print("---Block decomposition")
     elements = ripeti(extractBlocks,input_pixels,times_to_repeat)
     block_num = elements.num_elem()
     print("Number of blocks:", block_num)
     if DIM==2:
       print("Max sides:", elements.max_pair()) 
     else: #DIM==3
       print("Max sides:", elements.max_triplet())
     
     print("---Block moments (old)")
     ripeti(block_pre_old,elements,times_to_repeat)
     ripeti(block_mom_old,elements,times_to_repeat)

     print("---Block moments (new)")
     if (always or (block_num>max_x)):
       #No: done ripeti(block_pre_new,times_to_repeat)
       ripeti(block_mom_new,elements,times_to_repeat)
       stampaGestione()
     else:
       print("Not computed");
  

def readImageList(file_name):
  f = open(file_name,"r")
  L = f.read().split()
  f.close()
  print("Immagini:")
  for name in L: print("  "+name)
  print("")
  return L


# ------ PROGRAM

"""
Arguments on command line are:
- image dimension: 2 or 3
- optimization level: 0, 1, 2
- only for optimization level==2: "once" if one
   preprocessing stage for many images
- input image file name
- number of repetitions (opzional, default = 1) 
"""

import sys        
if __name__=="__main__":
  #print(sys.argv)
  global DIM, OPT, UNA
  ind = 3
  times_to_repeat = 1
  max_side = 0
  image = None
  UNA = True
  try:
    DIM = int(sys.argv[1])
    assert DIM in (2,3)
    print("DIM =",DIM)
    OPT = int(sys.argv[2])
    assert OPT in (0,1,2)
    print("OPT =",OPT)
    if OPT==2 and sys.argv[3]=='once': UNA = False; ind = 4
    if not UNA: print("One preprocessing for many images")
    image = sys.argv[ind]
    print("Image =",image)
    ind += 1
    if not UNA:
      max_side = int(sys.argv[ind])
      assert max_side > 2
      print("Max image side =",max_side)
      ind += 1
    if len(sys.argv)>ind:
      times_to_repeat = int(sys.argv[ind])
      assert times_to_repeat>0
      print("Repeating ",times_to_repeat,"times")
      ind += 1
    if len(sys.argv)>ind:
      print("Too many arguments")
      raise IndexError
  except:
    print("Error in arguments:")
    print("First argument must be image dimension (2 or 3)")
    print("Second argument must be optimization level (0, 1 or 2)")
    print("If optimization=2, third argument may be 'once' (optional)")
    print("Next argument must be input file")
    print("  Input file is one image,")
    print("  or a file containing a list of image names if third argument='once'")
    print("If third argument='once', next argument must be max side length")
    print("Last argument may be number of repetitions (optional, default 1)")
    DIM = None # to skip next code
    
  if DIM==2:
    from commons2D import readPixels as readInput
    from commons2D import orders, printMoments
    from quadtree import buildQuadtree as buildTree
    from momentTree2D import quadtreeMoments as tree_mom_old
    from momentTreeNew2D import preprocessing as tree_pre_new
    from momentTreeNew2D import quadtreeMoments as tree_mom_new
    from spiliotis2D import extractBlocks
    from momentBlock2D import preprocessing as block_pre_old
    from momentBlock2D import blockMoments as block_mom_old
    from momentBlockNew2D import preprocessing as block_pre_new
    from momentBlockNew2D import preprocessing_once as block_pre_once
    from momentBlockNew2D import blockMoments as block_mom_new
    from momentBlockNew2D import setOptimizationLevel, stampaGestione
  elif DIM==3:
    from commons3D import readCubes as readInput
    from commons2D import orders, printMoments
    from octree import buildOctree as buildTree
    from momentTree3D import octreeMoments as tree_mom_old
    from momentTreeNew3D import preprocessing as tree_pre_new
    from momentTreeNew3D import octreeMoments as tree_mom_new
    from spiliotis3D import extractBlocks
    from momentBlock3D import preprocessing as block_pre_old
    from momentBlock3D import blockMoments as block_mom_old
    from momentBlockNew3D import preprocessing as block_pre_new
    from momentBlockNew3D import preprocessing_once as block_pre_once
    from momentBlockNew3D import blockMoments as block_mom_new
    from momentBlockNew3D import setOptimizationLevel, stampaGestione
  print("DIM=",DIM,"OPT=",OPT,"UNA=",UNA,"image=",image,"times_to_repeat=",times_to_repeat)
  if DIM in (2,3):
    if UNA: 
       main_one_image(image, times_to_repeat, always=True)
    else:
       images = readImageList(image)
       main_many_images(images, max_side, times_to_repeat, always=True)
