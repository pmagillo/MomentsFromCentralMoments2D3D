
from commons import fatal_error, repeat_times
from commons import read_size_names, take_just_name

def perform_task(input, what, how, times_to_repeat=1):
   image = None
   decomp = None
   output = take_just_name(input)
   decomp_suffix = "_block3D.txt"
   
   if what.startswith("A"): # Need to build block decomposition
     from reading3D import adaptiveRead as readInput
     from spiliotis3D import extractBlocks
       
     print("---Image")
     image = readInput(input)
     if image==None: fatal_error("Wrong file content for stage A")

     print("Number of black voxels:", image.number(1))
     print("Image dimensions:", image.dimX,image.dimY,image.dimZ)
     
     print("---Block decomposition")
     decomp = repeat_times(extractBlocks, image, times_to_repeat)
     print("Number of blocks:", decomp.num_elem())
     
     if what=="A":
       name = output+decomp_suffix
       print("Saving block decomposition to "+name)
       decomp.print_me(name)
     print()

   else: # Need to read block decomposition
     from spiliotis3D import readBlocks
     if not input.endswith(decomp_suffix):
       print("Expected suffix",decomp_suffix)
       fatal_error("Wrong file format for stage B")   
     decomp = readBlocks(input)
     if decomp==None: 
       fatal_error("Wrong file content for stage B") 

   if what.endswith("BC"): # Need to compute moments
     print("---Block decomposition moments ("+how+")")
     if how=="old":
       from momentBlock3D import blockMoments
       from momentBlock3D import preprocessing
       repeat_times(preprocessing, decomp, times_to_repeat)
       moments = repeat_times(blockMoments, decomp, times_to_repeat)
     elif how=="new":
       from momentBlockNew3D import blockMoments
       from momentBlockNew3D import preprocessing
       from momentBlockNew2D import setOptimizationLevel
       setOptimizationLevel(2)
       repeat_times(preprocessing, decomp, times_to_repeat)
       moments = repeat_times(blockMoments, decomp, times_to_repeat)
     from commons3D import printMoments
     printMoments(moments)
     print("")

def perform_amortized(filename, what, how, times_to_repeat=1):
   """
   what must be "ABC" or "BC"
   how must be "old" or "new"
   """
   image = None
   decomp = None
   decomp_suffix = "_block3D.txt"
   
   #Input file contains size and list of input files
   max_side, filenames = read_size_names(filename)
   if how=="old":
     from momentBlock3D import preprocessing_once as preprocessing
   else:
     from momentBlockNew3D import preprocessing_once as preprocessing
     from momentBlockNew3D import setOptimizationLevel

   print("Stage B only once with max_side ",max_side)
   if how=="new": setOptimizationLevel(2)
   repeat_times(preprocessing, max_side, times_to_repeat)

   if what.startswith("A"): # Need to build block decomposition
     from reading3D import adaptiveRead as readInput
     from spiliotis3D import extractBlocks
   else: # Need to read block decomposition
     from spiliotis3D import readBlocks
   # Need to compute moments
   if how=="old":
      from momentBlock3D import blockMoments
      from momentBlock3D import preprocessing
   elif how=="new":
      from momentBlockNew3D import blockMoments
      from momentBlockNew3D import preprocessing
   
   # process all the input files      
   for input in filenames:
     print("Processing input file",input)
     if what.startswith("A"): # Need to build block decomposition
       print("---Image")
       image = readInput(input)
       if image==None: fatal_error("Wrong file content for stage A")

       print("Number of black pixels:", image.number(1))
       print("Image dimensions:", image.dimX,image.dimY,image.dimZ)
     
       print("---Block decomposition")
       decomp = repeat_times(extractBlocks, image, times_to_repeat)
       print("Number of blocks:", decomp.num_elem())     
       print()

     else: # Need to read block decomposition
       if not input.endswith(decomp_suffix):
         print("Expected suffix",decomp_suffix)
         fatal_error("Wrong file format for stage B")   
       decomp = readBlocks(input)
       if decomp==None: 
         fatal_error("Wrong file content for stage B") 
       print("")

     # Need to compute moments
     print("---Block decomposition moments ("+how+")")
     moments = repeat_times(blockMoments, decomp, times_to_repeat)
     from commons3D import printMoments
     printMoments(moments)
     print("")

if __name__=="__main__":
  import sys
  times_to_repeat = 1
  AMORTIZED = False
  INPUT = None
  all_stages = ["ABC","A","BC","ABonceC","BonceC"]
  all_ways = ["old","new"]
  try:
    assert len(sys.argv) in [4,5]
    indStages, indWay, indInput, indTimes = 1,2,3,4
    # stages to perform, amortized stage B or not
    STAGES = sys.argv[indStages]
    assert STAGES in all_stages
    if STAGES.find("once")>=0:
      AMORTIZED = True
      STAGES = STAGES.replace("once","")
    print("STAGES",STAGES,"with amortized stage B?",AMORTIZED)
    # in which way: old or new
    WAY = sys.argv[indWay]
    if STAGES!="A": 
      assert WAY in all_ways
      print("Mode of computation",WAY)
    # input to load
    INPUT = sys.argv[indInput]
    print("INPUT",INPUT)
    # how many repetitions
    if len(sys.argv)>indTimes:
      times_to_repeat = int(sys.argv[indTimes])
      assert times_to_repeat>0
      print("Repeating ",times_to_repeat,"times")
  except:
    print("Error in arguments. Expected arguments:")
    print("1: stages to perform (ABC,A,BC,ABonceC or BonceC)")
    print("2: way of computation (old or new) not used if stages==A")
    print("3: input file (an image for ABC or A decomposition for BC)")
    print("   if Bonce the input file contains max_side and list of input files")
    print("Last: number of repetitions (optional, default 1)")
    INPUT = None # to skip next code
    raise

  if INPUT != None:
    if not AMORTIZED:
      perform_task(INPUT, STAGES, WAY, times_to_repeat)
    else:
      perform_amortized(INPUT, STAGES, WAY, times_to_repeat)
      
 