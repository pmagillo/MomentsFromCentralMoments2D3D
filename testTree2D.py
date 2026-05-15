from commons import fatal_error, repeat_times
from commons import read_size_names, take_just_name
   
def perform_task(input, what, how, times_to_repeat=1):
   image = None
   tree = None
   output = take_just_name(input)   
   tree_suffix = "_tree2D.txt"
   
   if what.startswith("A"): # Need to build tree
     from reading2D import adaptiveRead as readInput
     from quadtree import buildQuadtree as buildTree
       
     print("---Image")
     image = readInput(input)
     if image==None: fatal_error("Wrong file content for stage A")

     print("Number of black pixels:", image.number(1))
     print("Image dimensions:", image.dimX,image.dimY)
     
     print("---Tree")
     tree = repeat_times(buildTree, image, times_to_repeat)
     print("Number of black leaves:", tree.num_elem())
     
     if what=="A":
       name = output+tree_suffix
       print("Saving tree to "+name)
       tree.print_me(name)
     print()

   else: # Need to read tree
     from quadtree import readQuadtree as readTree
     if not input.endswith(tree_suffix):
       print("Expected suffix",tree_suffix)
       fatal_error("Wrong file format for stage B")   
     tree = readTree(input)
     if tree==None: 
       fatal_error("Wrong file content for stage B") 

   if what.endswith("BC"): # Need to compute moments
     print("---Tree moments ("+how+")")
     if how=="old":
       from momentTree2D import quadtreeMoments as treeMoments
       moments = repeat_times(treeMoments, tree, times_to_repeat)
     elif how=="new":
       from momentTreeNew2D import quadtreeMoments as treeMoments
       from momentTreeNew2D import preprocessing
       repeat_times(preprocessing, tree, times_to_repeat)
       moments = repeat_times(treeMoments, tree, times_to_repeat)
     from commons2D import printMoments
     printMoments(moments)
     print("")

if __name__=="__main__":
  import sys
  times_to_repeat = 1
  INPUT = None
  try:
    assert len(sys.argv) in [4,5]
    indStages, indWay, indInput, indTimes = 1,2,3,4
    # stages to perform
    STAGES = sys.argv[indStages]
    assert STAGES in ["ABC","A","BC"]
    print("STAGES",STAGES)
    # in which way: old or new
    WAY = sys.argv[indWay]
    if STAGES!="A": 
      assert WAY in ["old","new"]
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
    print("1: stages to perform (ABC or A or BC)")
    print("2: way of computation (old or new) not used if stages==A")
    print("3: input file (an image for ABC or A, decomposition for BC)")
    print("Last: number of repetitions (optional, default 1)")
    INPUT = None # to skip next code
    #raise

  if INPUT != None:
     perform_task(INPUT, STAGES, WAY, times_to_repeat)
