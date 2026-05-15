"""
Some functions used by the main test programs both for
2D and 3D images.
"""

import sys

def fatal_error(message):
   print(message)
   sys.exit(1)

def repeat_times(funct, argom, times = 1):
   for t in range(times):
       result = funct(argom)
   return result  

def read_size_names(file_name):
   f = open(file_name,"r")
   names = f.read().split()
   f.close()
   size = int(names[0])
   names = names[1:]
   print("Names from file:")
   for name in names: print("  "+name)
   print("")
   return size, names

def take_just_name(filepath):
   # remove extension
   output = filepath[filepath.rfind('/')+1:] 
   # remove directory
   output = output[:output.find('.')] # eliminate extension
   return output