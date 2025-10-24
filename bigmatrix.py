
class  IntMatrix:
   """
   This data structure stores a 2-dimensional array of big integers. 
   The matrix is initialized as all -1 (meant as undefined),
   while meaningful values are supposed to be non-negative integers.

   The matrix is designed to store sums of powers in it.
   For k=0..maxK and n=1..maxN, matrix[k][n] contains the value of S_k(n), 
   where k is the order of moments to be computed, and n is a coordinate
   inside the image.
   This is defined as S_k(n) = n^k + (n-1)^k + (n-2)^k + ... 2^k + 1^k.
   """
   def __init__(self, max1, max2):
     # matrix of natural integer values, initially -1 (undefined)
     self.matrix = [ [-1 for j in range(max2)] for i in range(max1) ]
     # dimensions of the matrix
     self.dimX = max1
     self.dimY = max2

   def printMatrix(self, s=""):
      print("Matrix of dimension ",self.dimX," x ",self.dimY)
      for x in range(self.dimX):
         for y in range(self.dimY):
            print(s+"_"+str(x)+"("+str(y)+")= "+str(self.matrix[x][y]))

class PowerMatrix (IntMatrix):
   """
   Matrix storing S_k(n) for k=0...maxK, n=1..maxN
   """
   def __init__(self, maxK, maxN):
     IntMatrix.__init__(self, maxK+1,maxN+1)
     #print("Making matrix, maxK=",maxK,"maxN=",maxN, "counting=", counting)
     #set S_0(n) = n for all n 
     self.matrix[0][0] = 0
     for n in range(1, maxN+1):
        self.matrix[0][n] = n
     #set S_k(1) = 1 for all k   
     for k in range(1, maxK+1):
        self.matrix[k][0] = 0
        self.matrix[k][1] = 1
     #compute other values
     for n in range(2, maxN+1):
        p = 1; #init p as n^0 = 1
        for k in range(1, maxK+1):  
           p *= n # update p to n^k
           self.matrix[k][n] = p + self.matrix[k][n-1]

   def valueSum(self, k, n):
      return self.matrix[k][n]
