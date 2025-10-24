"""
Computation of moments of a 2D image on the block 
decomposition by Spiliotis and Mertzios,
with the new idea that precomputes central moments.

Optimization level 0:
Compute all central moments of rectangles with
dimX = 1...max width of a block
dimY = 1...max height of a block
and dimX>=dimY

Optimization level 1:
Compute all central moments of rectangles with
dimX = 1...max width of a block
dimY = 1...max height of a block
(given that max width>=max height, otherwise swap)

Optimization level 2:
Compute all central moments of rectangles with
dimX = 1...max width of a block
dimY = 1...min{8,max height of a block}
and the moments of larger blocks will be computed in
the traditional way
(this requires computing the bigmatrix as well).
 
Optimization for a set of images:
Compute all central moments of rectangles with
dimX = 1...given value
dimY = 1...8
and do it only once for many images;
the  the moments of larger blocks will be computed in
the traditional way
(this requires computing the bigmatrix as well).
"""

#---------------------MOMENTS-----------------------

from commons2D import orders

from bigmatrix import PowerMatrix         #APRILE
#Matrix storing precomputed sums of powers#APRILE
powers = None                             #APRILE

OPT_LEVEL = 2
LIMIT = 8

def setOptimizationLevel(level):
  assert level in [0,1,2]
  global OPT_LEVEL
  OPT_LEVEL = level

def setCentralMoments(max_side):
  """
  The argument can be a single integer or a pair of integers,
  which are the max length of a block (rectangle) in the two
  Cartesian directions of the plane.
  Compute and store all central moments of rectangles of size DX x DY
  with DX,DY in [1,max_side], DX>=DY.
  - for DX==DY==1 the only non-zero moment is mu_00 
  - otherwise the only non-zero moments are mu_00 and mu_20, mu_02
  """
  #print("Set central moments",max_side)
  # manage argument
  if type(max_side) is int:
    MoreMax,LessMax = max_side, max_side
  else:
    MoreMax,LessMax = max_side # already sorted, decreasing
  
  #print("***********LATI: ******",MoreMax,LessMax,"**************")  
  CentrMom00 = dict()
  CentrMom20 = dict()
  CentrMom02 = dict()
  # DX==DY==1: only mu_00=1
  CentrMom00[(1,1)] = 1
  CentrMom20[(1,1)] = 0
  CentrMom02[(1,1)] = 0
  
  # Compute auxiliary values: 
  # sum_half = 0.5^2 + 1.5^2 + 2.5^2 + ...
  # sum_full = 1^2   + 2^2   + 3^2   + ...
  sum_half=[0,0.25]
  sum_full=[0]
  for f in range(2,MoreMax+1):
    sum_half.append(sum_half[f-1]+((f-0.5)**2))
  for f in range(1,MoreMax+1):
    sum_full.append(sum_full[f-1]+(f**2))
  
  #print("Somme mezze", sum_half)
  #print("Somme intere", sum_full)
    
  # DY==1 and DX>1
  #print("x in [1,",MoreMax,"] e y=1")
  for edgeX in range(1,MoreMax+1):
     #print("APRILE (a) key ",(edgeX,1))
     CentrMom00[(edgeX,1)] = edgeX
     CentrMom02[(edgeX,1)] = 0
     if edgeX%2==0:
       CentrMom20[(edgeX,1)] = 2*sum_half[edgeX//2]
     else:
       CentrMom20[(edgeX,1)] = 2*sum_full[edgeX//2]

  # DY>1 and DX>=DY
  #print("y in [2,",LessMax,"] e x in [y,",MoreMax,"]")
  for edgeY in range(2,LessMax+1):
      for edgeX in range(edgeY,MoreMax+1):
         #print("APRILE (b) key",(edgeX,edgeY))
         CentrMom00[(edgeX,edgeY)] = edgeX*edgeY
         if edgeX%2==0: sum_x = sum_half
         else: sum_x = sum_full
         if edgeY%2==0: sum_y = sum_half
         else: sum_y  = sum_full
         CentrMom20[(edgeX,edgeY)] = 2*edgeY*sum_x[edgeX//2]
         CentrMom02[(edgeX,edgeY)] = 2*edgeX*sum_y[edgeY//2]
         # check
         """
         M20, M02 = 0, 0
         baric = ((edgeX-1)/2, (edgeY-1)/2)
         for i in range(edgeX):
           for j in range(edgeY):
             xi = i-baric[0]
             yj = j-baric[1]
             M20 += (xi*xi)
             M02 += (yj*yj)

         if CentrMom20[(edgeX,edgeY)] != M20:
           print("ERRORE!!!!!!!!!!! ",edgeX,edgeY,": mom_2,0 ",CentrMom20[(edgeX,edgeY)], " diverso dal vero ",M20)
           print("  20: ",2,"*",edgeY,"* sommax[",(edgeX//2),"] che vale ", sum_x[edgeX//2])

         if CentrMom02[(edgeX,edgeY)] != M02:
           print("ERRORE!!!!!!!!!!! ",edgeX,edgeY,": mom_0,2 ",CentrMom02[(edgeX,edgeY)], " diverso dal vero ",M02)
           print("  02: ",2,"*",edgeX,"* sommay[",(edgeY//2),"] che vale ", sum_y[edgeY//2]) 
           #print("  ", 2*edgeX*sum_y[edgeY//2])
         """
  return (CentrMom00, CentrMom20, CentrMom02)

def preprocessing(ibr):
  #assert isinstance(ibr,BW_BlockImage2D)

  # precompute central moments for new method
  global CC00, CC20, CC02
  # precompute matrix for traditional method
  global powers
  
  # manage optimization level
  if OPT_LEVEL==0:
    maximum = max(ibr.max_pair())
    CC00, CC20, CC02 = setCentralMoments(maximum)
  elif OPT_LEVEL==1:
    maximum = sorted(ibr.max_pair(),reverse=True)
    CC00, CC20, CC02 = setCentralMoments(maximum)
  elif OPT_LEVEL==2:
    maximum = sorted(ibr.max_pair(),reverse=True)
    maximum[1] = min(maximum[1],LIMIT)
    CC00, CC20, CC02 = setCentralMoments(maximum)  
    powers = PowerMatrix( 3, ibr.origsize )
  

def preprocessing_once(max_side):
  # precompute matrix for traditional method
  global powers
  powers = PowerMatrix( 3, max_side )

  # precompute central moments for new method
  global CC00, CC20, CC02
  if OPT_LEVEL==0:
    CC00, CC20, CC02 = setCentralMoments(max_side)
  else: # 1,2
    CC00, CC20, CC02 = setCentralMoments((max_side,LIMIT))

def blockMoments(ibr):
  """
  Compute all moments m_{p,q} for p,q>=0 and p+q<=3
  of a 2D image given as a set of blocks
  """

  # orders is the global variable imported from commons3D
  # initialize moments to be computed
  MM = {key:0 for key in orders}
  
  for p,q in orders:
        MM[(p,q)] = 0 
  #print("  num blocchi",len(ibr.block))

  global NUOVO
  global VECCHIO 
  NUOVO,VECCHIO = 0,0 #APRILE

  for b in ibr.block: # cycle on blocks

     if OPT_LEVEL>1:
       dimens = (b.x1-b.x0+1, b.y1-b.y0+1) #APRILE
       if min(dimens)>LIMIT: #APRILE faccio al modo vecchio
         #print("VECCHIO MODO",dimens,ordered)
         for p,q in orders:
           if p==0: mx = dimens[0]
           else:
             mx = powers.valueSum(p, b.x1)
             if b.x0>0:
                mx -= powers.valueSum(p, b.x0-1)
           if q==0: my = dimens[1]
           else:
             my = powers.valueSum(q, b.y1)
             if b.y0>0:
                my -= powers.valueSum(q, b.y0-1)
           if (mx or my): MM[(p,q)] += (mx*my)
         VECCHIO += 1
         continue
     #FINE APRILE   

     #print("NUOVO MODO",dimens,ordered)
     NUOVO += 1

     #print("Momento di ",b, " di ",b.pixel_num(), " pixel")
     # barycenter
     xx = 0.5*(b.x1+b.x0)
     yy = 0.5*(b.y1+b.y0)
     # retrieve central moments of block
     if (b.x1-b.x0)>=(b.y1-b.y0):
        key = (b.x1-b.x0+1, b.y1-b.y0+1)
        central00 = CC00[key]
        central20 = CC20[key]
        central02 = CC02[key]
     else:
        key = (b.y1-b.y0+1,b.x1-b.x0+1)
        central00 = CC00[key]
        central20 = CC02[key]
        central02 = CC20[key]
     #print(' chiave ',key)
     #print(" mom centr 00 20 02: ",central00,central20,central02)
     # compute moments from central ones
     m00 = central00
     m10 = xx*central00
     m01 = yy*central00
     m11 = yy*m10
     m20 = central20 + xx*m10
     m02 = central02 + yy*m01
     m30 = 3*xx*m20 -2*xx*xx*m10
     m03 = 3*yy*m02 -2*yy*yy*m01
     m21 = yy*m20
     m12 = xx*m02
     #print(" mom blocco 00 20 02: ",m00,m20,m02)
     # update image moments
     MM[(0,0)] += m00
     MM[(1,0)] += int(m10)
     MM[(0,1)] += int(m01)
     MM[(1,1)] += int(m11)
     MM[(2,0)] += int(m20)
     MM[(0,2)] += int(m02)
     MM[(3,0)] += int(m30)
     MM[(0,3)] += int(m03)
     MM[(2,1)] += int(m21)
     MM[(1,2)] += int(m12)
      
  return MM
   
# da chiamare subito dopo blockMoments
def stampaGestione():
   print("N. Blocks processed with new and with traditional way",NUOVO,VECCHIO)
   #print("Blocchi gestiti col nuovo e col vecchio",NUOVO,VECCHIO)

#---------------------MAIN-----------------------

from spiliotis2D import extractBlocks
from commons2D import main
import sys

if __name__ == "__main__":
   main(sys.argv, extractBlocks, preprocessing, blockMoments, "====2D Blocks, new method.")
   print("Blocchi gestiti col nuovo e col vecchio",NUOVO,VECCHIO)
