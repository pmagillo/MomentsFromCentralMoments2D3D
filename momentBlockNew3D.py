"""
Computation of moments of a 3D image on the block 
decomposition by Spiliotis and Mertzios,
with the new idea that precomputes central moments.

Optimization level 0:
Compute all central moments of cuboids with
dimX = 1...max width of a block
dimY = 1...max height of a block
dimZ = 1...max length of a block
and dimX>=dimY>=DimZ

Optimization level 1:
Compute all central moments of rectangles with
dimX = 1...max width of a block
dimY = 1...max height of a block
dimZ = 1...max length of a block
(given that max width>=max height>=max length, otherwise swap)

Optimization level 2:
Compute all central moments of rectangles with
dimX = 1...max width of a block
dimY = 1...min{8,max height of a block}
dimZ = 1...min{2,max length of a block}
and the moments of larger blocks will be computed in
the traditional way
(this requires computing the bigmatrix as well).

Optimization for a set of images:
Compute all central moments of rectangles with
dimX = 1...given value
dimY = 1...8
dimZ = 1...2
and do it only once for many images;
the  the moments of larger blocks will be computed in
the traditional way
(this requires computing the bigmatrix as well).
"""

#---------------------MOMENTS-----------------------

from commons3D import orders

from bigmatrix import PowerMatrix         #APRILE
#Matrix storing precomputed sums of powers#APRILE
powers = None                             #APRILE

OPT_LEVEL = 2
LIMIT_Y, LIMIT_Z = 8,2

def setOptimizationLevel(level):
  assert level in [0,1,2]
  global OPT_LEVEL
  OPT_LEVEL = level

def setCentralMoments(max_side):
  """
  The argument can be a single integer or a terne of integers,
  which are the max length of a block (rectangle) in the three
  Cartesian directions of the space.
  Compute and store all central moments of rectangles of size DX x DY x DZ
  with DX,DY,DZ in [1,max_side], DX>=DY>=DZ..
  - for DX==DY==DZ==1 the only non-zero moment is mu_000 
  - otherwise the only non-zero moments are mu_000 and mu_200, mu_020, mu_002
  """
  
  # manage argument
  if type(max_side) is int:
    MoreMax,MidMax,LessMax = max_side, max_side, max_side
  else:
    MoreMax,MidMax,LessMax = max_side # already sorted, decreasing
  
  CentrMom000 = dict()
  CentrMom200 = dict()
  CentrMom020 = dict()
  CentrMom002 = dict()

  #DX==DY==DZ==1
  CentrMom000[(1,1,1)] = 1
  CentrMom200[(1,1,1)] = 0
  CentrMom020[(1,1,1)] = 0
  CentrMom002[(1,1,1)] = 0
  
  # Compute auxiliary values: 
  # sum_half = 0.5^2 + 1.5^2 + 2.5^2 + ...
  # sum_full = 1^2   + 2^2   + 3^2   + ...
  sum_half=[0,0.25]
  sum_full=[0]
  for f in range(2,MoreMax+1):
    sum_half.append(sum_half[f-1]+((f-0.5)**2))
  for f in range(1,MoreMax+1):
    sum_full.append(sum_full[f-1]+(f**2))

  # DX>=1 and DY==DZ==1
  for edgeX in range(1,MoreMax+1):
     CentrMom000[(edgeX,1,1)] = edgeX
     CentrMom020[(edgeX,1,1)] = 0
     CentrMom002[(edgeX,1,1)] = 0
     if edgeX%2==0:
       CentrMom200[(edgeX,1,1)] = 2*sum_half[edgeX//2]
     else:
       CentrMom200[(edgeX,1,1)] = 2*sum_full[edgeX//2]
     #print("a) set ",(edgeX,1,1))

  # DZ==1 and DX>=DY>1
  for edgeY in range(2,MidMax+1):
    for edgeX in range(edgeY,MoreMax+1):
      CentrMom000[(edgeX,edgeY,1)] = edgeX*edgeY
      CentrMom002[(edgeX,edgeY,1)] = 0
      if edgeX%2==0: sum_x = sum_half
      else: sum_x = sum_full
      if edgeY%2==0: sum_y = sum_half
      else: sum_y = sum_full
      CentrMom200[(edgeX,edgeY,1)] = 2*edgeY*sum_x[edgeX//2]
      CentrMom020[(edgeX,edgeY,1)] = 2*edgeX*sum_y[edgeY//2]
      #print("b) set ",(edgeX,edgeY,1))

  # DX>=DY>=DZ>1
  for edgeZ in range(2,LessMax+1):
    for edgeY in range(edgeZ,MidMax+1):
      for edgeX in range(edgeY,MoreMax+1):
        CentrMom000[(edgeX,edgeY,edgeZ)] = edgeX*edgeY*edgeZ
        if edgeX%2==0: sum_x = sum_half
        else: sum_x = sum_full
        if edgeY%2==0: sum_y = sum_half
        else: sum_y = sum_full
        if edgeZ%2==0: sum_z = sum_half
        else: sum_z = sum_full
        CentrMom200[(edgeX,edgeY,edgeZ)] = 2*edgeY*edgeZ*sum_x[edgeX//2]
        CentrMom020[(edgeX,edgeY,edgeZ)] = 2*edgeX*edgeZ*sum_y[edgeY//2]
        CentrMom002[(edgeX,edgeY,edgeZ)] = 2*edgeX*edgeY*sum_z[edgeZ//2]
        #print("c) set ",(edgeX,edgeY,edgeZ))

  return (CentrMom000, CentrMom200, CentrMom020, CentrMom002)

def preprocessing(ibr):
  #assert isinstance(ibr,BW_BlockImage3D)

  # precompute central moments for new method
  global CC000, CC200, CC020, CC002
  # precompute matrix for traditional method
  global powers
  
  # manage optimization level
  if OPT_LEVEL==0:
    maximum = max(ibr.max_triplet())
    CC000, CC200, CC020, CC002 = setCentralMoments(maximum)
  elif OPT_LEVEL==1:
    maximum = sorted(ibr.max_triplet(),reverse=True)
    CC000, CC200, CC020, CC002 = setCentralMoments(maximum)
  elif OPT_LEVEL==2:
    maximum = sorted(ibr.max_triplet(),reverse=True)
    maximum[1] = min(maximum[1],LIMIT_Y)
    maximum[2] = min(maximum[2],LIMIT_Z)
    CC000, CC200, CC020, CC002 = setCentralMoments(maximum)  
    powers = PowerMatrix( 3, ibr.origsize )

def preprocessing_once(max_side):
  # precompute matrix for traditional method
  global powers
  powers = PowerMatrix( 3, max_side )
  
  # precompute central moments for new method
  global CC000, CC200, CC020, CC002
  if OPT_LEVEL==0:
    CC000, CC200, CC020, CC002 = setCentralMoments(max_side)
  else: # 1,2
    CC000, CC200, CC020, CC002 = setCentralMoments((max_side,LIMIT_Y,LIMIT_Z))
  
def blockMoments(ibr):
  """
  Compute all moments m_{p,q,r} for p,q,r>=0 and p+q+r<=3
  of a 3D image given as a set of blocks
  """
  #assert isinstance(ibr,BW_BlockImage3D)

  # orders is the global variable imported from commons3D
  # initialize moments to be computed
  MM = {key:0 for key in orders}
  
  #print()
  for p,q,r in orders:
     MM[(p,q,r)] = 0 
  
  global NUOVO
  global VECCHIO 
  NUOVO,VECCHIO = 0,0 #APRILE
  
  for b in ibr.block: # cycle on blocks

     if OPT_LEVEL>1:
       dimens = (b.x1-b.x0+1, b.y1-b.y0+1, b.z1-b.z0+1) #APRILE
       ordered = sorted(dimens)
       if (ordered[1]>LIMIT_Y) or (ordered[0]>LIMIT_Z): #APRILE faccio al modo vecchio
         #print("VECCHIO MODO",dimens,ordered)
         for p,q,r in orders:
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
           if r==0: mz = dimens[2]
           else:
             mz = powers.valueSum(r, b.z1)
             if b.z0>0:
                mz -= powers.valueSum(r, b.z0-1)
           if (mx or my or mz): MM[(p,q,r)] += (mx*my*mz)
         VECCHIO += 1
         continue
     #FINE APRILE   
     
     #print("NUOVO MODO",dimens,ordered)
     NUOVO += 1
     # barycenter
     xx = 0.5*(b.x1+b.x0)
     yy = 0.5*(b.y1+b.y0)
     zz = 0.5*(b.z1+b.z0)
     #print("Momenti di ",b, " di ",b.pixel_num(), " pixel, baricentro ",(xx,yy,zz))

     # retrieve central moments of block
     if (b.x1-b.x0)>=(b.y1-b.y0) and (b.y1-b.y0)>=(b.z1-b.z0):
               #print("  key xyz")
               key = (b.x1-b.x0+1, b.y1-b.y0+1, b.z1-b.z0+1)
               central000 = CC000[key]
               central200 = CC200[key]
               central020 = CC020[key]
               central002 = CC002[key]
     elif (b.x1-b.x0)>=(b.z1-b.z0) and (b.z1-b.z0)>=(b.y1-b.y0):
               #print("  key xzy,  020:=002 e 022:=020 ")
               key = (b.x1-b.x0+1, b.z1-b.z0+1, b.y1-b.y0+1)
               central000 = CC000[key]
               central200 = CC200[key]
               central020 = CC002[key]
               central002 = CC020[key]
     elif (b.y1-b.y0)>=(b.x1-b.x0) and (b.x1-b.x0)>=(b.z1-b.z0):
               #print("  key yxz,  200:=020 e 020:=200 ")
               key = (b.y1-b.y0+1, b.x1-b.x0+1, b.z1-b.z0+1)
               central000 = CC000[key]
               central200 = CC020[key]
               central020 = CC200[key]
               central002 = CC002[key]
     elif (b.y1-b.y0)>=(b.z1-b.z0) and (b.z1-b.z0)>=(b.x1-b.x0):
               #print("  key yzx,  200:=002 e 020:=200 e 002:=020")
               key = (b.y1-b.y0+1, b.z1-b.z0+1, b.x1-b.x0+1)
               central000 = CC000[key]
               central200 = CC002[key]
               central020 = CC200[key]
               central002 = CC020[key]
     elif (b.z1-b.z0)>=(b.x1-b.x0) and (b.x1-b.x0)>=(b.y1-b.y0):
               #print("  key zxy,  200:=020 e 020:=002 e 002:=200")
               key = (b.z1-b.z0+1, b.x1-b.x0+1, b.y1-b.y0+1)
               central000 = CC000[key]
               central200 = CC020[key]
               central020 = CC002[key]
               central002 = CC200[key]
     else: # (b.z1-b.z0)>=(b.y1-b.y0) and (b.y1-b.y0)>=(b.x1-b.x0)
               #print("  key zyx,  200:=002 e 002:=200")
               key = (b.z1-b.z0+1, b.y1-b.y0+1, b.x1-b.x0+1)
               central000 = CC000[key]
               central200 = CC002[key]
               central020 = CC020[key]
               central002 = CC200[key]
     ##print("Blocco",b,"chiave",key)
     # compute moments from central ones
     m000 = central000
     m100 = xx*central000
     m010 = yy*central000
     m001 = zz*central000 
     #print("\t m000=",m000,"\t m100=",m100,"\t m010=",m010,"\t m001=",m001)
     
     m011 = zz*m010
     m101 = xx*m001
     m110 = yy*m100
     #print("\t m011=",m011,"\t m0101=",m0101,"\t m110=",m110)
          
     m200 = central200 + xx*m100
     m020 = central020 + yy*m010
     m002 = central002 + zz*m001
     #print("\t m200=",m200,"\t m020=",m020,"\t m002=",m002)

     m021 = zz*m020
     m210 = yy*m200
     m120 = xx*m020
     m102 = xx*m002
     m012 = yy*m002
     m201 = zz*m200

     m300 = 3*xx*m200 -2*xx*xx*m100
     m030 = 3*yy*m020 -2*yy*yy*m010
     m003 = 3*zz*m002 -2*zz*zz*m001

     m111 = xx*m011

     MM[(0,0,0)] += m000
     MM[(1,0,0)] += int(m100)
     MM[(0,1,0)] += int(m010)
     MM[(1,1,0)] += int(m110)
     MM[(2,0,0)] += int(m200)
     MM[(0,2,0)] += int(m020)
     MM[(2,1,0)] += int(m210)
     MM[(1,2,0)] += int(m120)
     MM[(3,0,0)] += int(m300)
     MM[(0,3,0)] += int(m030)
     MM[(0,0,3)] += int(m003)

     MM[(0,0,1)] += int(m001)
     MM[(1,0,1)] += int(m101)
     MM[(0,1,1)] += int(m011)
     MM[(1,1,1)] += int(m111)
     MM[(2,0,1)] += int(m201)
     MM[(0,2,1)] += int(m021)
    
     MM[(0,0,2)] += int(m002)
     MM[(1,0,2)] += int(m102)
     MM[(0,1,2)] += int(m012)

  return MM
   
# da chiamare subito dopo blockMoments
def stampaGestione():
   print("N. Blocks processed with new and with traditional way",NUOVO,VECCHIO)
   #print("Blocchi gestiti col nuovo e col vecchio",NUOVO,VECCHIO)

#---------------------MAIN-----------------------

from spiliotis3D import extractBlocks, checkBlocks
from commons3D import main
import sys

if __name__ == "__main__":
   main(sys.argv, extractBlocks, preprocessing, blockMoments, "====3D Blocks, new method.")
   print("Blocchi gestiti col nuovo e col vecchio",NUOVO,VECCHIO)
   