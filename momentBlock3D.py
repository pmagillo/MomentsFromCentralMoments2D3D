"""
Computation of moments of a 3D image on the block 
decomposition by Spiliotis and Mertzios.
"""

#---------------------MOMENTS-----------------------

from bigmatrix import PowerMatrix
from commons3D import orders

#Matrix storing precomputed sums of powers
powers = None

def preprocessing(ibr):
  #assert isinstance(ibr,BW_BlockImage2D)
  global powers
  powers = PowerMatrix( 3, ibr.origsize )

def blockMoments(ibr):
  """
  Compute all moments m_{p,q,r} for p,q,r>=0 and p+q+r<=3
  of a 3D image given as a set of blocks
  """

  # orders is the global variable imported from commons3D
  MM = {key:0 for key in orders}
  
  for p,q,r in orders:
        #print("calcolo momenti ordine ",p,q,r)
        # value of moment
        MM[(p,q,r)] = 0 
        #print("  num blocchi",len(ibr.block))
        #print()
        for b in ibr.block: # cycle on blocks
            #print("Momento ord ",(p,q,r), " di ",b, " di ",b.pixel_num(), " pixel")
            if p==0: mx = b.x1-b.x0+1
            else:
              mx = powers.valueSum(p, b.x1)
              if b.x0>0:
                 mx -= powers.valueSum(p, b.x0-1)
            if q==0: my = b.y1-b.y0+1
            else:
              my = powers.valueSum(q, b.y1)
              if b.y0>0:
                 my -= powers.valueSum(q, b.y0-1)
            if r==0: mz = b.z1-b.z0+1
            else:
              mz = powers.valueSum(r, b.z1)
              if b.z0>0:
                 mz -= powers.valueSum(r, b.z0-1)
            #print("   Blocco ",b," ha contributi=",(mx,my,mz), mx*my*mz)
            """
            if (p,q,r)==(0,0,0) and mx*my*mz != b.pixel_num():
                print(" Su x: dovrei avere ",b.x1-b.x0+1," e ho ",mx)
                if b.x0>1: print("   = ",powers.valueSum(p, b.x1)," - ",powers.valueSum(p, b.x0-1))
                else: print("   = ",powers.valueSum(p, b.x1))
                print(" Su y: dovrei avere ",b.y1-b.y0+1," e ho ",my)
                if b.y0>1: print("   = ",powers.valueSum(q, b.y1)," - ",powers.valueSum(q, b.y0-1))
                else: print("   = ",powers.valueSum(q, b.y1))
                print(" Su z: dovrei avere ",b.z1-b.z0+1," e ho ",mz)
                if b.z0>1: print("   = ",powers.valueSum(r, b.z1)," - ",powers.valueSum(r, b.z0-1))
                else: print("   = ",powers.valueSum(r, b.z1))
            """    
            if (mx or my or mz):
                 MM[(p,q,r)] += (mx*my*mz)
  return MM
   
#---------------------MAIN-----------------------

from spiliotis3D import extractBlocks
from commons3D import main
import sys

if __name__ == "__main__":
   main(sys.argv, extractBlocks, preprocessing, blockMoments, "====3D Blocks, old method.")
