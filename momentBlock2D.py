"""
Computation of moments of a 2D image on the block 
decomposition by Spiliotis and Mertzios.
"""

#---------------------MOMENTS-----------------------

from bigmatrix import PowerMatrix
from commons2D import orders

#Matrix storing precomputed sums of powers
powers = None

def preprocessing(ibr):
  #assert isinstance(ibr,BW_BlockImage2D)
  global powers
  powers = PowerMatrix( 3, ibr.origsize )

def blockMoments(ibr):
  """
  Compute all moments m_{p,q} for p,q>=0 and p+q<=3
  of a 2D image given as a set of blocks
  """
  # orders is the global variable imported from commons2D
  MM = {key:0 for key in orders}

  for p,q in orders:
        #print("calcolo momenti ordine ",p,q)
        # value of moment
        MM[(p,q)] = 0 
        #print("  num blocchi",len(ibr.block))
        for b in ibr.block: # cycle on blocks
            #print("Momento ord ",(p,q), " di ",b, " di ",b.pixel_num(), " pixel")
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
            #print("Momento ord ",(p,q), " Blocco ",b," contrib= ",(mx,my), mx*my)
                
            if (mx or my):
                 MM[(p,q)] += (mx*my)
  return MM
   
#---------------------MAIN-----------------------

from spiliotis2D import extractBlocks
from commons2D import main
import sys

if __name__ == "__main__":
   main(sys.argv, extractBlocks, preprocessing, blockMoments, "====2D Blocks, old method.")
