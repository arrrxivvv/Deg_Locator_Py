import math
import numpy as np
import functools

def cmplxZeros( shape ):
	return np.zeros(shape,dtype=np.complex_)
	
def printArrPres( arr ):
	print( arr )
	print( np.round( arr, 2) )
	print( np.real( np.round( arr, 2) ) )
	
def sign( num ):
	return int( math.copysign(1, num) )
# sign = functools.partial(math.copysign, 1) 
	