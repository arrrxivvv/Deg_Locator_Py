import numpy as np
from numpy import linalg as LA
import utils
import time

# num_it = 50**3
num_it = 50**3
N=10

time_start = time.time()

for it in range(0,num_it,1):
	mat = np.random.rand( N,N )
	LA.eig( mat )
	

time_end = time.time()

print( "time: ", time_end - time_start )