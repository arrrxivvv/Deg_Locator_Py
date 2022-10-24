import Umat
import Chern
import utils
import functools

import time

import numpy as np

u=1.5
U_mat_QWZ = functools.partial( Umat.U_mat_QWZ_raw, u=u )
N=2

N = 3
V1 = 5
V2 = 5
tau = 1
#ang_divide = 100

U_mat_flq = functools.partial( Umat.U_mat_flq_raw, N=N, V1=V1, V2=V2, tau=tau )

angdivide = np.arange(10,100,10)

chernLst = utils.cmplxZeros( ( angdivide.size,N ) )

t_start = time.time()

for ind in range(0,angdivide.size):
	chernLst[ind] = Chern.Chern( N, U_mat_flq, angdivide[ind], isFlq = True )
	print( angdivide[ind] )
	utils.printArrPres( chernLst[ind] )

t_end = time.time()

print( "runtime: ", t_end-t_start )

utils.printArrPres( chernLst )