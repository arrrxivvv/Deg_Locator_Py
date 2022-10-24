import Umat
import Chern
import utils
import functools

import time
import argparse

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-o', type=str, nargs=2)

args = parser.parse_args()

filenamesLst = args.o

# u=1.5
u_step = 0.01
u_min = 1.98
u_max = 2.02
u_lst = np.arange(u_min,u_max,u_step)
# U_mat_QWZ = functools.partial( Umat.U_mat_QWZ_raw, u=u )
N=2

N = 2
V1 = 5
V2 = 5
tau = 1
ang_divide = 10

U_mat_flq = functools.partial( Umat.U_mat_flq_raw, N=N, V1=V1, V2=V2, tau=tau )

# angdivide = np.arange(10,100,10)

chernLst = utils.cmplxZeros( ( u_lst.size,N ) )

t_start = time.time()

for ind in range(0,u_lst.size):
	u = u_lst[ind]
	U_mat_QWZ = functools.partial( Umat.U_mat_QWZ_raw, u=u )
	chernLst[ind] = Chern.Chern( N, U_mat_QWZ, ang_divide, isFlq = False )
	print( u_lst[ind] )
	utils.printArrPres( chernLst[ind] )



t_end = time.time()

print( "runtime: ", t_end-t_start )

utils.printArrPres( chernLst )

if filenamesLst != None:
	np.savetxt( filenamesLst[0], [u_lst], delimiter=',',fmt='%f' )
	np.savetxt( filenamesLst[1], np.real( np.round( chernLst, 2 ) ), delimiter=',',fmt='%f' )