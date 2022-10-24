# main script
import argparse
import time

import utils
import Chern
import Umat

import math
import cmath
import numpy as np
from numpy import linalg as LA
import functools
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

def printArrPres( arr ):
	print( arr )
	print( np.round( arr, 2) )
	print( np.real( np.round( arr, 2) ) )

parser = argparse.ArgumentParser()
parser.add_argument('-o', type=str, nargs=4)

args = parser.parse_args()

filenamesLst = args.o

N = 5
V1 = 10
V2 = 10
tau = 1
ang_divide = 20

U_mat_flq = functools.partial( Umat.U_mat_flq_raw, N=N, V1=V1, V2=V2, tau=tau )

u=3.4
U_mat_QWZ = functools.partial( Umat.U_mat_QWZ_raw, u=u )


U_mat_fun = U_mat_flq

#U_mat

V_step = 1
V_min = 1
V_max = 10 + V_step

V1_lst = np.arange( V_min, V_max, V_step )
V2_lst = np.arange( V_min, V_max, V_step )

chernLst = utils.cmplxZeros( ( V1_lst.shape ) + (V2_lst.shape) + (N,) )
quasiELst = utils.cmplxZeros( ( V1_lst.shape ) + (V2_lst.shape) + (ang_divide+1,) + (ang_divide+1,) + (N,) )
eVecLst = utils.cmplxZeros( ( V1_lst.shape ) + (V2_lst.shape) + (ang_divide+1,) + (ang_divide+1,) + (N,) + (N,) )

tic = time.time()

U_mat_lst = Umat.U_mat_vec( Umat.U_flq_elem, N, tau, V1_lst, V2_lst, ang_divide )

toc = time.time()

print( toc-tic )

tic = time.time()

eValLst,eVecLst = LA.eig(U_mat_lst)

toc = time.time()

print( toc-tic )

ang2pi = 2 * math.pi
ang_ln = ang_divide + 1
ang_lst = np.linspace( 0, ang2pi, ang_ln )

U_mat_flq = functools.partial( Umat.U_mat_flq_raw, N=N, V1=V1, V2=V2, tau=tau )
U_mat_fun = U_mat_flq

tic = time.time()
for indV1 in range(0,V1_lst.size):
	V1 = V1_lst[indV1]
	for indV2 in range(0,V2_lst.size):
		# V = V_lst[indV]
		V2 = V2_lst[indV2]
		for ind1 in range(ang_ln):
			ang1 = ang_lst[ind1]
			for ind2 in range(ang_ln):
				ang2 = ang_lst[ind2]
				U_mat_temp = U_mat_fun( ang1, ang2 )

toc = time.time()

print( toc-tic )

tic = time.time()
for indV1 in range(0,V1_lst.size):
	V1 = V1_lst[indV1]
	for indV2 in range(0,V2_lst.size):
		# V = V_lst[indV]
		V2 = V2_lst[indV2]
		for ind1 in range(ang_ln):
			ang1 = ang_lst[ind1]
			for ind2 in range(ang_ln):
				ang2 = ang_lst[ind2]
				eValLstTmp,eVecLstTmp = LA.eig(U_mat_temp)

toc = time.time()

print( toc-tic )

U_mat_flq = functools.partial( Umat.U_mat_flq_raw_vec, N=N, V1=V1, V2=V2, tau=tau )
tic = time.time()
for indV1 in range(0,V1_lst.size):
	V1 = V1_lst[indV1]
	for indV2 in range(0,V2_lst.size):
		# V = V_lst[indV]
		V2 = V2_lst[indV2]
		for ind1 in range(ang_ln):
			ang1 = ang_lst[ind1]
			for ind2 in range(ang_ln):
				ang2 = ang_lst[ind2]
				U_mat_temp = U_mat_fun( ang1, ang2 )
				# eValLstTmp,eVecLstTmp = LA.eig(U_mat_temp)

toc = time.time()

for indV1 in range(0,V1_lst.size):
	V1 = V1_lst[indV1]
	for indV2 in range(0,V2_lst.size):
		# V = V_lst[indV]
		V2 = V2_lst[indV2]
		print(V1, V2)
		U_mat_flq = functools.partial( Umat.U_mat_flq_raw, N=N, V1=V1, V2=V2, tau=tau )
		(chernLstTmp, quasiELstTmp, eVecLstTmp) = Chern.Chern( N, U_mat_flq, ang_divide, isFlq = True, isBandMatch = True )
		chernLst[indV1,indV2] = chernLstTmp
		quasiELst[indV1,indV2] = quasiELstTmp
		eVecLst[indV1,indV2] = eVecLstTmp
		printArrPres( chernLst[indV1,indV2] )

# V1=1
# V2=1
# U_mat_flq = functools.partial( U_mat_flq_raw, N=N, V1=V1, V2=V2, tau=tau )
# print( Chern.Chern( N, U_mat_flq, ang_divide ) )

t_end = time.time()

print( "runtime: ", t_end - t_start )

# fig = plt.figure()
# plt.plot( V_lst, chernLst )
# fig.show()

if filenamesLst != None:
	np.save( filenamesLst[0], [V1_lst, V2_lst] )
	np.save( filenamesLst[1], np.real( np.round( chernLst, 2 ) ) )
	np.save( filenamesLst[2], quasiELst )
	np.save( filenamesLst[3], eVecLst )
	# np.savetxt( filenamesLst[0], [V1_lst, V2_lst], delimiter=',',fmt='%f' )
	# np.savetxt( filenamesLst[1], np.real( np.round( chernLst, 2 ) ), delimiter=',',fmt='%f' )