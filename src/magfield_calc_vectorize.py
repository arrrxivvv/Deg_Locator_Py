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
parser.add_argument('-o', type=str, nargs=2)

args = parser.parse_args()

filenamesLst = args.o

N = 2
tau = 1
ang_divide = 100

# U_mat_flq = functools.partial( Umat.U_mat_flq_raw, N=N, V1=V1, V2=V2, tau=tau )

u=3.4
U_mat_QWZ = functools.partial( Umat.U_mat_QWZ_raw, u=u )


# U_mat_fun = U_mat_flq

t_start = time.time()

#U_mat

V_step = 0.01
V1_min = 7.7
V1_max = 8.2 + V_step
V2_min = 4
V2_max = 5 + V_step

V1_lst = np.arange( V1_min, V1_max, V_step )
V2_lst = np.arange( V2_min, V2_max, V_step )
vLst = (V1_lst, V2_lst)

chernLst = utils.cmplxZeros( ( V1_lst.shape ) + (V2_lst.shape) + (N,) )
# quasiELst = utils.cmplxZeros( ( V1_lst.shape ) + (V2_lst.shape) + (ang_divide+1,) + (ang_divide+1,) + (N,) )
# eVecLst = utils.cmplxZeros( ( V1_lst.shape ) + (V2_lst.shape) + (ang_divide+1,) + (ang_divide+1,) + (N,) + (N,) )

# tic = time.time()

# U_mat_lst = Umat.U_mat_vec( Umat.U_flq_elem, N, tau, V1_lst, V2_lst, ang_divide )

# toc = time.time()

# print( toc-tic )

# tic = time.time()

# eValLst,eVecLst = LA.eig(U_mat_lst)

# eVecLst = np.swapaxes( eVecLst, -1, -2 )

# toc = time.time()

# print( toc-tic )

# Umat_test = Umat.U_mat_3sin( (math.pi/2,0), 0,0 )

ang2pi = 2 * math.pi
ang_ln = ang_divide + 1
ang_lst = np.linspace( 0, ang2pi, ang_ln )

uLst = np.arange(-1,1,0.2)

tic = time.time()

# U_mat_lst = Umat.U_mat_QWZ_vec( ang_lst,ang_lst,uLst )
U_mat_lst = Umat.U_mat_3sin( ang_lst,ang_lst,ang_lst )

toc = time.time()
print( toc-tic )

tic = time.time()

eValLst,eVecLst = LA.eig(U_mat_lst)

eVecLst = np.swapaxes( eVecLst, -1, -2 )
eValLst = np.squeeze( eValLst )
eVecLst = np.squeeze( eVecLst )

toc = time.time()

print( toc-tic )

tic = time.time()

ang_dim = 3
quasiELst, eVecLst = Chern.EigenSort( eValLst, eVecLst, N, ang_dim, isFlq=False )

toc = time.time()

print( toc-tic )

tic = time.time()

BfieldLst = Chern.BerryField( eVecLst, N, ang_dim )

cLst = np.sum( BfieldLst, axis=(1,2) )

toc = time.time()

print( toc-tic )


# toc = time.time()

# print( toc-tic )

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
	np.save( filenamesLst[0], BfieldLst )
	np.save( filenamesLst[1], cLst )