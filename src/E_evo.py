import numpy as np
import matplotlib.pyplot as plt
import utils
import math
import cmath
import argparse
import functools
import time

N = 5
ang_divide = 20
chernFile = "Clst_n5_V_10_test.npy"
vFile = "V_n5_V_10_test.npy"
Efile = "Elst_n5_V_10_test.npy"
vecFile = "vec_n5_V_10_test.npy"

chernLst = np.load( chernFile )
vLst = np.load( vFile )
Elst = np.load( Efile )
vecLst = np.load( vecFile )

HmatLst = utils.cmplxZeros( vecLst.shape )

def HmatFun( ang1, ang2, V1, V2, N ): 
	Hmat = utils.cmplxZeros( (N, N) )
	kLst = np.array( list( range(N) ), dtype=np.complex_ )
	HmatDiag = -V1 * np.cos( ( 2*math.pi*kLst + ang1 ) / N )	
	return np.diag( HmatDiag )

def HmatVec( V1lst, V2lst, N, ang_divide ): 
	Hmat = 0
	k1lst = np.arange(N)
	k2lst = np.arange(N)
	ang_ln = ang_divide + 1
	ang_lst = np.linspace( 0, ang2pi, ang_ln )
	V1, V2, ang1, ang2, k1, k2 = np.meshgrid( V1lst, V2lst, ang_lst, ang_lst, k1lst, k2lst )
	U_ans = HmatElem( V1, V2, N, ang1, ang2, k1, k2 )
	return U_ans

def HmatElem( V1, V2, N, ang1, ang2, k1, k2 ): 
	Hmat = -V1 * np.cos( ( 2*math.pi*k1 + ang1 ) / N ) * (k1==k2)
	return Hmat
	
def HmatFunOld( V1, V2, N, ang1, ang2 ): 
	Hmat = utils.cmplxZeros( (N, N) )
	for k in range(N):
		Hmat[k,k] = -V1 * math.cos( ( 2*math.pi*k + ang1 ) / N )
	return Hmat

parser = argparse.ArgumentParser()
parser.add_argument('-o', type=str, nargs=2)

args = parser.parse_args()

filenamesLst = args.o

k_num = N
ang2pi = 2 * math.pi
ang_ln = ang_divide + 1
# ang_step = ang2pi / ang_divide
ang_lst = np.linspace( 0, ang2pi, ang_ln )
ang1grd, ang2grd = np.meshgrid( ang_lst, ang_lst )
# ang_lst = ang_step * np.arange( ang_ln )
# ang_lst = utils.cmplxZeros(ang_ln);
# for ind in range(ang_ln):
	# ang_lst[ind] = ind * ang_step
	
tic = time.time()

for itV1 in range( vLst.shape[1] ):
	V1 = vLst[0, itV1]
	for itV2 in range( vLst.shape[1] ):
		V2 = vLst[1, itV2]
		# HmatFunV = functools.partial( HmatFun, V1=V1, V2=V2, N=N )
		HmatLstTmp = HmatVec( [V1], [V2], N, ang_divide )
		vecLstTmp = vecLst[itV1, itV2]
		HmatLst[itV1, itV2] = np.conj( vecLstTmp ) @ HmatLstTmp @ np.swapaxes( vecLstTmp, -1, -2 )

toc = time.time()
print( toc-tic )

evoNum = 50

EkinLst = np.diagonal( HmatLst, axis1=-2, axis2=-1 )

EevoLst = utils.cmplxZeros( Elst.shape[0:4] + (evoNum,) )

T = np.absolute( Ekinlst[:,:,:,:,0] )
cmLst = np.exp( - EkinLst / T[...,np.newaxis] )
cmLst = cmLst / np.sum( cmLst, axis = -1 )[...,np.newaxis]
cmLstCol = np.expand_dims( cmLst, axis = -1 )
cmLstRow = np.expand_dims( cmLst, axis = -2 )
ElstGrd1 = np.repeat( Elst[:,:,:,:,:,np.newaxis], Elst.shape[-1], axis=-1 )
ElstGrd2 = np.swapaxes( ElstGrd1, -1, -2 )

tic = time.time()

for itEvo in range( evoNum ):
	expE = np.exp( 1j * (ElstGrd1-ElstGrd2) * itEvo )
	EevoLst[ :,:,:,:, itEvo ] = np.squeeze( cmLstRow @ (HmatLst*expE) @ cmLstCol )
	
toc = time.time()

print( toc - tic )

cmLst = utils.cmplxZeros( Elst.shape )
cmLst[:,:,:,:,0] = 1
cmLstCol = np.expand_dims( cmLst, axis = -1 )
cmLstCol = np.swapaxes( np.conj( vecLst ), -1, -2 ) @ cmLstCol
cmLstRow = np.swapaxes( np.conj( cmLstCol ), -1, -2 )
EevoLst2 = utils.cmplxZeros( Elst.shape[0:4] + (evoNum,) )

tic = time.time()

for itEvo in range( evoNum ):
	expE = np.exp( 1j * (ElstGrd1-ElstGrd2) * itEvo )
	EevoLst2[ :,:,:,:, itEvo ] = np.squeeze( cmLstRow @ (HmatLst*expE) @ cmLstCol )
	
toc = time.time()

print( toc - tic )

# for itV1 in range( vLst.shape[1] ):
	# V1 = vLst[0, itV1]
	# for itV2 in range( vLst.shape[1] ):
		# V2 = vLst[1, itV2]
		# for itang1 in range( ang_ln ):
			# ang1 = ang_lst[itang1]
			# for itang2 in range( ang_ln ):
				# ang2 = ang_lst[itang2]
				# ElstTmp = Elst[itV1, itV2, itang1, itang2]
				# HmatTmp = HmatLst[itV1, itV2, itang1, itang2]
				# T = ElstTmp[0]
				# for itEvo in range( evoNum ):
					# EevoTmp = 0 + 0j
					# En1, Em1 = np.meshgrid( ElstTmp, ElstTmp )
					# expE = np.exp( 1j * (En1-Em1) * itEvo )
					# EevoTmp = np.vdot( En1, np.dot( expE*HmatTmp, Em1 ) )
					# EevoLst[ itV1, itV2, itang1, itang2, itEvo ] = EevoTmp


if filenamesLst != None:
	np.save( filenamesLst[0], EevoLst )
	np.save( filenamesLst[1], EevoLst2 )