import numpy as np
import matplotlib.pyplot as plt
import utils
import math
import cmath
import argparse

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

def HmatFun( V1, V2, N, ang1, ang2 ): 
	Hmat = utils.cmplxZeros( (N, N) )
	kLst = np.array( list( range(N) ), dtype=np.complex_ )
	HmatDiag = -V1 * np.cos( ( 2*math.pi*kLst + ang1 ) / N )	
	return np.diag( HmatDiag )
	
def HmatFunOld( V1, V2, N, ang1, ang2 ): 
	Hmat = utils.cmplxZeros( (N, N) )
	for k in range(N):
		Hmat[k,k] = -V1 * math.cos( ( 2*math.pi*k + ang1 ) / N )
	return Hmat

parser = argparse.ArgumentParser()
parser.add_argument('-o', type=str, nargs=1)

args = parser.parse_args()

filenamesLst = args.o

k_num = N
ang2pi = 2 * math.pi
ang_ln = ang_divide + 1
# ang_step = ang2pi / ang_divide
ang_lst = np.linspace( 0, ang2pi, ang_ln )
# ang_lst = ang_step * np.arange( ang_ln )
# ang_lst = utils.cmplxZeros(ang_ln);
# for ind in range(ang_ln):
	# ang_lst[ind] = ind * ang_step
	
for itV1 in range( vLst.shape[1] ):
	V1 = vLst[0, itV1]
	for itV2 in range( vLst.shape[1] ):
		V2 = vLst[1, itV2]
		for itang1 in range( ang_ln ):
			ang1 = ang_lst[itang1]
			for itang2 in range( ang_ln ):
				ang2 = ang_lst[itang2]
				Hmat = HmatFun( V1, V2, N, ang1, ang2 )
				vecLstTmp = vecLst[ itV1, itV2, itang1, itang2]
				HmatLst[ itV1, itV2, itang1, itang2] = np.dot( np.conj( vecLstTmp ), np.dot( Hmat, np.transpose( vecLstTmp ) ) )

evoNum = 10

EevoLst = utils.cmplxZeros( Elst.shape[0:4] + (evoNum,) )

T = np.absolute( Elst[:,:,:,:,0] )
cmLst = np.exp( - Elst / T[...,np.newaxis] )
cmLst = cmLst / np.sum( cmLst, axis = -1 )[...,np.newaxis]
cmLstCol = np.expand_dims( cmLst, axis = -1 )
cmLstRow = np.expand_dims( cmLst, axis = -2 )
ElstGrd1 = np.repeat( Elst[:,:,:,:,:,np.newaxis], Elst.shape[-1], axis=-1 )
ElstGrd2 = np.swapaxes( ElstGrd1, -1, -2 )

for itEvo in range( evoNum ):
	expE = np.exp( 1j * (ElstGrd1-ElstGrd2) * itEvo )
	EevoLst[ :,:,:,:, itEvo ] = np.squeeze( cmLstRow @ (HmatLst*expE) @ cmLstCol )

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