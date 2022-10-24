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
ang_step = ang2pi / ang_divide
ang_lst = ang_step * np.array( list( range(ang_ln) ) )
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
				for itn in range( N ):
					vecn = vecLst[ itV1, itV2, itang1, itang2, itn ]
					for itm in range( N ):
						vecm = vecLst[ itV1, itV2, itang1, itang2, itm ]
						HmatLst[ itV1, itV2, itang1, itang2, itn, itm ] = np.vdot( vecn, np.dot( Hmat, vecm ) )

cmLst = utils.cmplxZeros( (N,) )



# cmLst[0] = 1

evoNum = 10

EevoLst = utils.cmplxZeros( Elst.shape[0:4] + (evoNum,) )

for itV1 in range( vLst.shape[1] ):
	V1 = vLst[0, itV1]
	for itV2 in range( vLst.shape[1] ):
		V2 = vLst[1, itV2]
		for itang1 in range( ang_ln ):
			ang1 = ang_lst[itang1]
			for itang2 in range( ang_ln ):
				ang2 = ang_lst[itang2]
				ElstTmp = Elst[itV1, itV2, itang1, itang2]
				HmatTmp = HmatLst[itV1, itV2, itang1, itang2]
				T = ElstTmp[0]
				cmLst = np.exp( - ElstTmp / T )
				cmLst = cmLst / np.sum( cmLst )
				for itEvo in range( evoNum ):
					EevoTmp = 0 + 0j
					for itn in range(N):
						for itm in range(N):
							EevoTmp += cmLst[itn].conjugate() * cmLst[itm] * cmath.exp( 1j*(ElstTmp[itn] - ElstTmp[itm]) *itEvo ) * HmatTmp[itn,itm]
					EevoLst[ itV1, itV2, itang1, itang2, itEvo ] = EevoTmp


if filenamesLst != None:
	np.save( filenamesLst[0], EevoLst )