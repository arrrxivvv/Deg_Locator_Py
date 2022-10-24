# Chern number calculation functions #
import utils
import math
import cmath
import numpy as np
from numpy import linalg as LA
	
def PairLst(list1, list2): 
      
    paired_list = tuple(zip(list1, list2))  
    return paired_list 

def SortPairLst( pairLst ):
	def getKey(item):
		return item[0]
	return sorted(pairLst,key=getKey)
	
def RollE_Vec( eLst, vecLst, bndShft )
	bndShftDir = utils.sign(bndShft)
	movedPos = (bndShftDir-1)/2
		for rep in range( 0, bndShft, utils.sign(bndShft) )
			eLst = np.roll( quasiELstTmp, bndShftDir )
			eLst[movedPos] += -bndshft * 2 * math.pi
			eVecLstTmp = np.roll( eVecLstTmp, bndShftDir, axis=0 )
	return (eLst, vecLst)
	
def SortABbyA( key, object ):
	pairedLst = PairLst( key, object )
	sortedLst = SortPairLst( pairedLst )
	return list(zip(*sortedLst)) 

def Chern( N, U_mat_fun, ang_divide, isFlq = True, isBandMatch = True ):
	k_num = N
	ang2pi = 2 * math.pi
	ang_ln = ang_divide + 1
	ang_step = ang2pi / ang_divide
	ang_lst = utils.cmplxZeros(ang_ln);
	for ind in range(ang_ln):
		ang_lst[ind] = ind * ang_step

	#U_mat_lst = utils.cmplxZeros( shape = (ang_ln, ang_ln, k_num, k_num) )
	U_mat_temp = utils.cmplxZeros( shape = (k_num,k_num) )
	quasiELst = utils.cmplxZeros( (ang_ln, ang_ln, N) )
	#quasiELstTheo = utils.cmplxZeros( (ang_ln, ang_ln, N) )
	eVecLst = utils.cmplxZeros( (ang_ln, ang_ln, N, N) )
	
	bndShft = 0
	for ind1 in range(ang_ln):
		ang1 = ang_lst[ind1]
		for ind2 in range(ang_ln):
			ang2 = ang_lst[ind2]
			U_mat_temp = U_mat_fun( ang1, ang2 )
			
			eValLstTmp,eVecLstTmp = LA.eig(U_mat_temp)
			eVecLstTmp = np.transpose(eVecLstTmp)
			if isFlq:
				quasiELstTmp = np.log(eValLstTmp) / (-1j)
			else:
				quasiELstTmp = eValLstTmp
			quasiELstImTmp = np.imag(quasiELstTmp)
			quasiELstTmp = np.real(quasiELstTmp)
			#quasiELstTmp = eValLstTmp
			
			# (quasiELst[ind1][ind2], eVecLst[ind1][ind2]) = SortABbyA( quasiELstTmp, eVecLstTmp )			
			( quasiELstTmp, eVecLstTmp ) = SortABbyA( quasiELstTmp, eVecLstTmp )

			# (quasiELstTmp, eVecLstTmp) = RollE_Vec( quasiELstTmp, eVecLstTmp )
			
			if isFlq and isBandMatch and ind2 != 0 and N > 2:
				E0_prev = quasiELst[ind1][ind2-1][0]
				E0_this = quasiELstTmp[0]
				E1_prev = quasiELst[ind1][ind2-1][1]
				Emax_prev = quasiELst[ind1][ind2-1][N-1]
				Emax_this = quasiELstTmp[N-1]
				Emax1_prev = quasiELst[ind1][ind2-1][N-2]
				if ( E0_this - E0_prev ) > 1/2 * ( E1_prev - E0_prev ): 
					quasiELstTmp = np.roll( quasiELstTmp, 1 )
					quasiELstTmp[0] -= 2 * math.pi
					eVecLstTmp = np.roll( eVecLstTmp, 1, 1 )
				elif (Emax_prev - Emax_this) > 1/2 * (Emax_prev - Emax1_prev):
					quasiELstTmp = np.roll( quasiELstTmp, -1 )
					quasiELstTmp[N-1] += 2 * math.pi
					eVecLstTmp = np.roll( eVecLstTmp, -1, 0 )
			
			quasiELst[ind1][ind2] = quasiELstTmp
			eVecLst[ind1][ind2] = eVecLstTmp

	link1Lst = utils.cmplxZeros( (ang_ln, ang_ln-1) )
	link2Lst = utils.cmplxZeros( (ang_ln-1, ang_ln) )	
			
	chernLst = utils.cmplxZeros(N)
	uProdLst = utils.cmplxZeros( (ang_ln,ang_ln,N) )
	f12Lst = utils.cmplxZeros( (ang_ln,ang_ln,N) )
	for indN in range(N):
		f12Sum = 0
		for ind1 in range(ang_ln):
			for ind2 in range(ang_ln-1):
				dotProd = np.vdot( eVecLst[ind1][ind2][indN], eVecLst[ind1][ind2+1][indN] )
				link1Lst[ind1][ind2] = dotProd / abs( dotProd )	
		for ind1 in range(ang_ln-1):
			for ind2 in range(ang_ln):
				dotProd = np.vdot( eVecLst[ind1][ind2][indN], eVecLst[ind1+1][ind2][indN] )
				link2Lst[ind1][ind2] = dotProd / abs( dotProd )	
		for ind1 in range(ang_ln-1):
			for ind2 in range(ang_ln-1):
				f12 = cmath.log( link1Lst[ind1][ind2] / link1Lst[ind1+1][ind2] * link2Lst[ind1][ind2+1] / link2Lst[ind1][ind2] )
				f12Sum = f12Sum + f12
				
				uProdLst[ind1][ind2][indN] = link1Lst[ind1][ind2] / link1Lst[ind1+1][ind2] * link2Lst[ind1][ind2+1] / link2Lst[ind1][ind2]
				f12Lst[ind1][ind2][indN] = f12
		chernLst[indN] = f12Sum / (2 * math.pi * 1j )

	if isFlq:
		chernLst += 1 / N
	return chernLst

