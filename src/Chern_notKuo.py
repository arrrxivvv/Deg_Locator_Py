# Chern number calculation functions #
import utils
import math
import cmath
import copy
import numpy as np
from numpy import linalg as LA
	
def PairLst(list1, list2): 
      
    paired_list = tuple(zip(list1, list2))  
    return paired_list 

def SortPairLst( pairLst ):
	def getKey(item):
		return item[0]
	return sorted(pairLst,key=getKey)
	
def RollE_Vec( eLst, vecLst, bndShft ):
	bndShftDir = utils.sign(bndShft)
	movedPos = int( (bndShftDir-1)/2 )
	for rep in range( 0, bndShft, utils.sign(bndShft) ):
		eLst = np.roll( eLst, bndShftDir )
		eLst[movedPos] = eLst[movedPos] -bndShftDir * 2 * math.pi
		vecLst = np.roll( vecLst, bndShftDir, axis=0 )
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
	
	bndShftMem = 0
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

			(quasiELstTmp, eVecLstTmp) = RollE_Vec( quasiELstTmp, eVecLstTmp, bndShftMem )
			
			if isFlq and isBandMatch and N > 2 and not(ind1==0 and ind2==0):
				indsThis = [ind1,ind2]
				if  (ind2 != 0):
					prevDir = 1
				else:
					prevDir = 0
				indsPrev = copy.deepcopy( indsThis )
				indsPrev[prevDir] -= 1
				indsThis = tuple(indsThis)
				indsPrev = tuple(indsPrev)
				
				for movedPos in range(0, -1-1, -1):
					bndShft = 2*movedPos + 1
					E0_prev = quasiELst[indsPrev][movedPos]
					E0_this = quasiELstTmp[movedPos]
					E1_prev = quasiELst[indsPrev][movedPos+bndShft]
					if bndShft*( E0_this - E0_prev ) > bndShft * 0.5 * ( E1_prev - E0_prev ): 
						(quasiELstTmp, eVecLstTmp) = RollE_Vec( quasiELstTmp, eVecLstTmp, bndShft )
						bndShftMem += bndShft
						break
			
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

