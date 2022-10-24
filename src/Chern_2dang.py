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
	
def ShftOrNot( ShftFunOpt, quasiELst, quasiELstTmp, ind1, ind2, movedPos, bndShft ):
	indsThis = [ind1,ind2]
	if  (ind2 != 0):
		prevDir = 1
	else:
		prevDir = 0
	indsPrev = copy.deepcopy( indsThis )
	indsPrev[prevDir] -= 1
	indsThis = tuple(indsThis)
	indsPrev = tuple(indsPrev)
	
	quasiELstPrev = quasiELst[indsPrev]

	return ShftFunOpt( quasiELstPrev, quasiELstTmp, movedPos, bndShft )
	
def ShftOrNot2pi( quasiELstPrev, quasiELstTmp, movedPos, bndShft ):
	sumThis = np.sum( quasiELstTmp )
	sumPrev = np.sum( quasiELstPrev )
	
	return ( abs( sumThis - sumPrev ) > ( sumThis - ( sumPrev + bndShft*2*math.pi ) ) )
		
	
	
def ShftOrNot0order( quasiELstPrev, quasiELstTmp, movedPos, bndShft ):
	# bndShft = 2*movedPos + 1
	oppositePos = -1 - movedPos
	E0_prev = quasiELstPrev[movedPos]
	E0_this = quasiELstTmp[movedPos]
	E0_shft = quasiELstTmp[oppositePos] - bndShft * 2*math.pi
	
	return ( abs( E0_this - E0_prev ) > abs( E0_shft - E0_prev ) )
	
def ShftOrNotSum( quasiELstPrev, quasiELstTmp, movedPos, bndShft ):
	# bndShft = 2*movedPos + 1
	quasiELstShft = RollE_Vec( quasiELstTmp, quasiELstTmp, bndShft )[0]
	
	diff = abs( np.sum( np.subtract( quasiELstTmp, quasiELstPrev ) ) )
	diff_shft = abs( np.sum( np.subtract( quasiELstShft, quasiELstPrev ) ) )
	
	return ( diff_shft < diff )

def ShftOrNotHalfGap( quasiELstPrev, quasiELstTmp, movedPos, bndShft ):
	E0_prev = quasiELstPrev[movedPos]
	E0_this = quasiELstTmp[movedPos]
	E1_prev = quasiELstPrev[movedPos+bndShft]
	
	return bndShft*( E0_this - E0_prev ) > bndShft * 0.5 * ( E1_prev - E0_prev )

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
				for movedPos in range(0, -1-1, -1):
					bndShft = 2*movedPos + 1
					if ShftOrNot( ShftOrNotSum, quasiELst, quasiELstTmp, ind1, ind2, movedPos, bndShft ):
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
	return (chernLst, quasiELst, eVecLst)
	
def ChernEigen( N, ang_divide, eValLst, eVecLst, isFlq = True, isBandMatch = True ):
	k_num = N
	ang2pi = 2 * math.pi
	ang_ln = ang_divide + 1
	ang_step = ang2pi / ang_divide
	ang_lst = np.linspace( 0, ang2pi, ang_ln )
	quasiELst = utils.cmplxZeros( (ang_ln, ang_ln, N) )
	
	bndShftMem = 0
	for ind1 in range(ang_ln):
		ang1 = ang_lst[ind1]
		for ind2 in range(ang_ln):
			ang2 = ang_lst[ind2]
			
			eValLstTmp = eValLst[ind1, ind2]
			eVecLstTmp = eVecLst[ind1, ind2]
			if isFlq:
				quasiELstTmp = np.log(eValLstTmp) / (-1j)
			else:
				quasiELstTmp = eValLstTmp
			quasiELstImTmp = np.imag(quasiELstTmp)
			quasiELstTmp = np.real(quasiELstTmp)			
			
			( quasiELstTmp, eVecLstTmp ) = SortABbyA( quasiELstTmp, eVecLstTmp )

			(quasiELstTmp, eVecLstTmp) = RollE_Vec( quasiELstTmp, eVecLstTmp, bndShftMem )
			
			if isFlq and isBandMatch and N > 2 and not(ind1==0 and ind2==0):				
				for movedPos in range(0, -1-1, -1):
					bndShft = 2*movedPos + 1
					if ShftOrNot( ShftOrNotSum, quasiELst, quasiELstTmp, ind1, ind2, movedPos, bndShft ):
						(quasiELstTmp, eVecLstTmp) = RollE_Vec( quasiELstTmp, eVecLstTmp, bndShft )
						bndShftMem += bndShft
						break
			
			quasiELst[ind1][ind2] = quasiELstTmp
			eVecLst[ind1,ind2] = np.array( eVecLstTmp )

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
	return (chernLst, quasiELst, eVecLst)
	
def EigenSort( eLst, vecLst, dim_param )
	shape_param = eLst.shape[0:dim_param+1]
	ln_param = np.product( shape_param )
	
	for itLn in range(0,ln_param):
		id = 

