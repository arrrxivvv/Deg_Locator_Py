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
	
def ShftOrNot( ShftFunOpt, quasiELst, quasiELstTmp, indTuple, movedPos, bndShft ):
	indsThis = np.asarray( indTuple )
	non0pos = np.nonzero( indsThis )
	if not np.any( indsThis!=0 ):
		return False
	else:
		prevDir = non0pos[-1]
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
			
			if isFlq and isBandMatch and N > 2:				
				for movedPos in range(0, -1-1, -1):
					bndShft = 2*movedPos + 1
					if ShftOrNot( ShftOrNotSum, quasiELst, quasiELstTmp, [ind1, ind2] , movedPos, bndShft ):
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
			
			if isFlq and isBandMatch and N > 2:				
				for movedPos in range(0, -1-1, -1):
					bndShft = 2*movedPos + 1
					if ShftOrNot( ShftOrNotSum, quasiELst, quasiELstTmp, [ind1, ind2], movedPos, bndShft ):
						(quasiELstTmp, eVecLstTmp) = RollE_Vec( quasiELstTmp, eVecLstTmp, bndShft )
						bndShftMem += bndShft
						break
			
			quasiELst[ind1,ind2] = quasiELstTmp
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
	
def EigenSort( eValLst, eVecLst, N, dim_param, isFlq = True ):
	shape_param = eValLst.shape[0:dim_param]
	ln_param = np.product( shape_param )
	
	bndShftMem = 0
	for itLn in range(0,ln_param):
		ind = np.unravel_index( itLn, shape_param )
		eValLstTmp = eValLst[ind]
		eVecLstTmp = eVecLst[ind]
		quasiELst = utils.cmplxZeros( eValLst.shape )
		if isFlq:
			quasiELstTmp = np.log(eValLstTmp) / (-1j)
		else:
			quasiELstTmp = eValLstTmp
		quasiELstImTmp = np.imag(quasiELstTmp)
		quasiELstTmp = np.real(quasiELstTmp)			
		
		( quasiELstTmp, eVecLstTmp ) = SortABbyA( quasiELstTmp, eVecLstTmp )

		(quasiELstTmp, eVecLstTmp) = RollE_Vec( quasiELstTmp, eVecLstTmp, bndShftMem )
		
		if isFlq:				
			for movedPos in range(0, -1-1, -1):
				bndShft = 2*movedPos + 1
				if ShftOrNot( ShftOrNotSum, quasiELst, quasiELstTmp, ind, movedPos, bndShft ):
					(quasiELstTmp, eVecLstTmp) = RollE_Vec( quasiELstTmp, eVecLstTmp, bndShft )
					bndShftMem += bndShft
					break
		
		quasiELst[ind] = quasiELstTmp
		eVecLst[ind] = np.array( eVecLstTmp )
	
	return eValLst, eVecLst

def BerryField( eVecLst, N, dim_param ):
	shape_param = eVecLst.shape[0:dim_param]
	ln_param = np.product( shape_param )
	dimLst = np.arange( len(shape_param) )
	
	linkLst = [0]*dim_param
	for it in range(0,dim_param):
		linkSize = shape_param
		# linkSize[it] -= 1
		linkLst[it] = utils.cmplxZeros( linkSize )
	
	numF = int( dim_param * (dim_param-1) / 2 )
	fLst = [0]*numF
	for dim in range(0, dim_param):
		linkLst[dim][:] = 0
		eVecLstSh = np.roll( eVecLst, -1, axis=dim )
		dotProd = np.multiply( np.conj(eVecLstSh), eVecLst ).sum(-1)
		linkLst[dim] = dotProd / np.abs(dotProd)
	itF = 0
	for dim1 in range(0, dim_param):
		for dim2 in range(dim1+1, dim_param):
			link1 = linkLst[dim1]
			link2 = linkLst[dim2]
			link1sh = np.roll( link1, -1, dim2 )
			link2sh = np.roll( link2, -1, dim1 )
			uProd =  link1 / link1sh * link2sh / link2
			fLstRaw = np.log( link1 / link1sh * link2sh / link2 ) / (2 * math.pi * 1j)
			slcInd = tuple( [slice(0,-1)]*dim_param )
			fLst[itF] = fLstRaw[slcInd]
			itF += 1
	return fLst