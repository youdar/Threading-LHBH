import numpy as nm
import string as st
import pylab as pl
from ThreadClasses import *
from Protein_properties_v3 import *
from AminoData import *
import copy as cp
import time
from ReadDate import *


def build_path_matrix_v3(row,col,Structure_rules):
	'''
	Build an empty mpath matrix
	'''
	tmp1 =[]
	mp = []
	base_length = len(Structure_rules)
	for i in range(row+1):
		for j in range(col):
			tmp2 = ScoreMatrix()
			jj = (j+1)%base_length - 1
			# 0:Not a corner, 1:can cut a corner, 3:after corner cut
			if j>1: 
				tmp2.CornerFlag = Structure_rules[jj].CornerFlag 	# 0:not cut, 1:corner cut allowed
			else:
				tmp2.CornerFlag = 0
			# 0:not a corner, 1:side start, 2:side end, -3:physical corner
			tmp2.CornerMarker = Structure_rules[jj].CornerMarker	
			tmp2.PointInOut = Structure_rules[jj].PointInOut
			tmp2.PosNum = Structure_rules[jj].PosNum
			tmp2.Special = Structure_rules[jj].Special	   
			tmp2.Turn = Structure_rules[jj].Turn
			tmp2.DiSulDis = Structure_rules[jj].DiSulDis
			tmp2.ResSize = Structure_rules[jj].ResSize			
			tmp2.ResSizeSigma = Structure_rules[jj].ResSizeSigma	
			tmp1 = tmp1+[tmp2]
		mp.extend([tmp1])
		tmp1=[]
	return mp
# End of build_path_matrix_v3

def StartCol(seq1,base_length):
	'''
	Looking for the starting position
	Note that the column in the matrix is starts with 0 
	seq1 is the sequence with the loop and corner cuting marks
	'''
	countFlag = True
	i = 0
	col = 0
	while seq1[i] != '|':
		if countFlag: col += 1
		if i == 50: 
			print 'Problem with StartCol'
			break
		if seq1[i] == '(': 
			countFlag = False
			col -= 1
		if seq1[i] == ')': countFlag = True
		i += 1
	col = base_length - col
	return col

def GetThreadScores(seq,seq1,Structure_rules,ScoresP,BaseData):
	'''
	Returns the score of a know thread seq1
	Scores : the position score (as a class variable)
	tot : the total score
	seq : is the sequence without the thread marks
	seq1 : the actual thread, with loop, corners and corner cuts
	'''
	# get the list of amino acid letters
	AAlist = AminoList()
	# list of symbols allowed in a thread
	ThreadMarks = ['|','-','(',')']
	# Dictionary of amino acids and their properties
	aminotable = AminoTable()
	# the size of our scoring matrix
	row = len(seq)			
	# basic repeated length in the structure. for LBH-18 6, for LBH-15 5
	base_length = len(Structure_rules)
	# calculate column number according to possible number of corners cut
	col = end_col(row,Structure_rules)
	# initializing a list of lists that will store the best routes in mpath matrix
	mpath = build_path_matrix_v3(row,col,Structure_rules)
	# Check at what column the thread starts
	icol = StartCol(seq1,base_length) 	
	irow = 0
	NoLoopFlag = True
	PScore = 0.0		# The side chain score
	i = 0			# Counts the steps in the seq1
	tot = 0.0

	# Score data collection
	# Collected for information on score distribution
	HphobicS = 0
	ChargeS = 0
	PolarS = 0
	SpecialS = 0
	SideChainVolS = 0
	CornerCutS = 0
	TurnVolS = 0
	LoopHydroS = 0
	SideCahinBondsS = 0
	LoopPenalty = 0
	LoopNotOnCornerS = 0
	TotalScoreS = 0


	for t in seq1:
		if t== '(': NoLoopFlag = False	# we've riched a loop
		if NoLoopFlag:
			if t not in ThreadMarks:
				AA = aminotable[t]
				Mpos = mpath[irow][icol]
				if t == 'C': turnLenght = mpath[irow][icol].DiSulDis
				else: turnLenght = mpath[irow][icol].Turn
				tempS = threadSidescore(t,seq1[0:i],turnLenght,ScoresP)	# Side chain bond score
				tot += tempS
				SideCahinBondsS += tempS
				#PScore += ThreadVolumeScore(seq1[0:i],BaseData)*ScoresP.TurnVolPenaly
				tempS = AA_VolumeScore(seq1[0:i+1],BaseData)*ScoresP.TurnVolPenaly
				tot += tempS
				TurnVolS += tempS
				tmpScr,total_scr = PosScore(AA,Mpos,ScoresP)			# position score
				HphobicS += tmpScr.hydrophobic_scr
				ChargeS += tmpScr.charge_scr
				PolarS += tmpScr.polar_scr
				SideChainVolS += tmpScr.size_scr
				SpecialS += tmpScr.special_scr
				tot += total_scr
				icol += 1 
				# if AA.name in ['H','K','R','D','E','Y','Q','S','N']:
				#print '~~~~~~~~~~~~~~~~~~~~~~~'
				#print AA.name,t
				#print seq1[0:i]
				#print 'side score',threadSidescore(t,seq1[0:i],turnLenght,ScoresP),PScore
				#print 'turn vol',AA_VolumeScore(seq1[0:i+1],BaseData)*ScoresP.TurnVolPenaly
				#print 'pos score',total_scr
				#print 'tot',tot

		if not NoLoopFlag:	# we are in a loop
			# Score for loops
			if t in AAlist: 
				# print '++++++++++++================='
				AA = aminotable[t]
				tot += AA.loop * ScoresP.SpecialRatio	# Adding special score for selected AA
				temp1 = AA.Hydrophobic * ScoresP.HydroRatio * ScoresP.LoopWeight	# adding the hydrophobicity score of AA on loop
				tot -= temp1
				temp2 = AA.polar * ScoresP.PolarRatio * ScoresP.LoopWeight	# adding the polarity score of AA on loop
				tot += temp2
				temp3 = AA.charge * ScoresP.ChargeRatio * ScoresP.LoopWeight	# adding the charge score of AA on loop
				tot += temp3
				LoopHydroS += (temp2 + temp3 - temp1)


		if t in AAlist: irow += 1        
		if t == '-': 					# Corner cut
			icol += 1					# Continue counting
			tempS = ScoresP.CornerCutPenalty		# panelty for corner cutting
			tot -= tempS
			CornerCutS -= tempS

		if t == ')': 
			NoLoopFlag = True		# End of loop
			temp1 = ScoresP.LoopPenalty / 3.0		# Small panelty for loops
			tot -= temp1
			LoopPenalty -= temp1
			temp2 = (mpath[irow][icol].CornerMarker in [0,2]) * ScoresP.LoopPenalty 	# apply loop penalty if not in corner
			tot -= temp2
			LoopNotOnCornerS -= temp2

		i += 1		

	# print 'tot un normalized',tot,float(row)
	tot += PScore		# the total score
	nrows = float(row)
	tot = tot / nrows	# Normelizing the to the number of amino acid

	# Print the scores
	print 'Total Score: ',tot
	print 'H-phobic Score: ',HphobicS/nrows
	print 'Charge Score: ',ChargeS/nrows
	print 'Polar Score: ',PolarS/nrows
	print 'Special Score: ',SpecialS/nrows
	print 'Side Chain Volume Score: ',SideChainVolS/nrows
	print 'Corner Cut Score: ',CornerCutS/nrows
	print 'Turn Volume Score: ',TurnVolS/nrows
	print 'Loop Hydro Score: ',LoopHydroS/nrows
	print 'Side Cahin Bonds Score: ',SideCahinBondsS/nrows
	print 'Loop Penalty: ',LoopPenalty/nrows
	print 'Loop Not On Corner Score: ',LoopNotOnCornerS/nrows

	# Return the totoal score
	return tot

def end_col(row,Structure_rules):
	# Getting the size of mpath (columns)
	# It looks to me that col = row + (base_length - 1) + int(row+1)/int(base_length-1)   should do the  
	# same, at least for the LHBH-18
	# This function over complecation is to make it possible to expand the program 
	# to other structure types
	add_col = 0
	base_length = len(Structure_rules)
	for i in range(row):
		j = (i+1)%base_length - 1	# the relative location (mod) in the base stracture
		tmp1 = (Structure_rules[j].CornerFlag > 1)	# 3:after corner cut
		tmp2 = (j>0)								# we are not immidiatly after a corner
		if tmp1 and tmp2: add_col += 1
		col = row + (base_length - 1) + add_col
	return col
# end end column calculation



def result_output_v4(seq,highscore,threadlist1,dScore):
	'''
	Printing results to file file_name_out.txt and to the screen
	'''

	# Collection of the best threads
	threadlist = [x.thread for x in threadlist1]	
	# Collection of the best scores, sorting and keeping only unique values
	threadlistScore = [x.score for x in threadlist1]	
	threadlistScore.sort()
	threadlistScore = uniqueList(threadlistScore)
	# Numebr of threads with scores up to dScore points difference from the highest score
	nThreads = len(threadlist)		

	# open file 'file_name_out.txt' for writing
	#f_name = file_name
	#f_path = file_path
	#filestr = f_path + 'out_' + f_name + '.txt'
	#f = open(filestr, 'w')	

	# setting strings for printing
	highscoreStr = 'Highest score (per amino acid):%.2f' % highscore
	bThreads = 'Numebr of threads with scores up to %g points difference from the highest score: %g' %(dScore,nThreads)

	# print to screen
	#print 'See %s.txt for results \n' % file_name
	print highscoreStr   
	print bThreads

	# Presentation of results:
	print 'Collection of those best threads'
	while threadlistScore != []:	
		for thrd in threadlist1:
			if thrd.score == threadlistScore[-1]:
				# Printing all threads		    
				cThreads = '%.1f :: %s' %(thrd.score, thrd.thread)
				print cThreads 
		del(threadlistScore[-1])
	print 'end'

# End resultoutput


def mpathCalc(Structure_rules,seq,BaseData,ScoresP):

	shortl = BaseData.shortl
	longl = BaseData.longl

	# the size of our scoring matrix
	row = len(seq)			
	# Add extra colomns to allow corner cutting
	base_length = len(Structure_rules)
	# Checking how many corners can be cut
	col = end_col(row,Structure_rules)
	# initializing a list of lists that will store the best routes in mpath matrix
	mpath = build_path_matrix_v3(row,col,Structure_rules)
	# mpath2 = build_path_matrix_v3(row,col,Structure_rules)
	# Scoring - creating mpath
	for irow in range(0,row):

		# Evaluating the possible range of colomns for each row
		# for each row, valid columns valuse are between startcol and endcol  
		startcol = Start_col(irow,shortl,longl,BaseData.maxLoopFlag)
		endcol = end_col((irow + 1),Structure_rules) # - (irow == 0)
		# Start scoring
		for icol in range(startcol,endcol):
			# startrow is the lowest row from the previous colomn, considering possible loops
			startrow = Start_row(irow,shortl,longl)	
			# Geting the score
			mpath[irow][icol] = score_v6(irow,icol,startrow,mpath,seq,BaseData,ScoresP)

	# Get the high score from the last row
	yy = nm.array([x.withLoops.Score for x in mpath[row-1]])      
	highScore = max(yy)

	return mpath,row,col,highScore

# End mpathCalc

def PosScore(AA,Mpos,ScoresP):
	'''
	Calculate the location score NOT including the loop scores   
	tmpScr is the score as a class variable
	total_scr is the sum of tmpScr
	'''
	tmpScr = ScoreFunction()
	#tmpScr.hydrophobic_scr = AA.Hydropathy * Mpos.PointInOut
	tmpScr.hydrophobic_scr = AA.Hydrophobic * Mpos.PointInOut
	tmpScr.polar_scr = -AA.polar * Mpos.PointInOut
	tmpScr.charge_scr = -AA.charge * Mpos.PointInOut
	tmpScr.size_scr = ProbabilityScore(AA.SideChainVol,Mpos.ResSize,Mpos.ResSizeSigma)
	# tmpScr.size_scr = ProbabilityScore(AA.ResWeight,Mpos.ResSize,Mpos.ResSizeSigma)
	tmpScr.corner_scr = AA.corner  * (Mpos.CornerMarker > 0) 		# AA.corner = 0 at current setup
	try:
		tmpScr.special_scr = AA.SpecialRes[Mpos.Special]
	except:
		tmpScr.special_scr = 0

	tmpScr,total_scr = SumScore(tmpScr,ScoresP)

	return tmpScr,total_scr

# End PosScore

def SumScore(RawScore,ScoresP):
	'''
	Calculates the weighted scores
	Not including loop scores
	'''
	hydrophobic_scr = RawScore.hydrophobic_scr * ScoresP.HydroRatio
	polar_scr = RawScore.polar_scr * ScoresP.PolarRatio
	charge_scr = RawScore.charge_scr * ScoresP.ChargeRatio 
	size_scr = RawScore.size_scr * ScoresP.SizeRatio 
	corner_scr = RawScore.corner_scr * ScoresP.CornerRatio 			# AA.corner = 0 at current setup
	special_scr = RawScore.special_scr * ScoresP.SpecialRatio

	# Calculate the weighted scores
	tmpScr = ScoreFunction()
	tmpScr.hydrophobic_scr = hydrophobic_scr
	tmpScr.polar_scr = polar_scr
	tmpScr.charge_scr = charge_scr
	tmpScr.size_scr = size_scr
	tmpScr.corner_scr = corner_scr
	tmpScr.special_scr = special_scr

	tot1 = hydrophobic_scr + polar_scr + charge_scr
	tot2 = size_scr + corner_scr + special_scr
	tot = tot1 + tot2

	return tmpScr,tot

# End SumScore

def SideChainScore(AA1,AA2,ScoresP):
	'''
	Test if side bond is possible:

	sideBond :
	's' - disulfide bond
	'h' - hydrogen bond

	returns-
	BondScore:
	0 no bond
	1 we have hydrogen bond
	1 * s_factor we have disulfide bond

	AA1/2 are amino acids

	n1 and n2 are the numbers that indicate if we can have hydrogen bond and give it a value to
	bond contribution.
	n = -1 when the amino acid (AA) residue has an N donor, short residue.
	n = -2 when the AA residue has an O acceptor, short residue.
	n = -3 when the AA residue has an N donor, long residue that able to bond across two turns.
	n = -5 when the AA residue has an O acceptor, long residue that able to bond across two turns.
	n = -7 when it is a short disulfide  bond
	n = -11 when it is a long disulfide  bond
	n = 0 when bond is not possible.
	n = 1 when this N or O participating in a bond.

	Amino acids that can have side-chain to side-chain bonds (that have N or O in their residue):
	1.	Tryptophan, Trp (W),  (Non-polar, neutral side chain) 			[ N]
	2.	Aspartic acid, Asp (D), (Polar, acidic, hydrophilic)			[O,O]
	3.	Glutamic acid, Glu (E), (Polar, acidic, hydrophilic)			[O,O]
	4.	Serine, Ser (S), (Polar, neutral acidity of side chain)			[O]
	5.	Tyrosine, Tyr (Y), (Polar, neutral acidity, slightly hydrophilic)	[O]
	6.	Asparagine, Asn (N), (Polar, neutral acidity, hydrophilic)		[O,N]
	7.	Glutamine, Gln (Q), (Polar, neutral acidity, hydrophilic)		[O,N]
	8.	Lysine, Lys (K), (Polar, basic side chain, hydrophilic)			[N]
	9.	Arginine, Arg (R), (Polar, strong basic acidity, hydrophilic (strong))	[N,N]
	10.	Histidine, His (H), (Polar, weakly basic acidity, hydrophilic)		[N, N]

	Amino acid participate only in one side-chain bond

	'''
	s_factor = ScoresP.DiS_Hyd_ration		# Disulfide bond is ~10 time stronger than H-bond
	s_ratio = ScoresP.SideBondRatio
	BondScore = 0.
	sideBond = 'n'
	# posible combinations
	s1 = AA1.n1 * AA2.n1
	s2 = AA1.n1 * AA2.n2
	s3 = AA1.n2 * AA2.n1
	s4 = AA1.n2 * AA2.n2
	#Possible results
	HydroBond = [2,5,6]
	DisulfBond = [49,77,121]

	t1 = (s1 in HydroBond) + (s2 in HydroBond) + (s3 in HydroBond) + (s4 in HydroBond) 
	t2 = (s1 in DisulfBond) + (s2 in DisulfBond) + (s3 in DisulfBond) + (s4 in DisulfBond)

	if t2 > 0:
		sideBond = 's'
		BondScore = 1. * s_factor
	if t1 > 0:
		sideBond = 'h'
		BondScore = 1.

	return BondScore * s_ratio

def get_AA(FullThread,AADist):
	'''
	Get the Amino acid that is at AADist backwards from current 
	location at the thread FullThead
	'''
	aaTable = AminoTable()
	aaList = AminoList()
	# we want to go back from current location
	FlipThread = FullThread[::-1]
	CountFlag = (AADist > 0)
	jj = 0 	# Count the steps to the correct amino acid
	ll = len(FlipThread)	# Length of thread
	if CountFlag:
		for i in FlipThread:
			ll -= 1
			jj += 1
			if i == ')': CountFlag = False	# Loop starts
			if i == '(': 			# Loop ends
				CountFlag = True
				AADist += 1
			if i == '|': AADist += 1
			if CountFlag:
				if (AADist > 1) and (ll > 0): AADist -= 1
				else: 
					NewAA = i
					break	# Stop when AADist ==1

	if AADist == 1:
		# Getting the Amino Acid
		if NewAA in aaList: AA = aaTable[NewAA]
		# if we are exactly above a corner cut put 'G'
		else: AA = aaTable['G']
	else: AA = aaTable['G']	
	return AA
# End function get_AA

def threadSidescore(AA,thread,AADist,ScoresP):
	'''
	AA is the amino acid we are adding the the thread
	thread is the existing thread
	AADist is the lendth of a turn, in terms of amino acids
	'''
	aminotable = AminoTable()
	AA1 = aminotable[AA] 				# Current amino acid
	AA2 = get_AA(thread,AADist)				# amino acid a turn away
	scr = SideChainScore(AA1,AA2,ScoresP)		# side chain score

	i = 0					# odd, even counter
	while SideChainScore(AA1,AA2,ScoresP) != 0:
		i += 1
		AA1 = AA2
		AA2 = get_AA(thread,(i+1) * AADist)
	# keep score only if i is odd

	scr = scr * (i % 2)
	return scr
# End threadSidescore

def Start_col(irow,shortl,longl,LoopSpace):
	# Evaluating the possible range of colomns for each row
	if irow == 0: startcol = 0
	else:
		td = LoopSpace + longl + 1	# the combined loop and distance between loops length
		col1 = (LoopSpace + 1)*((irow - 1)/td)	# number of full cycles
		col2 = irow%td + (irow%td == 0) * td - longl
		col3 = (col2 < 1)*(1 + longl - (irow%td))
		startcol = col1 + col2 + col3

	# for each row, valid columns valuse are between startcol and endcol  
	return startcol
# End Start_col

def Start_row(irow,shortl,longl):
	# startrow is the lowest row from the previous colomn, considering possible loops
	return (irow-longl)*((irow-longl)>0)+((irow-longl)<1)-1	

# End Start_row





def AddAAtoThread(irow,icol,crow,ccol,AA,mpath,seq,Thread_data,maxLoopFlag):
	'''
	This function adds the amino acid AA to the existing thread which is in Thread_data
	'''  

	if irow !=0:
		impath = mpath[irow][icol]
		iimpath = mpath[crow][ccol]
		loopFlag = ((irow - crow) > 1)

		if impath.CornerMarker == 1:
			if iimpath.CornerMarker == 2: Thread_data.thread += '|'	# No corner cutting
			else: Thread_data.thread += '-|'				# corner cutting type 1
		elif iimpath.CornerMarker == 2: Thread_data.thread += '|-'	# corner cutting type 2
		if loopFlag:		# Add loop 
			Thread_data.thread = Thread_data.thread + '(' + seq[(crow + 1):irow] + ')'
			Thread_data.loopFlag = maxLoopFlag
		else: Thread_data.loopFlag = ((Thread_data.loopFlag - 1) > 0)*(Thread_data.loopFlag - 1)
	Thread_data.thread += AA   

	return Thread_data

# end AddAAtoThread

def ThreadsScoreTest(Thread_data,highScore):
	'''
	Delete all threads with less than highScore side chain bond
	'''
	NewThreaddataList = []
	tmp = Thread_data.threadData
	for i in tmp:
		if i.score == highScore: NewThreaddataList += [i]
	# NewThreaddataList = uniqueList(NewThreaddataList)
	Thread_data.threadData = NewThreaddataList
	return Thread_data

# End ThreadsScoreTest

def uniquethreadList(list1):
	'''
	Find a uniq list for list of variables from ThreadData() class
	'''
	list2 = []
	while list1 != []:
		clst = list1[0]
		list2 += [clst]
		tmpList = []
		for lst in list1:
			flg1 = (lst.thread != clst.thread)
			if flg1: tmpList += [lst]
		list1 = cp.copy(tmpList)
	return list2
# End uniquethreadList

def uniqueList(list1):
	'''
	Find a uniq list 
	'''
	list2 = []
	while list1 != []:
		clst = list1[0]
		list2 += [clst]
		tmpList = []
		for lst in list1:
			if lst != clst: tmpList += [lst]
		list1 = cp.copy(tmpList)
	return list2
# End uniquethreadList

def ProbabilityScore(x,DataMean,Sigma):
	'''
	This function uses normal distribution to evaluate the ratio
	between the probability of x with a given mean and sigma

	we use harmonic potential that is very similar to noraml distribution well 
	'''
	if Sigma != 0: p1 = nm.exp(-((x - DataMean)**2)/(2.0 * (Sigma**2)))
	else: p1 = 0.0


	return p1
# end of Probability score function


def score_v6(irow,icol,startrow,mpath,seq,BaseData,ScoresP):
	'''
	Returning the a cell of mpath which is updated with the score, threads and loop counter information
	for an amino acid at a specific location, taking into account the previous best score
	'''

	newMpath = cp.deepcopy(mpath[irow][icol])	# Creating a copy of the current mpath cell     
	AAcid = seq[irow]				# Current amino acid letter
	aminotable = AminoTable()
	seq_list = [aminotable[x] for x in seq]		# list of amino acid objects

	# Check if corner cutting is allowed. True when corner cutting allowed   
	if newMpath.CornerFlag == 3: cornercutflag = True	# corner cut allowed
	else: cornercutflag = False				# corner cut not allowed

	# getting this location score
	AA = seq_list[irow]
	# tmpScr is the score as a class variable, tot is the sum of tmpScr
	# this is the score for the AA in the current position
	tmpScr,totPosScore = PosScore(AA,newMpath,ScoresP)
	# Normelizing the score to the number of AA
	#tmpScr = tmpScr/(irow +1.0)
	#totPosScore = totPosScore/(irow +1.0)

	# Creating mat1 and mat2. mat2 is used when corner cuting is allowed
	# mat1/2 are used to collect the cells information mpath
	mat1 = nm.array([])
	dat1 = []
	mat2 = nm.array([])
	dat2 = []
	mat3 = nm.array([])
	dat3 = []  

	if irow > 0:
		# creating mat1/2 for 2nd to last rows
		for i in range(startrow,(irow-1)): 
			# Collecting score2's up to irow-2
			# getting all the side chain score with the current amino acid      
			mat1,dat1 = Score_collection_v3(mat1,dat1,irow,icol,i,(icol-1),totPosScore,
						                    AAcid,AA,ScoresP,BaseData,seq,mpath,True)
			if cornercutflag:		# if corner cutting allowed
				mat2,dat2 = Score_collection_v3(mat2,dat2,irow,icol,i,(icol-2),totPosScore,
								                AAcid,AA,ScoresP,BaseData,seq,mpath,True)
				# if corner cutting is not allowed - create equal length mat2 with bad scores
				# to ensure they will never be selected
			else:
				mat2 = nm.append(mat2,-1300)
				dat2 = dat2 + [MatrixPosList()]

		# Creating mat1/2 for the 1st row
		# 1st raw is allowed to have history with loops. Loop flag need not be zero

		mat1,dat1 = Score_collection_v3(mat1,dat1,irow,icol,(irow-1),(icol-1),totPosScore,
				                        AAcid,AA,ScoresP,BaseData,seq,mpath,False)

		if cornercutflag:		# if corner cutting allowed
			mat2,dat2 = Score_collection_v3(mat2,dat2,irow,icol,(irow-1),(icol-2),totPosScore,
						                    AAcid,AA,ScoresP,BaseData,seq,mpath,False)
		else:
			# if corner cutting is not allowed - create equal length mat2 with bad scores
			# to ensure they will never be selected
			mat2 = nm.append(mat2,-1300)
			dat2 = dat2 + [MatrixPosList()]

		# Collecting the no-loop data score1 : irow-1
		mat3,dat3 = Score_collection_v3(mat3,dat3,irow,icol,(irow-1),(icol-1),totPosScore,
				                        AAcid,AA,ScoresP,BaseData,seq,mpath,True)
		if cornercutflag:
			mat3,dat3 = Score_collection_v3(mat3,dat3,irow,icol,(irow-1),(icol-2),totPosScore,
						                    AAcid,AA,ScoresP,BaseData,seq,mpath,True)

	else:
		tmp = MatrixPosList()
		tmp.threadData[0].score = totPosScore
		tmp.threadData[0].thread = AAcid
		mat1 = nm.append(mat1,0)
		dat1 += [tmp]
		mat2 = nm.append(mat2,-1500)
		dat2 += [tmp]
		mat3 = nm.append(mat3,0)
		dat3 += [tmp]

	lmat = len(mat1)
	xi = lmat - BaseData.shortl
	xf = lmat-1
	if xi < 2:	# is the shortes loop is longer than mat
		xi = 1
		xf = lmat
	mat1[xi:xf] = -1000   	# assigning low score to forbiden threading options
	mat2[xi:xf] = -1000   	# assigning low score to forbiden threading options


	# calculating the total average high score
	mat = nm.append(mat1,mat2)	
	mat = nm.append(mat,mat3)

	dat = dat1 + dat2 + dat3
	hiscr = mat.max()		# get the highest value in 'mat'
	resultsDat = []


	for i in range(len(mat)):
		if abs(mat[i] - hiscr)<= BaseData.tolerance: 
			resultsDat += dat[i].threadData 	# collecting all treads that are within the tolerance

	resultsDat = uniquethreadList(resultsDat)


	# Collecting the threads with high score and no loop history from dat
	NoLoopDat1 = []	# no-loop threads within resultDat
	NoLoopDat2 = []	# no-loop threads within dat3/mat3 that have local (mat3) max score
	NoLoopDat = []	# general no-loop thread results

	mat3Score = mat3.max()		# the highest score from (row-1),(col-1)

	for i in resultsDat:		
		if i.loopFlag == 0:
			NoLoopDat1 += [i]

	for i in range(len(mat3)):		# including results with loop history
		if abs(mat3[i] - mat3Score)<= BaseData.tolerance:
			NoLoopDat2 += dat3[i].threadData

	if NoLoopDat1 == []:		# there are no no-loop threads with high score
		NoLoopDat = cp.copy(NoLoopDat2)
	else:
		NoLoopDat = cp.deepcopy(NoLoopDat1)
		mat3Score = hiscr 

	NoLoopDat = uniquethreadList(NoLoopDat)

	# score2 final value
	scr2 = hiscr
	scr1 = mat3Score 
	# setting the values in newMpath

	newMpath.NoLoop.threadData = cp.deepcopy(NoLoopDat)
	newMpath.NoLoop.Score = scr1
	newMpath.withLoops.threadData = cp.deepcopy(resultsDat)
	newMpath.withLoops.Score = scr2
	newMpath.UsageFlag = 0

	return newMpath

# End of score function

def Score_collection_v3(mat,dat,irow,icol,CurrentRow,CurrnetCol,totPosScore,
                        AAcid,AA,ScoresP,BaseData,seq,mpath2,LoopHistory):
	'''
	this function updates mat, an array with the good scores and dat, a list of mpath 
	cells coresponding to the scores

	mat, dat will be updated
	irow,icol is the current position in mpath
	currentRow/currentCol is the position we moving from, to irow/icol
	AAcid is the amino acid we are adding
	turnlenght is the number of amino acid in a full turn
	ScoresP is the scores ratio between the scoring parameters
	maxLoopFlag is the minimum number of steps we need to take between loops
	seq is the sequence in a string form
	mpath - the scoring matrix
	LoopHistory - flag that is true if we using a loop - if we use a loop we should look at 
	the NoLoop data
	'''
	aminotable = AminoTable()
	maxLoopFlag = BaseData.maxLoopFlag		# number of step needed between loops   
	maxCornerFlag = BaseData.maxCornerFlag	# number of step needed between corners  

	if LoopHistory: tmp = cp.deepcopy(mpath2[CurrentRow][CurrnetCol].NoLoop)
	else: tmp = cp.deepcopy(mpath2[CurrentRow][CurrnetCol].withLoops)

	if AAcid == 'C': turnLenght = mpath2[irow][icol].DiSulDis
	else: turnLenght = mpath2[irow][icol].Turn	
	highScore = -10000.0 		# collect the highest score
	newmat = nm.array([])
	newdat = []
	for i in range(len(tmp.threadData)):
		LoopLength = 0
		LoopScore = 0
		CornerCutFlag = (abs(icol - CurrnetCol)>1)			# check if should apply corner cut penalty
		t = cp.deepcopy(tmp.threadData[i])	
		tempThread = t.thread + '-' * CornerCutFlag			# add '-' if we cut corner
		newScore = threadSidescore(AAcid,tempThread,turnLenght,ScoresP)	# add the side chain score	
		newScore += (irow>0)*mpath2[CurrentRow][CurrnetCol].UsageFlag 	# making sure to use valid paths
		loopflag = (abs(irow - CurrentRow)>1)				# do we have a loop
		if loopflag:		# Calculating the score of amino acid on the loop
			loopAA_hydro = [aminotable[x].Hydrophobic for x in seq[(CurrentRow+1):irow]]
			loopAA_charge = [aminotable[x].charge for x in seq[(CurrentRow+1):irow]]
			loopAA_polar = [aminotable[x].polar for x in seq[(CurrentRow+1):irow]]
			loopAB = [aminotable[x].loop for x in seq[(CurrentRow+1):irow]] 	# Adding special loop score (this term is 0 in current version)
			LoopScore = ScoresP.SpecialRatio * sum(loopAB)						# loopAB is 0 so this term is zero as well
			LoopScore -= (sum(loopAA_hydro) * ScoresP.HydroRatio * ScoresP.LoopWeight)
			LoopScore += (sum(loopAA_polar) * ScoresP.PolarRatio * ScoresP.LoopWeight)
			LoopScore += (sum(loopAA_charge) * ScoresP.ChargeRatio * ScoresP.LoopWeight)
			LoopScore -= ScoresP.LoopPenalty / 3.0
			LoopLength = float(irow - CurrentRow - 1)


		loopPenaltyFlag = loopflag*(mpath2[irow][icol].CornerMarker in [0,2]) 	# check if should apply loop penalty
		newScore -= loopPenaltyFlag * ScoresP.LoopPenalty			# applying loop penalty
		newScore -= CornerCutFlag * ScoresP.CornerCutPenalty
		newScore += totPosScore
		t = AddAAtoThread(irow,icol,CurrentRow,CurrnetCol,AAcid,mpath2,seq,t,maxLoopFlag)
		newScore += AA_VolumeScore(t.thread,BaseData)*ScoresP.TurnVolPenaly
		t.score = (t.score * (CurrentRow + 1) + LoopScore + newScore)/(irow + 1.0) 	# Keeping the score normelized
		if t.score > highScore : highScore = cp.copy(t.score)	# tracking the highest score
		#t = AddAAtoThread(irow,icol,CurrentRow,CurrnetCol,AAcid,mpath2,seq,t,maxLoopFlag)
		tmp.threadData[i] = cp.deepcopy(t)

	tmp.Score = highScore
	newmat = nm.append(mat,tmp.Score)
	newdat = dat + [tmp]

	return newmat,newdat

# end Score_collection

def PlotScoreMatrix (mpath,ProtName):
	'''
	Plotting the scoring matrix:  score vs. col at the last row
	'''
	top_row = nm.size(mpath,0)-2
	col = nm.size(mpath,1)
	yy = [x.withLoops.Score * (x.withLoops.Score > 0) for x in mpath[top_row]]
	y = nm.array(yy)
	x = nm.array(range(col))

	# Plotting
	TitleString = 'Matrix Scores for ' + ProtName
	xLabelStr = 'Scoring Matrix Column'
	yLabelStr = 'Scores per amino acid'

	pl.title(TitleString)
	pl.xlabel(xLabelStr)
	pl.ylabel(yLabelStr)
	pl.plot(x,y)
	pl.ylim([0,110])
	pl.show()

# end PlotScoreMatrix 

def LHBqualityPanelty(Thread,seq):
	'''
	this function return a penalty to the score in relation to the number of
	abnormalities compare to LHB structure.
	This is to eliminate false positive identification when scanning protein files
	'''
	seqLength = len(seq)
	nSides = len([x for x in Thread if x == '|']) + 1
	nDefects = len([x for x in Thread if (x == '(') or (x == '-')])
	f = 0.0

	if nSides > 2:
		f = nSides*(nm.exp((nDefects/float(nSides))**3)-1.0)
	if f < 1.5 : f = 0	# Ignor small abnormalities 


	return f
# end LHBqualityPanelty

def get_seq_score(Structure_rules,Partial_seq,BaseData,ScoresP,tol=0):
	'''
	Get the high score form the matrix last row
	'''
	mpath,row,col,highScore = mpathCalc(Structure_rules,Partial_seq,BaseData,ScoresP)  
	Threa_list = Get_thread_list(mpath,row,col,highScore,tol)
	threadlist = [x.thread for x in Threa_list]
	return highScore,threadlist

# end get_seq_score

def Seq_file_score_LBH18(seq,Structure_rules,ScoresP,shortl=1,longl=2,pseq_length=72,minLoopDis=17,tol=0,Hi_score=60):
	'''
	This function scans a protein sequence, seq, for possible LHB-18 regions
	It outputs (x,scores_list)
	where x is the amino acid number and the scores_list is the coresponding scores

	shortl=1,longl=2,minLoopDis=16 : default parameters for scoring
	pseq_length is length of the of sub sequence on which we run the scoring algorithm
	Hi_score is the score above which we want to collect the threads
	'''

	# init
	BaseData = BaseInputData(Structure_rules)
	BaseData.shortl = shortl
	BaseData.longl = longl
	BaseData.maxLoopFlag = minLoopDis
	build_Reverse_in_sidechain_list(BaseData,Structure_rules)
	thread_list = []
	ScoreThreadList = []

	# parce the sequence
	side_length = 6				# side length in the LHB structure model
	n_iter = (len(seq)- pseq_length)/side_length + 1	# number of iterations
	scores_array = nm.zeros(n_iter)		# Array for results collection
	x = nm.array(range(n_iter))

	for i in x:
		ii = i * side_length
		p_seq = seq[ii:(ii+pseq_length)]	# sending two turns for scoring
		score,threads = get_seq_score(Structure_rules,p_seq,BaseData,ScoresP,tol)     
		scores_array[i] = score			# Collect the scores
		if score > Hi_score:
			ScoreThreadList += [[ii+1,threads,score]]
	# Evaluating LHB regions
	x2 = nm.array([])
	for i in scores_array:
		if i>70:x2 = nm.append(x2,100)
		else:x2 = nm.append(x2,10)
	x *= side_length		# addjusting the valuse of x to represent amino acid number

	return x,scores_array,ScoreThreadList
# end of Seq_file_score_LBH18

def Get_thread_list(mpath,row,col,highscore,tol):
	# getting the results of the top scoring threads
	irow = row -1
	threadlist1 = []
	ScoreList = []
	for icol in range(col): ScoreList += [mpath[irow][icol].withLoops.Score]

	for icol in range(col):
		if abs(ScoreList[icol] - highscore) <= tol:
			threadlist1 += mpath[irow][icol].withLoops.threadData

	threadlist1 = uniquethreadList(threadlist1)  
	return threadlist1
#end of Get_thread_list

def Convert_to_turn_volume(TurnResWeight):
	'''
	Converting the information of volume per side to volume per turn.

	Note that there is an overlapping between the sides
					v1 = s1+s2+s3, v2 = s2+s3+s4 ...
	'''
	tempArray = cp.copy(TurnResWeight)
	nn = len(tempArray)-2
	newArray = nm.zeros(nn)
	for i in range(nn):
		newArray[i] = nm.sum(tempArray[i:(i+3)])
	return newArray
# end of Convert_to_turn_volume

def build_Reverse_in_sidechain_list(BaseData,Structure_rules):
	'''
	This function create a list of all side chain that are pointing into the structure.

	it is used when calculating full turn volume.

	We reverse the order because we are interested only at the last turn and we reverse
	the sequence when calculating that volume

	when we send a list of AA for a volume check, the first AA is in the position of the first AA
	with side chain pointing in. so we need to change the list of side chain pointing in accordingly
	'''
	l = BaseData.base_length + 1
	ll = BaseData.base_length
	sl = BaseData.base_length
	rev_in_list = []			# full turn, reverse order, inward side chain position
	n = Structure_rules[0].Turn / sl	# number of unique side in a full turn
	for i in range(n):
		for p in Structure_rules:
			if p.PointInOut == 1: rev_in_list += [l - p.PosNum + ll*i]
	rev_in_list.sort()
	i = rev_in_list[0]
	rev_in_list = [x - i  for x in rev_in_list]		# substruct the smalest term
	BaseData.Reverse_in_sidechain = rev_in_list

	in_list = []			# # one side, inward side chain position
	for p in Structure_rules:
		if p.PointInOut == 1: in_list += [p.PosNum]
	in_list.sort()
	BaseData.Side_in_sidechain = in_list

# End build_Reverse_in_sidechain_list

def AA_VolumeScore(thread,BaseData):
	'''
	When the new Amino Acid is in the first pointing in position, we need to check the turn volume contribution
	'''
	base_length = BaseData.base_length
	ResPointIn = BaseData.Side_in_sidechain
	DisPosIn = BaseData.Side_in_sidechain[1] - BaseData.Side_in_sidechain[0]
	RevPosList = BaseData.Reverse_in_sidechain
	meanVol = float(BaseData.TurnVolume)
	SigmaVol = BaseData.TurnVolSigma  
	tl = BaseData.TurnLength
	aminotable = AminoTable()
	r_thread = thread[::-1]
	x1 = []			# Collect reverse list without corner marks
	x2 = []			# Collect reverse list with corner marks
	l = 0			# thread length
	# Collect the thread of a full turn without loops
	NotLoopFlag = True
	for a in r_thread:
		if a == ')' : NotLoopFlag = False
		if a =='|' : NotCornerFlag = False
		else: NotCornerFlag = True
		if NotLoopFlag: 
			x2 += [a]
			if NotCornerFlag: 
				x1 += [a]
				l += 1
		if a =='(' : NotLoopFlag = True
		if l == tl: break
	x2 = st.join(x2,'')			# this is the reverse string without loops
	if x2 == '':
		icol = -1
	elif x2[0] == '-' : icol = -1	# x2 starts with corner cut
	else: icol = x2.find('|')		# the location of the first |

	TotVol = 0.				# total volume, full turn
	VolScr = 0.				# Volume penalty contribution
	# check if we should calculate score
	if (icol in ResPointIn):		# we only taking contribution from first 
		if len(x1) == tl: 	# making sure we have full turn
			for i in RevPosList:
				if x1[i] != '-':
					TotVol += aminotable[x1[i]].SideChainVol
			i = 0.
			for x in x1:
				if x == '-': i += 1
			TotVol += i * 27.		# compensating for corner cuts
			#VolScr = nm.exp(-((TotVol - meanVol)**2)/(2.0 * (SigmaVol**2))) - 1.
			VolScr = -(TotVol/meanVol - 1.)**2 

	return VolScr

# End AA_VolumeScore 

def convertPathNames(pathName):

	newfilePath = ''
	for i in pathName:
		if i == '\\':
			newfilePath += '//'
		else:
			newfilePath += i

	newfilePath += '//'

	return newfilePath

# End convertPathNames

def ThreadSeq(seq):
	'''
	This function is identical to the betatesting just in a function form
	seq : the sequence to be threaded
	seq should be in this format -  'LCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQGTHSQWNKPSKPKTNMKHMAGAA'
	'''
	# Dictionary of amino acids
	aminotable = AminoTable()	

	# scoring parameters ratios
	ScoresP = ScoreParameters()

	# Read the protein sequence
	print seq
	print 'Sequence length = ',len(seq)
	# End of reading the seq
	StartTime = time.clock()
	# Read structure rules
	Structure_rules,s_file_path,s_file_name = structure_file_info()
	print 'Finish reading structure rules \n'
	# End structure rules reading

	# read some program restrictions
	print 'Enter data'
	# Read the structure information
	BaseData = readProgParameters(Structure_rules)
	# create a list of all side chain that are pointing into the structure in BaseData
	build_Reverse_in_sidechain_list(BaseData,Structure_rules)

	# print warning  - sequence too short    
	if len(seq)<35:
		print '************  sequence too short ***************' 
	# ##############################################
	# Get the score matrix 		
	print 'Start scoring and thread calculation \n'
	mpath,row,col,highscore = mpathCalc(Structure_rules,seq,BaseData,ScoresP)
	print 'End Scoring \n'
	# getting the results of the top scoring threads
	irow = row -1
	threadlist1 = []
	ScoreList = []
	for icol in range(col): ScoreList += [mpath[irow][icol].withLoops.Score]

	for icol in range(col):
		if abs(ScoreList[icol] - highscore) <= BaseData.dScore:
			threadlist1 += mpath[irow][icol].withLoops.threadData

	threadlist1 = uniquethreadList(threadlist1)       
	dtime = (time.clock() - StartTime)
	print 'Calculation time is %.2f sec' % dtime
	if len(threadlist1)>0:
		# Printing results
		result_output_v4(seq,highscore,threadlist1,BaseData.dScore)
	else:
		print 'No Threads, Program Error'

	ttl = file_name[0:5]
	PlotScoreMatrix (mpath,ttl)
	# Return nothiong, just print results
	return

# End ThreadSeq

def StructureRulls(StructureName):
	'''
	This function return the scoring rules used in structure testing software

	Point in: 1
	Point out: -1

	CHANGES ********************************************************
	Jul-22-2009: By Youval: Change residue weight to residue volume
	****************************************************************

	the data is in an matrix, the 11 columns in the array are:

	Position,Point	in/out,Can-cut corner,Corner marker,Residue	prefered size,
	sigma,Turn	Residue Volume,Sigma,Special feature,Turn length,Disulfide lenght	

	The number of rows in the matrix represents the basic structure repeated unit length
	LBH-18: 6 rows
	LBH-15: 5 rows

	'''
	if StructureName == 'LBH-18':
		DataMatrix = [[1,-1,0,1,0,0,487.23,61.48,1,18,0],
				      [2,-1,3,3,0,0,487.23,61.48,2,18,36],
				      [3,1,0,0,63.83,30.09,487.23,61.48,3,18,0],	
				      [4,-1,0,0,0,0,487.23,61.48,4,18,0],
				      [5,1,0,0,99.14,21.22,487.23,61.48,5,18,15],
				      [6,-1,1,2,0,0,487.23,61.48,6,18,0]]
	elif StructureName == 'LBH-15':
		DataMatrix = [[1,-1,3,1,0,0,405.47,50.19,7,15,0],
				      [2,1,0,0,45.44,15.97,405.47,50.19,8,15,3],
				      [3,-1,0,0,0,0,405.47,50.19,9,15,0],
				      [4,1,1,0,90.29,31.35,405.47,50.19,10,15,12],
				      [5,-1,0,2,0,0,405.47,50.19,11,15,0]]

	else:
		DataMatrix = [[0,0,0,0,0,0,0,0,0,0,0]]

	tmp = []
	for line in DataMatrix[:]:
		a = ThreadStructure()
		a.PosNum = int(line[0])
		a.PointInOut = float(line[1])
		a.CornerFlag = int(line[2])
		a.CornerMarker = int(line[3])
		a.ResSize = float(line[4])
		a.ResSizeSigma = float(line[5])	
		a.TurnVolume = float(line[6])
		a.TurnVolSigma = float(line[7])
		a.Special = int(line[8])
		a.Turn = int(line[9])	
		a.DiSulDis = int(line[10])
		tmp = tmp + [a]

	return tmp
# End of structure data function