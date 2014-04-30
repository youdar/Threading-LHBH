class ScoreParameters:
	def __init__(self):
		self.HydroRatio = 60.			# [55,55.]
		self.PolarRatio = 35.+ 50.		# 90[21,25,35,35.]
		self.ChargeRatio = 88.			# 50[85,50,88,95.]
		self.SizeRatio = 175.			# 
		self.DiS_Hyd_ration = 30.0		# [30.]not varied in minuit
		self.SpecialRatio = 110.0		# [110.]not varied in minuit
		self.CornerRatio = 60.0			# [60.]not varied in minuit
		self.SideBondRatio = 35.0		# [60,80,60]
		self.LoopPenalty = 500.			# 650[162,280,300.]penalty for loops when not at the corner
		self.CornerCutPenalty = 62.0	# [97,120.]penalty for cutting corners
		self.TurnVolPenaly = 250.0		# [60,420,430.]penalty for turn volume too small or too big
		self.LoopWeight = 0.75			# Weight of hydrophobicy when in a loop vs. in the LBH

# end class ScoreParameters

class ThreadStructure:
	def __init__(self,Special=0):
		self.PosNum = 0
		self.PointInOut = 0			# Side chain points out: -1 point in: 1 
		self.ResSize = 0			# mean of residue weight of amino acid in this position
		self.ResSizeSigma = 0			# standard diviation from that mean
		self.TurnVolume = 0			# mean of residues weight in a full turn
		self.TurnVolSigma = 0			# Standad diviation from that mean
		self.CornerFlag = 0			# 0:Not a corner, 1:can cut a corner, 3:after corner cut 
		self.CornerMarker = 0			# 0:not a corner, 1:side start, 2:side end
		self.Turn = 0				# Full turn length
		self.DiSulDis = 0			# Distance for disulfide bond
		self.Special = Special			# Uses dictionary from AA special instructions 
												# a number that will be matched with a score

class ThreadData:
	def __init__(self):
		self.score = 0
		self.loopFlag = 0			# count steps till loops are allowed again
		self.cornerflag = 0			# count steps till corner cut is allowed again
		self.thread = ''

class MatrixPosList:
	def __init__(self,row=0,col=0,loopFlag=0,Score=0,threads=[]):
		self.Score = Score					# Store best score
		self.threadData = [ThreadData()]	# list of threads with threads and loopflag data
		self.threadsScores = []				# used when updating scores

class ScoreMatrix:
	def __init__(self,Special=0):
		self.NoLoop = MatrixPosList()		# Best score,pos without loop
		self.withLoops = MatrixPosList()	# Overall Best score
		self.UsageFlag = -5550				# 0:used, -5550:not used (score panelty)
		self.CornerFlag = 0					# 0:Not a corner, 1:can cut a corner, 3:after corner cut
		self.CornerMarker = 0				# 0:not a corner, 1:side start, 2:side end
		self.PointInOut = 0					# 1: point in, -1: point out
		self.ResSize = 0					# mean of residue size of amino acid in this position
		self.ResSizeSigma = 0				# standard diviation from that mean
		self.PosNum = 0						# the relative position in the characteristic structure
		self.Turn = 0						# Full turn length
		self.DiSulDis = 0					# Distance for disulfide bond
		self.Special = Special				# Uses dictionary from AA special instructions 
											# a number that will be matched with a score
class ScoreFunction:
	def __init__(self):
		self.hydrophobic_scr = 0.0
		self.polar_scr = 0.0
		self.charge_scr = 0.0
		self.size_scr = 0.0
		self.corner_scr = 0.0
		self.special_scr = 0.0 
		self.loop_scr = 0.0
	def __div__(self,n):
		nn = float(n)
		self.hydrophobic_scr = self.hydrophobic_scr/nn
		self.polar_scr = self.polar_scr/nn
		self.charge_scr = self.charge_scr/nn
		self.size_scr = self.size_scr/nn
		self.corner_scr = self.corner_scr/nn
		self.special_scr = self.special_scr/nn 
		self.loop_scr = self.loop_scr/nn
		return self

class pathstack:

	def __init__(self,row=-1,col=-1,list1=[],list2=[],nPath=[],pointer=0,loopcounter=0):
		self.pos = [row,col]
		self.list1 = list1
		self.list2 = list2
		self.nPath = nPath
		self.pointer = pointer
		self.loopcounter = loopcounter

	def GetNextPathBranch(self):
		n = len(self.nPath)
		if n > (self.pointer + 1):
			flg = False	# stop looking
			self.pointer += 1	        
		else:
			flg = True	  
		return flg
	# End GetNextPathBranch

	def GetNextPath(self):
		'''
		Get the next point in the thread, if exist
		'''
		n = len(self.nPath)
		if n > 0:
			flg = True	# Can get another path
			[row,col] = self.nPath[self.pointer]	    
		else:
			flg = False
			[row,col] = [-1,-1]
			print 'No Path is possible from this point : ',self.pos	
		return flg,row,col
	# End GetNextPath


	def Empty(self):
		return (self.nPath != [])

# Class End

class BaseInputData:
	def __init__(self,Structure_rules):
		self.shortl = 0
		self.longl = 0
		self.dScore = 0
		self.tolerance = 0
		self.maxLoopFlag = 0					# Count the steps till next loop is allowed
		self.maxCornerFlag = 0					# Count the steps till next corner cut is allowed
		self.Allowednumpath = 0
		self.TurnVolume = Structure_rules[1].TurnVolume	# mean of residues weight in a full turn
		self.TurnVolSigma = Structure_rules[1].TurnVolSigma	# Standad diviation from that mean
		self.TurnLength = Structure_rules[1].Turn
		self.base_length = len(Structure_rules)			# For a triangle structure, this will the side length
		self.Side_in_sidechain = []				# the list of AA position with inward side chain
		self.Reverse_in_sidechain = []				# the number of the inward side chains, in full turn, in reverse sequnce

# end class BaseInputData
