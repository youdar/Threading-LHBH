'''
The main functions of this program are:

1) Thread a protein sequence using:  ThreadSeq(seq)
   seq is string of amino acids in the form  seq = 'SELNIYQYGGGNSALALQTDARNSDLTITQHGGGNGADVGQGSDDSSIDLTQ'
   You can use more options: ThreadSeq(seq,ScoreTolerance=2,ProteinName = '')
   ScoreTolerance is the threads score range that is being collected in the scoring process
   ProteinName is an optional name that will be added to the output plot

2) Scaning a protein sequnce using: AnalyzeProteinSeq(segLength,LoopLength,seq,StructureSTR)
   StructureSTR is the structure type. The default is 'LBH-18' but you can change it by setting StructureSTR = 'LBH-15'
   seq is string of amino acids in the form  seq = 'SELNIYQYGGGNSALALQTDARNSDLTITQHGGGNGADVGQGSDDSSIDLTQ'
   segLength is an integer representing the test segment length. for example segLength = 54
   LoopLength is the maximum allowed loop length. For example LoopLength = 4

   The output file and folder need to be set in the function before execution.
   for example:
   OutFileName = 'plotdata'
   WinFilePath = 'C:\Users\youval\Documents\Education\Research\Paper 2 - Threading program\Plots'

   The function will create a txt file with the the results of the scan.
   The file path is writen for windows machines, you might need to write the file path in a different
   format if using Mac or Linux machine
   Threads are collected but not saved on the file

3) Score and see score weight distribution for a thread using: ThreadScores(ThreadSTR,StructureSTR)
   ThreadSTR is the thread being scored. It should  be a string like:  ThreadSTR = 'TTTTKG|ENFTET|-DVKMM|ERVVEQ|MCITQY|E(RES)QAYYQ|RGSSM'
   StructureSTR is the structure type. The default is 'LBH-18' but you can change it by setting StructureSTR = 'LBH-15'

   The output of this function is displied on the terminal, it not bying saved on a file


REMEMBER that LBH-15 might need more scoring parameter adjusting

When Scanning a protein sequence, look at the output txt file to see at what amino acid number we start getting good scores.
You can than take the high scoring sequence regions of the protein and thread them

Date of last update - Nov-4-2010
Updated by Youval Dar :  youval@gmail.com
'''

import numpy as nm
import string as st
import pylab as pl
import copy as cp
import time

class ScoreParameters:
    def __init__(self):
        self.HydroRatio = 60.			# [55]
        self.PolarRatio = 35.+ 50.		# Add the HydroRatio [35 + 55]
        self.ChargeRatio = 88.			# 88 [85]
        self.SizeRatio = 175.			# 170 [175]
        self.DiS_Hyd_ration = 30.0		# not varied in minuit
        self.SpecialRatio = 110.0		# not varied in minuit  [110]
        self.CornerRatio = 60.0			# not varied in minuit
        self.SideBondRatio = 35.0		# 30 [35]
        self.LoopPenalty = 500.			# 300 penalty for loops when not at the corner [650]
        self.CornerCutPenalty = 62.0	# 90  penalty for cutting corners [50]
        self.TurnVolPenaly = 250.0		# 300 penalty for turn volume too small or too big [120] 250
        self.LoopWeight = 0.75			# 0.8 Weight of hydrophobicy when in a loop vs. in the LBH [.45] 		

# end class ScoreParameters

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


class ThreadStructure:
    def __init__(self,Special=0):
        self.PosNum = 0
        self.PointInOut = 0			# Side chain points out: -1 point in: 1 
        self.ResSize = 0			# mean of residue weight of amino acid in this position
        self.ResSizeSigma = 0		# standard diviation from that mean
        self.TurnVolume = 0			# mean of residues weight in a full turn
        self.TurnVolSigma = 0		# Standad diviation from that mean
        self.CornerFlag = 0			# 0:Not a corner, 1:can cut a corner, 3:after corner cut 
        self.CornerMarker = 0		# 0:not a corner, 1:side start, 2:side end
        self.Turn = 0				# Full turn length
        self.DiSulDis = 0			# Distance for disulfide bond
        self.Special = Special		# Uses dictionary from AA special instructions 
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


class AminoAcid:
    def __init__(self,name='AA'):
        self.name = name
        self.name3L = ''
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.ResWeight = 0 		# residue weight (weight - 56) backbone weight is 56
        self.ResVol = 0			# Residue volume from http://prowl.rockefeller.edu/aainfo/volume.htm
        self.SideChainVol = 0		# Side Chain volume is evaluated as ResVol - 0.9 Gly.ResVol
        self.Hydropathy = 0		# Hydropathy index
        self.n1 = 0			
        self.n2 = 0 
        # n values
        # -1 when the amino acid (AA) residue has an N donor, short residue.	
        # -2 when the AA residue has an O acceptor, short residue.
        # -3 when the AA residue has an N donor, long residue that able to bond across two turns.
        # -5 when the AA residue has an O acceptor, long residue that able to bond across two turns.
        # -7 when it is a Cystine(C) 
        # 0 when bond is not possible.
        # 1 when this N or O participating in a bond.   
        # A residu can only participate in one side-chain bond. So when a bond is created
        # for example with n1, n1 get the bond value and n2 will be assigned 0


    def __mul__ (self,other):
        # Evaluating side chain interaction
        Prod = self.donor * other.acceptor
        return Prod


# ############   Non Polar, HydroPhobic   ###########
class Ala(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'A')
        # Alanine	 
        # ###
        # CH3-CH(NH2)-COOH
        #
        #
        # Molecular weight 89.09 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.616 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 1.8 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point (when protonation accure) pH 6.01
        # pKa( alpha-COOH) 2.35 
        # pKa( alpha-NH2) 9.87
        # CAS # 56-41-7
        # PubChem ID 5950
        #
        self.name3L = 'ALA'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0	
        self.Hydropathy = 1.8
        self.ResWeight = 33
        self.ResVol = 88.6
        self.SideChainVol = 88.6-54.1


class Val(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'V')
        # Valine
        # #########
        # (CH3)2-CH-CH(NH2)-COOH
        #
        #
        # Essential AA (cannot be synthesized by humans)
        # Molecular weight 117.15 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.825 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 4.2 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.00
        # pKa( alpha-COOH) 2.39
        # pKa( alpha-NH2) 9.74
        # CAS # 72-18-4
        # PubChem ID 1182
        #
        self.name3L = 'VAL'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.Hydropathy = 4.2
        self.ResWeight = 61
        self.ResVol = 140.0
        self.SideChainVol = 140-54.1


class Leu(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'L')
        # Leucine	 
        # #############
        # (CH3)2-CH-CH2-CH(NH2)-COOH
        #
        #
        # Essential AA
        # Molecular weight 131.18 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.943 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 3.8 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.01
        # pKa( alpha-COOH) 2.33
        # pKa( alpha-NH2) 9.74
        # CAS # 61-90-5
        # PubChem ID 6106
        #
        self.name3L = 'LEU'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.Hydropathy = 3.8
        self.ResWeight = 75
        self.ResVol = 166.7
        self.SideChainVol = 166.7-54.1


class Ile(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'I')
        # Isoleucine
        # ###############
        # CH3-CH2-CH(CH3)-CH(NH2)-COOH
        #
        #
        # Essential AA
        # Molecular weight 131.18 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.943 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 4.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.05
        # pKa( alpha-COOH) 2.33
        # pKa( alpha-NH2) 9.74
        # CAS # 61-90-5
        # PubChem ID 6106
        #
        self.Hydropathy = 4.5
        self.ResWeight = 75
        self.name3L = 'ILE'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 166.7
        self.SideChainVol = 166.7-54.1

class Phe(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'F')
        # Phenylalanine
        # ######
        # Ph-CH2-CH(NH2)-COOH
        # The residue Ph-CH2 : C6H5-CH2  benzyl
        #
        # Essential AA
        # Molecular weight 165.19 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 1 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 2.8 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.49
        # pKa( alpha-COOH) 2.20
        # pKa( alpha-NH2) 9.31
        # CAS # 63-91-2
        # PubChem ID 994
        #
        self.Hydropathy = 2.8
        self.ResWeight = 109
        self.name3L = 'PHE'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 189.9
        self.SideChainVol = 189.9-54.1


class Trp(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'W')
        # Tryptophan 
        # ##############
        # Ph-NH-CH=C-CH2-CH(NH2)-COOH
        # |________|
        #
        # contains an indole functional group.
        # aromatic heterocyclic organic compound
        # It has a bicyclic structure, consisting of a six-membered benzene ring fused to 
        # a five-membered nitrogen-containing pyrrole ring
        #
        # Essential AA
        # Molecular weight 204.23 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.878 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -0.9 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.89
        # pKa( alpha-COOH) 2.46
        # pKa( alpha-NH2) 9.41
        # CAS # 73-22-3
        # PubChem ID 6305
        #
        self.Hydropathy = -0.9
        self.ResWeight = 148
        self.name3L = 'TRP'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 227.8
        self.SideChainVol = 227.8-54.1

class Met(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'M')
        # Methionine 
        # ############
        # CH3-S-(CH2)2-CH(NH2)-COOH
        # sulfur-containing residue
        # methyl donor  R-CH3
        # methionine is incorporated into the N-terminal position of all proteins 
        # in eukaryotes and archaea during translation, although it is usually removed 
        # by post-translational modification
        #
        # Essential AA
        # Molecular weight 149.21 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.738 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 1.9 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.74
        # pKa( alpha-COOH) 2.13
        # pKa( alpha-NH2) 9.28
        # CAS # 63-68-3
        # PubChem ID 876
        #
        self.Hydropathy = 1.9
        self.ResWeight = 93
        self.name3L = 'MET'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 162.9
        self.SideChainVol = 162.9-54.1


class Pro(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'P')
        # Proline
        # **********
        # NH-(CH2)3-CH-COOH
        # |_________|
        # Side chain bond to C alpha
        # exceptional conformational rigidity 
        # usually solvent-exposed.
        # lacks a hydrogen on the amide group, it cannot act as a hydrogen bond donor, 
        # only as a hydrogen bond acceptor.
        #
        # Molecular weight 115.13 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.711 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -1.6 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.30
        # pKa( alpha-COOH) 1.95
        # pKa( alpha-NH2) 10.64
        # CAS # 147-85-3
        # PubChem ID 614
        #
        self.Hydropathy = -1.6
        self.ResWeight = 59
        self.name3L = 'PRO'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0,3:-3,4:-3,5:-3,6:-2}	# special value scores
        self.n1 = 0			
        self.n2 = 0 
        self.ResVol = 112.7
        self.SideChainVol = 112.7-54.1

# ############ Non Polar Uncharged	###########

class Gly(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'G')	
        # 
        # NH2-CH2-COOH
        # 
        # Molecular weight 75.07 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.501 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -0.4 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.06
        # pKa( alpha-COOH) 2.35
        # pKa( alpha-NH2) 9.78
        # CAS # 56-40-6
        # PubChem ID 750
        self.Hydrophobic = 0		# 1: Hydrophobic, 0: Hydrophilic
        self.Hydropathy = -0.4
        self.ResWeight = 19
        self.name3L = 'GLY'
        self.SpecialRes = {0:0,3:-6,5:-6}	# special value scores
        self.ResVol = 60.1
        self.SideChainVol = 60.1-54.1

# ############ Polar Uncharged	###########

class Ser(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'S')
        # Serine
        # ######
        # HO-CH2-CH(NH2)-COOH
        # 
        # Molecular weight 105.09 Da
        # Ploar
        # Acidity - Natural
        # Hydrophobicity 0.359 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -0.8 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.68
        # pKa( alpha-COOH) 2.19
        # pKa( alpha-NH2) 9.21
        # CAS # 56-45-1
        # PubChem ID 617
        #
        self.Hydropathy = -0.8
        self.ResWeight = 49
        self.name3L = 'SER'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 89.0
        self.SideChainVol = 89-54.1


class Thr(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'T')
        # Threonine
        # ##########
        # CH3-CH(OH)-CH(NH2)-COOH
        # bearing an alcohol group
        #
        # Essential AA
        # Molecular weight 119.12 Da
        # Ploar
        # Acidity - Natural
        # Hydrophobicity 0.450 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -0.7 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.60
        # pKa( alpha-COOH) 2.09
        # pKa( alpha-NH2) 9.10
        # CAS # 72-19-5
        # PubChem ID 6288
        #
        self.Hydropathy = -0.7
        self.ResWeight = 63
        self.name3L = 'THR'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 116.1
        self.SideChainVol = 116.1-54.1

class Cys(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'C')
        # Cysteine
        # ######
        # HS-CH2-CH(NH2)-COOH
        # thiol (R-S-H) side chain 
        # Has Sulfur in side chain
        #
        # Molecular weight 121.16 Da
        # Ploar
        # Acidity - Natural
        # Hydrophobicity 0.680 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 2.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.05
        # pKa( alpha-COOH) 1.92
        # pKa( alpha-NH2) 10.70
        # CAS # 59-90-4
        # PubChem ID 5862
        #
        self.Hydropathy = 2.5
        self.ResWeight = 65
        self.name3L = 'CYS'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.n1 = -7			
        self.n2 = 0
        self.SpecialRes = {0:0}	# special value scores
        self.ResVol = 108.5
        self.SideChainVol = 108.5-54.1



class Tyr(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'Y')
        # Tyrosine  
        # ###########
        # HO-p-Ph-CH2-CH(NH2)-COOH
        # 
        # Molecular weight 181.19 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.880 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -1.3 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.64
        # pKa( alpha-COOH) 2.20
        # pKa( alpha-NH2) 9.21
        # CAS # 60-18-4
        # PubChem ID 1153
        #
        self.Hydropathy = -1.3
        self.ResWeight = 125
        self.name3L = 'TYR'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 193.6
        self.SideChainVol = 193.6-54.1

class Asn(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'N')
        # Asparagine
        # ##########
        # H2N-CO-CH2-CH(NH2)-COOH
        # N Donor - NH2
        # 
        # has carboxamide as the side chain's functional group(R-CO-NH2)
        # side chain can form hydrogen bond interactions with the peptide backbone
        # often found near the beginning and the end of alpha-helices, 
        # and in turn motifs in beta sheets.
        # 
        # Molecular weight 132.12 Da
        # Ploar
        # Acidity - Natural
        # Hydrophobicity 0.236 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.41
        # pKa( alpha-COOH) 2.14
        # pKa( alpha-NH2) 8.72
        # CAS # 70-47-3
        # PubChem ID 236
        #
        self.Hydropathy = -3.5
        self.ResWeight = 76
        self.name3L = 'ASN'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = -1			
        self.n2 = -2
        self.ResVol = 114.1
        self.SideChainVol = 114.1-54.1


class Gln(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'Q')
        # Glutamine  
        # #############
        # H2N-CO-(CH2)2-CH(NH2)-COOH
        # N Donor - NH2 
        #
        # Molecular weight 146.14 Da
        # Ploar
        # Acidity - Natural
        # Hydrophobicity 0.251 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.65
        # pKa( alpha-COOH) 2.17
        # pKa( alpha-NH2) 9.13
        # CAS # 56-85-9
        # PubChem ID 5950
        #
        self.Hydropathy = -3.5 
        self.ResWeight = 90
        self.name3L = 'GLN'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.n1 = -1			
        self.n2 = -2
        self.SpecialRes = {0:0}	# special value scores
        self.ResVol = 143.8
        self.SideChainVol = 143.8-54.1

#  ##########   Polar Acidic   ###########

class Asp(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'D')
        # Aspartic acid
        # ########
        # HOOC-CH2-CH(NH2)-COOH
        # 
        # Molecular weight 133.10 Da
        # Ploar
        # Acidity - Acidic
        # Hydrophobicity 0.028 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 2.85
        # pKa( alpha-COOH) 1.99
        # pKa( alpha-NH2) 9.90
        # CAS # 56-84-8
        # PubChem ID 5960
        #
        self.Hydropathy = -3.5 
        self.ResWeight = 77
        self.name3L = 'ASP'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 1
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop self.loop = 0
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0,3:-1,5:-1}	# Special characteristic of residue (change to {0:0,2:-1,4:-1} for LBH-15)
        self.n1 = -2			
        self.n2 = 0
        self.ResVol = 111.1
        self.SideChainVol = 111.1-54.1

class Glu(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'E')
        # ###########
        # HOOC-(CH2)2-CH(NH2)-COOH
        # 
        # Molecular weight 147.13 Da
        # Ploar
        # Acidity - Acidic
        # Hydrophobicity 0.043 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 3.15
        # pKa( alpha-COOH) 2.10
        # pKa( alpha-NH2) 9.47
        # CAS # 56-86-0
        # PubChem ID 611
        #
        self.Hydropathy = -3.5 
        self.ResWeight = 91
        self.name3L = 'GLU'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 1
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop self.loop = 0
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0,3:-1,5:-1}		# Special characteristic of residue (change to {0:0,2:-1,4:-1} for LBH-15)
        self.n1 = -2			
        self.n2 = 0
        self.ResVol = 138.4
        self.SideChainVol = 138.4-54.1


# ##############  Polar Basic   #############

class Lys(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'K')
        # Lysine
        # ##########
        # H2N-(CH2)4-CH(NH2)-COOH
        # often participates in hydrogen bonding 
        # N Donor - NH2
        #
        # Essential AA
        # Molecular weight 146.19 Da
        # Ploar
        # Acidity - Basic
        # Hydrophobicity 0.283 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.9 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.90
        # pKa( alpha-COOH) 2.16
        # pKa( alpha-NH2) 9.06
        # CAS # 56-87-1
        # PubChem ID 866
        #
        self.Hydropathy = -3.9
        self.ResWeight = 90
        self.Hydropathy = -3.9
        self.ResWeight = 90
        self.name3L = 'LYS'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 1
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop self.loop = 0
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0,3:-1,5:-1}		# Special characteristic of residue  (change to {0:0,2:-1,4:-1} for LBH-15)
        self.n1 = -1			
        self.n2 = 0
        self.ResVol = 168.6
        self.SideChainVol = 168.6-54.1


class Arg(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'R')
        # Arginine
        # ###################
        # HN=C(NH2)-NH-(CH2)3-CH(NH2)-COOH
        # N Donor - NH2
        # 
        # Molecular weight 174.20 Da
        # Ploar
        # Acidity - Basic (strong)
        # Hydrophobicity 0.000 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -4.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point (when protonation accure) pH 10.76
        # pKa( alpha-COOH) 1.82
        # pKa( alpha-NH2) 8.99
        # CAS # 74-79-3
        # PubChem ID 5950
        #
        self.Hydropathy = -4.5 
        self.ResWeight = 118
        self.name3L = 'ARG'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 1
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop self.loop = 1
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = -1			
        self.n2 = -3
        self.ResVol = 173.4
        self.SideChainVol = 173.4-54.1

class His(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'H')
        # Histidine 
        # ################
        # NH-CH=N-CH=C-CH2-CH(NH2)-COOH
        # |__________|
        # N Donor - NH
        # The imidazole side chain has two nitrogens with different properties
        # 
        # Molecular weight 155.15 Da
        # Ploar
        # Acidity - Basic (week)
        # Hydrophobicity 0.165 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.2 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 7.60
        # pKa( alpha-COOH) 1.80
        # pKa( alpha-NH2) 9.33
        # CAS # 71-00-1
        # PubChem ID 773
        #
        self.Hydropathy = -3.2
        self.ResWeight = 99
        self.name3L = 'HIS'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0.5
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = -1			
        self.n2 = 0
        self.ResVol = 153.2
        self.SideChainVol = 153.2-54.1




def AminoTable():
    '''
    Dictionary of amino acids
    '''
    aminotable = {'A':Ala(),'V':Val(),'L':Leu(),'I':Ile(),'F':Phe(),'W':Trp(),'G':Gly(),'S':Ser(),
                  'T':Thr(),'C':Cys(),'Y':Tyr(),'N':Asn(),'M':Met(),'P':Pro(),
                  'D':Asp(),'E':Glu(),'Q':Gln(),'K':Lys(),'R':Arg(),'H':His()}
    return aminotable

def AminoList():
    AAlist = ['A','V','F','I','L','P','M','D','E','K','R',
              'S','T','Y','C','N','Q','H','W','G']
    return AAlist

def AA_name_convert():
    Namelist = {'ALA':'A','VAL':'V','PHE':'F','ILE':'I','LEU':'L','PRO':'P','MET':'M','ASP':'D',
                'GLU':'E','LYS':'K','ARG':'R','SER':'S','THR':'T','TYR':'Y','CYS':'C',
                'ASN':'N','GLN':'Q','HIS':'H','TRP':'W','GLY':'G'}

    return Namelist

def AminoAcidtoNum():
    '''
    maps amino acid to numbers. there is no significans to the order
    '''
    AAtoNum = {'A':1,'V':2,'L':3,'I':4,'F':5,'W':6,'G':7,'S':8,
               'T':9,'C':10,'Y':11,'N':12,'M':13,'P':14,
               'D':15,'E':16,'Q':17,'K':18,'R':19,'H':20}
    return AAtoNum

def NumtoAminoAcid():
    '''
    maps amino acid to numbers. there is no significans to the order
    '''
    NumtoAA = {1:'A',2:'V',3:'L',4:'I',5:'F',6:'W',7:'G',8:'S',
               9:'T',10:'C',11:'Y',12:'N',13:'M',14:'P',
               15:'D',16:'E',17:'Q',18:'K',19:'R',20:'H'}
    return NumtoAA



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
    if StructureName in ['LBH-18','LHBH-18','lbh-18','lhbh-18','lbh18','LBH18','LHBH18','lhbh18']:
        DataMatrix = [[1,-1,0,1,0,0,487.23,61.48,1,18,0],
                      [2,-1,3,3,0,0,487.23,61.48,2,18,36],
                      [3,1,0,0,63.83,30.09,487.23,61.48,3,18,0],	
                      [4,-1,0,0,0,0,487.23,61.48,4,18,0],
                      [5,1,0,0,99.14,21.22,487.23,61.48,5,18,15],
                      [6,-1,1,2,0,0,487.23,61.48,6,18,0]]
    elif StructureName in ['LBH-15','LHBH-15','lbh-15','lhbh-15','lbh15','LBH15','LHBH15','lhbh15']:
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
    cccounter = 0	# corner cut counter
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
    repeatedCornerCuts = 0


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
    # check for repeated corner cuts. add penalty for 3 in a row
    for t in range(len(seq1)):
        if seq1[t] == '|':
            if (seq1[t+1] == '-') or (seq1[t-1] == '-'):
                cccounter += 1
            else:
                cccounter = 0
            if cccounter > 2:
                repeatedCornerCuts += ScoresP.CornerCutPenalty * 10

            
    tot -= repeatedCornerCuts	# add repeated curner cut penalty
        

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
    print 'Repeated loop penalty score: ',-repeatedCornerCuts/nrows

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
            if iimpath.CornerMarker == 2: 
                Thread_data.thread += '|'	# No corner cutting
                Thread_data.cornerflag = 0
            else: 
                Thread_data.thread += '-|'				# corner cutting type 1
                Thread_data.cornerflag += 1
        elif iimpath.CornerMarker == 2: 
            Thread_data.thread += '|-'	# corner cutting type 2
            Thread_data.cornerflag += 1
        if loopFlag:		# Add loop 
            Thread_data.thread = Thread_data.thread + '(' + seq[(crow + 1):irow] + ')'
            Thread_data.loopFlag = maxLoopFlag
        else: Thread_data.loopFlag = ((Thread_data.loopFlag - 1) > 0)*(Thread_data.loopFlag - 1)
    Thread_data.thread += AA   

    return Thread_data

# end AddAAtoThread

#def ThreadsScoreTest(Thread_data,highScore):
    #'''
    #Delete all threads with less than highScore side chain bond
    #'''
    #NewThreaddataList = []
    #tmp = Thread_data.threadData
    #for i in tmp:
        #if i.score == highScore: NewThreaddataList += [i]
    ## NewThreaddataList = uniqueList(NewThreaddataList)
    #Thread_data.threadData = NewThreaddataList
    #return Thread_data

## End ThreadsScoreTest

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
        newScore -= (t.cornerflag > 1) * CornerCutFlag * ScoresP.CornerCutPenalty * 10	# add penalty for 3 corner cuts in a row
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

def PlotScoreMatrix (mpath,ProtName=''):
    '''
    Plotting the scoring matrix:  score vs. col at the last row
    '''
    top_row = nm.size(mpath,0)-2
    col = nm.size(mpath,1)
    yy = [x.withLoops.Score * (x.withLoops.Score > 0) for x in mpath[top_row]]
    y = nm.array(yy)
    x = nm.array(range(col))

    # Plotting
    TitleString = 'Matrix Scores : ' + ProtName
    xLabelStr = 'Scoring Matrix Column'
    yLabelStr = 'Scores per amino acid'

    pl.title(TitleString)
    pl.xlabel(xLabelStr)
    pl.ylabel(yLabelStr)
    pl.plot(x,y)
    pl.ylim([0,110])
    pl.show()

# end PlotScoreMatrix 


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

def ThreadSeq(seq,ScoreTolerance=0,StructureType='LBH-18',ProteinName = ''):
    '''
    This function is identical to the betatesting just in a function form
    seq : the sequence to be threaded
    seq should be in this format -  'LCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQGTHSQWNKPSKPKTNMKHMAGAA'

    ProteinName is an optional string that will be put in the scan plot title
    '''
    # convert to upper case
    seq = seq.upper()
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
    # Current availble structures are LBH-18 and LBH-15
    Structure_rules = StructureRulls(StructureType)
    print 'Finish reading structure rules \n'
    # End structure rules reading

    # read some program restrictions
    print 'Enter data'
    # Read the structure information
    BaseData = readProgParameters(Structure_rules)
    # Set the tolerance for scanning
    BaseData.tolerance = ScoreTolerance
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

    PlotScoreMatrix (mpath,ProteinName)
    # Return nothiong, just print results
    return

# End ThreadSeq

def ThreadSeq15(seq,ScoreTolerance=0,ProteinName = ''):
    '''
    This function is identical to the betatesting just in a function form
    seq : the sequence to be threaded
    seq should be in this format -  'LCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQGTHSQWNKPSKPKTNMKHMAGAA'

    ProteinName is an optional string that will be put in the scan plot title
    '''
    # convert to upper case
    seq = seq.upper()
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
    # Current availble structures are LBH-15
    Structure_rules = StructureRulls('LBH-15')
    print 'Finish reading structure rules \n'
    # End structure rules reading

    # read some program restrictions
    print 'Enter data'
    # Read the structure information
    BaseData = readProgParameters(Structure_rules)
    # Set the tolerance for scanning
    BaseData.tolerance = ScoreTolerance
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

    PlotScoreMatrix (mpath,ProteinName)
    # Return nothiong, just print results
    return

# End ThreadSeq15

def AnalyzeProteinSeq(segLength,LoopLength,seq,StructureSTR='LBH-18'):
    '''
    You need to be in the same folder as the protein file to run this program   
    This program is for producing protein scan plot only   
    Valid vaues for StructureSTR = 'LBH-18'  or   'LBH-15

    segLength,LoopLength or two integers representing the seqments length used in testing 
    and the maximum loop length allowed
    '''
    # convert to upper case
    seq = seq.upper()

    Structure_rules = StructureRulls(StructureSTR)
    ScoresP = ScoreParameters()		# Scoring parameters


    # Reading the protien sequence
    #OutFileName = 'plotdata'
    #WinFilePath = 'C:\Users\youval\Documents\Education\Research\Paper 2 - Threading program\Plots'
    OutFileName = 'test2new3'
    WinFilePath = 'C:\Users\youval\Downloads\\test2'
    file_path = convertPathNames(WinFilePath)
    # file_name = 'proteinSequence.txt'
    OutputFile = file_path + OutFileName +'.txt'
    print seq
    # Extacting and calculating basic test info
    Hi_score = 60				# Score above which I'll pressent threads
    tol = 3						# The range of score from the run top score for which treads will be collected

    # Creating a file with plot data points
    #open file for sequence writing
    f = open(OutputFile, 'w')
    x,y,ScoreThreadList = Seq_file_score_LBH18(seq,Structure_rules,ScoresP,1,LoopLength,segLength,17,tol,Hi_score)

    # Disply and save results
    print '================================'
    print 'Results :  Position in the protein sequence  ,  Score'
    print '================================'

    for i in range(len(x)):
        f.write(str(x[i]))
        f.write(' ')
        f.write(str(y[i]))
        f.write('\n')
        print x[i],y[i]
    f.close()    

# End AnalyzeProteinSeq

def ThreadScores(ThreadSTR,StructureSTR='LBH-18'):
    '''
    this function calculates a thread score and the weights of the different parts of the score

    StructureSTR = 'LBH-18'  or   'LBH-15

    '''
    # Convert to upper case
    ThreadSTR = ThreadSTR.upper()
    # Get basic data
    aminotable = AminoTable()	# Dictionary of amino acids
    AAlist = AminoList()
    Structure_rules = StructureRulls(StructureSTR)
    # read some program ristrictions
    BaseData = BaseInputData(Structure_rules)
    build_Reverse_in_sidechain_list(BaseData,Structure_rules)
    # Scoring parameters ratios
    ScoresP = ScoreParameters()

    # Collect clean amino acid sequence and combine to 1 string
    seqA = [x for x in ThreadSTR if x in AAlist]	
    seqA = st.join(seqA,'')
    # Get the correct thread score
    FixScoreA = GetThreadScores(seqA,ThreadSTR,Structure_rules,ScoresP,BaseData)

    # Printing info
    print '============================'
    print 'Total score results'
    print '%.1f :: %s' % (FixScoreA,ThreadSTR)

# End ThreadScores

def readProgParameters(Structure_rules):
    BaseData = BaseInputData(Structure_rules)
    # read shortest allowed loop length
    BaseData.shortl = read_var(1,1,'Enter shortest allowed loop length')
    # read longest allowed loop length  
    BaseData.longl = read_var(4,1,'Enter longest allowed loop length')
    # What range of score we would like to get information on   
    BaseData.dScore = read_var(3.0,2,'Enter results score range')   
    # Tolerance in scoring - scores that are within the tolerance will be consider to be a good path score   
    #BaseData.tolerance = read_var(0.0,2,'Enter scoring tolerance')      
    # minimum steps needed before next loop is allowed   
    BaseData.maxLoopFlag = read_var(17,1,'Enter minimum number of steps between loops')  
    # minimum steps needed before next corner cut is allowed   

    return BaseData
# End of readProgParameters


def read_var(default_val,inpType,qStr):
    '''
    default_val is the defualt value of the input
    inpType 1:integer, 2:float, 3:string
    qStr is the question
    '''
    Qstr = qStr + ' [%s]: ' % str(default_val)
    inpStr = raw_input(Qstr)
    if inpStr == '':
        outVal =default_val
    elif inpType == 1:
        outVal = int(inpStr)
    elif inpType == 2:
        outVal = float(inpStr)
    elif inpType == 3:
        outVal = inpStr
    else:
        print 'Wrong input type'
    return outVal
# End of reading function
