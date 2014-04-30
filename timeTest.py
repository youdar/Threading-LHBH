import numpy as nm
import string as st
from TheaderCommonFunctions import *
from ReadDate import *
from ThreadClasses import *
from AminoData import *
import pylab as pl
import copy as cp
import time

aminotable = AminoTable()	
AAlist = AminoList()
Structure_rules,s_file_path,s_file_name = structure_file_info()
BaseData = BaseInputData(Structure_rules)
ScoresP = ScoreParameters()

def ThreadScores(seq,seq1,Structure_rules,ScoresP,BaseData,AAlist,ThreadMarks,aminotable,mpath):
    '''
    Returns the score of a know thread seq1
    Scores : the position score (as a class variable)
    tot : the total score
    '''
    
    icol = StartCol(seq1,base_length) 	# The know structure thread starting column
    irow = 0
    NoLoopFlag = True
    PScore = 0.0		# The side chain score
    i = 0			# Counts the steps in the seq1
    tot = 0.0

    for t in seq1:
        if t== '(': NoLoopFlag = False	# we've riched a loop
        if NoLoopFlag:
            if t not in ThreadMarks:
                AA = aminotable[t]
                Mpos = mpath[irow][icol]
                if t == 'C': turnLenght = mpath[irow][icol].DiSulDis
                else: turnLenght = mpath[irow][icol].Turn
                PScore += threadSidescore(t,seq1[0:i],turnLenght,ScoresP)	# Side chain bond score
                PScore += ThreadVolumeScore(seq1[0:i],BaseData)*ScoresP.TurnVolPenaly
                tmpScr,total_scr = PosScore(AA,Mpos,ScoresP)			# position score
                tot += total_scr
                icol += 1            
        if not NoLoopFlag:	# we are in a loop
            # Score for loops
            if t in AAlist: 
                AA = aminotable[t]
                tot += AA.loop * ScoresP.SpecialRatio	# Adding special score for selected AA
                tot -= AA.Hydrophobic * ScoresP.HydroRatio * ScoresP.LoopWeight	# adding the hydrophobicity score of AA on loop


        if t in AAlist: irow += 1        
        if t == '-': 					# Corner cut
            icol += 1					# Continue counting
            tot -= ScoresP.CornerCutPenalty		# panelty for corner cutting

        if t == ')': 
            NoLoopFlag = True		# End of loop
            tot -= ScoresP.LoopPenalty / 5.0		# Small panelty for loops
            tot -= (mpath[irow][icol].CornerMarker in [0,2]) * ScoresP.LoopPenalty 	# apply loop penalty if not in corner
        i += 1		

    tot += PScore		# the total score
    tot = tot / float(row)	# Normelizing the to the number of amino acid
    return tot


start_time = time.clock()
AAlist = AminoList()
ThreadMarks = ['|','-','(',')']
aminotable = AminoTable()
# the size of our scoring matrix
row = 72			
# Add extra colomns to allow corner cutting
base_length = len(Structure_rules)
# Testing how many corners can be cut
col = end_col(row,Structure_rules)
# initializing a list of lists that will store the best routes in mpath matrix
mpath = build_path_matrix_v3(row,col,Structure_rules)

for i in range(1000):  
    seq1 = 'AQVHIG|NFVEVK|-GSSIG|ENTKAG|HLTYIG|-NCEVG|SNVNFG|AGTITV|(NYDGKNK)YKTVIG|NNVFVG|SNSTII'
    seq = [x for x in seq1 if x in AAlist]	
    seq = st.join(seq,'')
    FixScoreA = ThreadScores(seq,seq1,Structure_rules,ScoresP,BaseData,AAlist,ThreadMarks,aminotable,mpath)
dtime = (time.clock() - start_time)
print 'Calculation time is %.2f sec' % dtime
print '%.1f :: %s' % (FixScoreA,seq1)

