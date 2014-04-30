# ScoreQuality_minuit.py is a program that reads protein sequence and its correct confermation
# reads structure information, calculates the scores of the scoring matrix using different scoring
# parameters ratios
# The program uses pyminuit to optimize the calculated score in a way that it will match the correct 
# confermation score


# Import modules
import numpy as nm
import string as st
from TheaderCommonFunctions import *
from ReadDate import *
from ThreadClasses import *
from AminoData import *
import pylab as pl
import copy as cp
import minuit

def GetScoreQuality(x2,x3,x4,x8,x9,x10,x11,x12):   
    # This is the function that needs to be minimized
    # it calculates the difference in the calculated score and the score
    # of the correct confermation
    #
    # we vary only 4 parameters and hold all the rest fixed

    #ScoresP.HydroRatio = x1
    ScoresP.PolarRatio = x2
    ScoresP.ChargeRatio = x3
    ScoresP.SizeRatio = x4
    #ScoresP.DiS_Hyd_ration = x5
    #ScoresP.SpecialRatio = x6
    #ScoresP.CornerRatio = x7		# not in use
    ScoresP.SideBondRatio = x8
    ScoresP.LoopPenalty = x9		# penalty for loops when not at the corner
    ScoresP.CornerCutPenalty = x10
    ScoresP.TurnVolPenaly = x11
    ScoresP.LoopWeight = x12

    # Get the score matrix 
    mpath,row,col,highScore = mpathCalc(Structure_rules,seq,BaseData,ScoresP)
    # Get the correct thread score
    FixScore = GetThreadScores(k_seq,k_seq1,Structure_rules,ScoresP,BaseData)
    # Collect score matrix scores
    yy = nm.array([x.withLoops.Score for x in mpath[row-1]])      
    ymax = max(yy)  
    f = abs((ymax - FixScore)/ymax) 
    print 'seq1',seq1
    print 'seq',seq
    print 'ymax',ymax
    print 'k_seq1',k_seq1
    print 'k_seq',k_seq
    print 'fixscore',FixScore
    return f

aminotable = AminoTable()	# Dictionary of amino acids


# Read structure rules
structur_type = read_var(2,1,'Collect data on (1)LBH15, (2)LBH18')
if structur_type not in [1,2]: structur_type = 2
if structur_type == 1:
    ttl = 'Data for LHB15'
    s_file_name = 'StructureFile-LeftBetaH-15.txt'
else:
    ttl = 'Data for LHB18'
    s_file_name = 'StructureFile-LeftBetaH-18.txt'
s_file_path = 'C://Users//youval//Documents//Education//work biophy//Structure and scoring rules//'    
Structure_rules = read_structure_rules(s_file_name, s_file_path)
# read some program ristrictions
BaseData = readProgParameters(Structure_rules)
build_Reverse_in_sidechain_list(BaseData,Structure_rules)
# Read the protein sequence
# You need to put in files test_pdb and test, the correct structure and the sequence
file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//test//'
file_name = 'test'
seq,seq1  = read_known_seq_file(file_name, file_path)
# End of reading the seq

# read known sequence information
file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//test//'
file_name = 'test_pdb'
k_seq,k_seq1  = read_known_seq_file(file_name, file_path)
print '--------------------------'
print seq
print '\n'
print k_seq1
print '--------------------------'
# scoring parameters ratios 
# setting initial values
ScoresP = ScoreParameters()
xx1 = ScoresP.HydroRatio 
xx2 = ScoresP.PolarRatio
xx3 = ScoresP.ChargeRatio 
xx4 = ScoresP.SizeRatio 
#x5 = ScoresP.DiS_Hyd_ration - Held fixed at default
#x6 = ScoresP.SpecialRatio - Held fixed at default
#x7 = ScoresP.CornerRatio - Not in use
xx8 = ScoresP.SideBondRatio 
xx9 = ScoresP.LoopPenalty 
xx10 = ScoresP.CornerCutPenalty 
xx11 = ScoresP.TurnVolPenaly
xx12 = ScoresP.LoopWeight

print 'start looking for optimal parameters'
# limit_x$ set the value limits for the parameter x$
# x$ = xx$ sets the starting point for x$
m = minuit.Minuit(GetScoreQuality,limit_x2=(20,90),
                  limit_x3=(50,150),limit_x4=(50,350),limit_x8=(10,100),limit_x9=(10,650),
                  limit_x10=(20,250),limit_x11=(10,1500),limit_x12=(0.3,0.9),
                  x2=xx2,x3=xx3,x4=xx4,x8=xx8,x9=xx9,x10=xx10,x11=xx11,x12=xx12,
                  err_x2=10,err_x3=10,err_x4=10,
                  err_x8=10,err_x9=10,err_x10=10,err_x11=10,err_x12=10)
# Strategy 0,1,2 is the level of accuracy
m.strategy = 1    #  	 0 = fast, 1 = default, 2 = thorough 
# tol is the tolerance we ask for in the values we look for

#print '\nm.limits',m.limits
#print '\nm.parameters',m.parameters
#print 'm.values',m.values
#print 'm.strategy',m.strategy
#print 'm.edm',m.edm
print '\nx1 = %.2f ScoresP.HydroRatio - Held fixed at default' % xx1
print 'x2 = %.2f ScoresP.PolarRatio' % xx2
print 'x3 = %.2f ScoresP.ChargeRatio' % xx3
print 'x4 = %.2f ScoresP.SizeRatio' % xx4
print 'x5 = ScoresP.DiS_Hyd_ration - Held fixed at default' 
print 'x6 = ScoresP.SpecialRatio - Held fixed at default' 
print 'x7 = ScoresP.CornerRatio - Held fixed at default' 
print 'x8 = %.2f ScoresP.SideBondRatio' % xx8
print 'x9 = %.2f ScoresP.LoopPenalty' % xx9
print 'x10 = %.2f ScoresP.CornerCutPenalty ' % xx10
print 'x11 = %.2f ScoresP.TurnVolPenaly' % xx11
print 'x12 = %.2f ScoresP.LoopWeight \n' % xx12


do_something = raw_input('Press any key to continue')

m.printMode = 1 # 0 = print nothing, 1 = print call-by-call parameter values, 
                # 2 = differences from starting point, 3 = differences from previous step
m.migrad() 	# MIGRAD algorithm, an optimized gradient-based minimum search. 
                # This function's behavior is affected by 
m.printMode = 0

print 'done with optimization'
