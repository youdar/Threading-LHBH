'''
    betatesting_v#.py: Reads a protein sequence, reads structure energy and stability rules.
    The progrm use dynamic programming method to thread the sequence to a left-hand beta and 
    returns threds with the best score 
'''

# Import modules
#from math import *
#from sys import *
import numpy as nm
import string as st
from TheaderCommonFunctions import *
from ReadDate import *
from ThreadClasses import *
from AminoData import *
import copy as cp
import time


aminotable = AminoTable()	# Dictionary of amino acids

# scoring parameters ratios
ScoresP = ScoreParameters()

# Read the protein sequence
seq,file_path,file_name = seq_file_info()
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
    # saving the result to a file
    result_output_v4(file_name,file_path,seq,highscore,threadlist1,BaseData.dScore)
    # keeping the terminal open
else:
    print 'No Threads, Program Error'

ttl = file_name[0:5]
PlotScoreMatrix (mpath,ttl)

raw_input()

    
