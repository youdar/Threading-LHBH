# ScorePlot.py is a program that reads protein sequence and its correct confermation
# reads structure information, creates a scoring matrix, calculates the score of the
# confermation and then plots the scores vs. the column in matrix and the score of the 
# correct confermation.


# Import modules
import numpy as nm
import string as st
from TheaderCommonFunctions import *
from ReadDate import *
from ThreadClasses import *
from AminoData import *
import pylab as pl
import copy as cp

aminotable = AminoTable()	# Dictionary of amino acids

# Read the protein sequence
seq,file_path,file_name = seq_file_info()
# End of reading the seq

# Read structure rules
Structure_rules,s_file_path,s_file_name = structure_file_info()
# End structure rules reading

# read some program ristrictions
BaseData = readProgParameters()

# scoring parameters ratios
ScoresP = ScoreParameters()


# read known sequence information
k_seq,k_seq1,k_file_path,k_file_name = known_seq_info()


# Printing info
prtStr1 = '\nComparing the files: %s, with %g amino acids' % (file_name,len(seq))
prtStr2 = 'with the files: %s, that have %g amino acids \n' % (k_file_name,len(k_seq))
print prtStr1
print prtStr2
print seq
print k_seq,'\n'
print 'The correct thread'
print k_seq1

# Get the score matrix 		
mpath,row,col = mpathCalc(Structure_rules,seq,BaseData,ScoresP)

# Get the correct thread score
Scores,FixScore = GetThreadScores(k_seq,k_seq1,Structure_rules,ScoresP)

yy = nm.array([x.withLoops.Score for x in mpath[row-1]])
imin = 0
imax = col
while yy[imin] == 0: imin += 1
while yy[imax-1] == 0: imax -= 1
y = yy[imin:imax]   
topMatrixScore = max(y)
lowMatrixScore = min(y)
x = nm.arange(imin,imax)
y2 = nm.array(len(x) * [FixScore])
HighY = max(topMatrixScore,FixScore)
LowY = min(lowMatrixScore,FixScore)

# test position on plot
TextY = (HighY + LowY)*.35

# the 'f' function that indicates how good is the scoring process
f = abs(topMatrixScore - FixScore)/topMatrixScore

# Plotting
TitleString = 'Matrix Scores and the Score of the correct thread'
textStr1 = 'Highest matrix score %.1f \n' % topMatrixScore
textStr2 = 'Correct thread score %.1f \n' % FixScore
textStr3 = 'f = %.2f (Good score is f -> 0)' % f
textStr = textStr1 + textStr2 + textStr3
xLabelStr = 'Scoring Matrix Column'
yLabelStr = 'Scores'

pl.title(TitleString)
pl.xlabel(xLabelStr)
pl.ylabel(yLabelStr)
pl.plot(x,y,x,y2)

pl.text(imin,TextY,textStr)
pl.show()
