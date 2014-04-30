'''
LHBproteinSearch.py is a program that contains search functions that
search a pdb file for regions of possible left-hand-beta-helix structure
'''

# Import modules
import numpy as nm
import string as st
import pylab as pl
import copy as cp
import time 
from TheaderCommonFunctions import *
from ReadDate import *
from ThreadClasses import *
from AminoData import *


def SearchLHB18(pbdFileName,shortl=1,longl=2,pseq_length=42,minLoopDis=16):
    '''
    This function scans a pdb file for possible LHB-18 regions
    It outputs a list where all LHB marked by 1 and other structures marked by 0
    It output the list of amino acids in the LHB regions and the possible threads

    pbdFileName : the name of the pdb file
    shortl=1,longl=2,minLoopDis=16 : default parameters for scoring
    pseq_length is length of the of sub sequence on which we run the scoring algorithm
    '''
    start_time = time.clock()
    print 'Searching for possible LHB18 structure' # Making sure we running LHB-18 search


    # General amino acid data, tables
    #AAlist = AminoList()
    #aminotable = AminoTable()
    # Structure rule file name and location
    s_file_name = 'StructureFile-LeftBetaH-18.txt'
    s_file_path = 'C://Users//youval//Documents//Education//work biophy//Structure and scoring rules//'       
    Structure_rules = read_structure_rules(s_file_name, s_file_path)
    # init
    BaseData = BaseInputData(Structure_rules)
    BaseData.shortl = shortl
    BaseData.longl = longl
    BaseData.maxLoopFlag = minLoopDis
    # Scoring parameters
    ScoresP = ScoreParameters()
    # reading protein sequence
    #pbdFileName += '.pdb'
    #file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//pdb files//'
    #file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//Diseas related proteins//p53//'
    file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//test//'
    file_name = 'testLBH'
    seq  = read_seq_file(file_name, file_path)	# Get amino acid sequence
    #seq = Get_pdb_seq_info(file_name, file_path)	# Get amino acid from pdb file
    seq = seq.upper()
    # parce the sequence

    side_length = 6				# side length in the LHB structure model
    n_iter = (len(seq)- pseq_length)/side_length + 1	# number of iterations
    scores_array = nm.zeros(n_iter)		# Array for results collection
    x = nm.array(range(n_iter))
    for i in x:
        ii = i * side_length
        p_seq = seq[ii:(ii+pseq_length)]	# sending two turns for scoring
        score = get_seq_score(Structure_rules,p_seq,BaseData,ScoresP)
        scores_array[i] = score			# Collect the scores

    # Calculation time
    dtime = (time.clock() - start_time)/60
    print 'Calculation time is %.1f min' % dtime
    # Evaluating LHB regions
    x2 = nm.array([])
    for i in scores_array:
        if i>70:x2 = nm.append(x2,100)
        else:x2 = nm.append(x2,10)
    # Plotting
    x *= side_length		# addjusting the valuse of x to represent amino acid number
    TitleString = 'Scores for %s - steps of %i amino acids' % (pbdFileName,side_length)
    xLabelStr = 'Position - amino acid start + %i, loops (%i : %i) ' % (pseq_length,shortl,longl)
    yLabelStr = 'Scores' 
    #pl.subplot(311)
    pl.figure(1)
    pl.title(TitleString)
    pl.xlabel(xLabelStr)
    pl.ylabel(yLabelStr)
    pl.plot(x,scores_array,'r--o')  	# r fro red, --o both dash line and dots
    pl.ylim((0,120))
    pl.xlim((0,(len(seq)+10)))
    #pl.figure(2)
    #pl.subplot(313)
    #TitleString = 'Estimated possible LBH regions'
    #pl.title(TitleString)
    #pl.xlabel(xLabelStr)
    #pl.plot(x,x2,'b')
    #pl.ylim((0,120))
    #pl.xlim((0,(len(seq)+10)))
    pl.show()

def SearchLHB15(pbdFileName,shortl=1,longl=2,pseq_length=42,minLoopDis=16):
    '''
    This function scans a pdb file for possible LHB-18 regions
    It outputs a list where all LHB marked by 1 and other structures marked by 0
    It output the list of amino acids in the LHB regions and the possible threads

    pbdFileName : the name of the pdb file
    shortl=1,longl=2,minLoopDis=16 : default parameters for scoring
    pseq_length is length of the of sub sequence on which we run the scoring algorithm
    '''
    start_time = time.clock()
    print 'Searching for possible LHB15 structure' # Making sure we running LHB-15 search

    # Structure rule file name and location
    s_file_name = 'StructureFile-LeftBetaH-15.txt'
    s_file_path = 'C://Users//youval//Documents//Education//work biophy//Structure and scoring rules//'       
    Structure_rules = read_structure_rules(s_file_name, s_file_path)
    # init
    BaseData = BaseInputData(Structure_rules)
    BaseData.shortl = shortl
    BaseData.longl = longl
    BaseData.maxLoopFlag = minLoopDis
    # Scoring parameters
    ScoresP = ScoreParameters()
    # reading protein sequence
    #pbdFileName += '.pdb'
    #file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//pdb files//'
    #file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//Diseas related proteins//p53//'
    file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//test//'
    file_name = 'testLBH'
    seq  = read_seq_file(file_name, file_path)	# Get amino acid sequence
    #seq = Get_pdb_seq_info(pbdFileName, file_path)	# Get amino acid from pdb file
    seq = seq.upper()
    # parce the sequence
    
    side_length = 5				# side length in the LHB structure model
    n_iter = (len(seq)- pseq_length)/side_length + 1	# number of iterations
    scores_array = nm.zeros(n_iter)		# Array for results collection
    x = nm.array(range(n_iter))
    for i in x:
        ii = i * side_length
        p_seq = seq[ii:(ii+pseq_length)]	# sending two turns for scoring
        score = get_seq_score(Structure_rules,p_seq,BaseData,ScoresP)
        scores_array[i] = score			# Collect the scores

    # Calculation time
    dtime = (time.clock() - start_time)/60
    print 'Calculation time is %.1f min' % dtime
    # Evaluating LHB regions
    x2 = nm.array([])
    for i in scores_array:
        if i>70:x2 = nm.append(x2,100)
        else:x2 = nm.append(x2,10)
    # Plotting
    x *= side_length		# addjusting the valuse of x to represent amino acid number
    TitleString = 'Scores for %s - steps of %i amino acids' % (pbdFileName,side_length)
    xLabelStr = 'Position - amino acid start + %i, loops (%i : %i) ' % (pseq_length,shortl,longl)
    yLabelStr = 'Scores' 
    #pl.subplot(311)
    pl.figure(1)
    pl.title(TitleString)
    pl.xlabel(xLabelStr)
    pl.ylabel(yLabelStr)
    pl.plot(x,scores_array,'r--o')  	# r fro red, --o both dash line and dots
    pl.ylim((0,90))
    pl.xlim((0,(len(seq)+10)))
    #pl.figure(2)
    #pl.subplot(313)
    #TitleString = 'Estimated possible LBH regions'
    #pl.title(TitleString)
    #pl.xlabel(xLabelStr)
    #pl.plot(x,x2,'b')
    #pl.ylim((0,120))
    #pl.xlim((0,(len(seq)+10)))
    pl.show()
    
    
SearchLHB18('Prion',1,2,7*6,17)
#SearchLHB15('Prion',1,2,12*5,14)
