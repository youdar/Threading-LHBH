'''
This program use known protein sequences to generate statistical 
data, that can use to determine scores for amino acid that points
into the structure

I gathered all known LHB helix threads in one file. only the helix portion 
of the protein. And use those threads to extract statistical data.

probability for a residue to be in a specific position
probability for a turn to have total number of atom pointing in
'''

import numpy as nm
import string as st
from TheaderCommonFunctions import *
from ReadDate import *
from ThreadClasses import *
from AminoData import *
import pylab as pl
import copy as cp


aminotable = AminoTable()				# Dictionary of amino acids
AAlist = AminoList()
FullAAList = AAlist + ['-']
ThreadMarks = ['|','-','(',')']
AAtoNum = AminoAcidtoNum()	# maps amino acid to numbers. there is no significans to the order
NumtoAA = NumtoAminoAcid()	# maps amino acid to numbers. there is no significans to the order
CharegdeAminoList  = ['D','E','K','R']
PosCharegdeAminoList  = ['K','R']
NegCharegdeAminoList  = ['D','E']
PolarAminoList = ['S','T','Y','H','C','N','Q','W']
HPobicAminoList = ['A','V','F','P','M','I','L']
ThreadMarks = ['|','-','(',')']


# Read the protein sequence
file_path = 'C://Documents and Settings//youval//My Documents//Education//work biophy//Protein seq//Collection of Known LHB//'
s_file_path = 'C://Users//youval//Documents//Education//work biophy//Structure and scoring rules//' 
# List of known threads LBH 15
Known_thread_files_LHB15 = ['1L0S_pdb','1m8n_pdb','1m8n-chain-B_pdb','1N4I_pdb']
# for testing ...
#Known_thread_files_LHB15 = ['test_pdb']
# List of known threads LBH 18
Known_thread_files_LHB18 = ['1G97_pdb','1HV9_pdb','1KRR_pdb','1LXA_pdb', '1OCX_pdb','1QRE_pdb','1SSQ_pdb',
                            '1TDT_pdb', '1V3W_pdb','1XAT_pdb','1XHD_pdb','3TDT_pdb']
# '1QRE_pdb' and '1QQ0_pdb' are the same
                                                     

#What are we working on?
structur_type = read_var(2,1,'Collect data on (1)LBH15, (2)LBH18')
if structur_type not in [1,2]: structur_type = 1
if structur_type == 1:
    ttl = 'Data for LHB15'
    Known_thread_files_LHB = Known_thread_files_LHB15
    s_file_name = 'StructureFile-LeftBetaH-15.txt'
else:
    ttl = 'Data for LHB18'
    Known_thread_files_LHB = Known_thread_files_LHB18
    s_file_name = 'StructureFile-LeftBetaH-18.txt'


Structure_rules = read_structure_rules(s_file_name, s_file_path)
base_length = len(Structure_rules)
turn_length = Structure_rules[0].Turn


# Collect the number of each AA in each position
AA_data_collect = nm.zeros([21,(base_length+1)])	
# residue position for reisdues that point in
ResPointIn = [x.PosNum for x in Structure_rules if x.PointInOut == 1] 
# colelct the residues weight at position that points into structure
ResWeightData = nm.zeros((base_length+1))	
# convert to a list
ResWeightData = [[x] for x in ResWeightData]
# counts the total number of amino acid at each position
ResPosCounter = nm.zeros((base_length+1))
# Weigth of residues in a turn
TurnResWeight = nm.zeros(2000)

iside = 0		# counts the total number of sides
print 'ResPointIn',ResPointIn
for file_name in Known_thread_files_LHB:
    seq,seq1  = read_known_seq_file(file_name, file_path)	# Read the protein sequence
    icol = StartCol(seq1,base_length)+1	# get the position of the first amino acid (Starts at 1)
    tcol = 0		# amino acid in a turn counter
    iturn = 0		# counts the turn
    turnflag = False	# will turn to true when we see the first corner
    turnWeight = 0	# collects the residues weight in a turn
    NotInLoop = True	# indicates to stop counting the columns when in a loop
    CountResFlag = False	# Indicate if we should count the residues pointing in
    ColFlag = False	# indication if col<2



    for AA in seq1:     
        if AA in FullAAList:			# FullAAList is AA list + '-'
            col = icol%base_length		# get the position number in the structure (1-6 in LBH18)
            col += base_length*(col == 0)	# if icol%base_length = 0 add base_length     
            if (AA in AAlist) and NotInLoop: 
                row = AAtoNum[AA]	# find the row in the data collection matrix, each row num represents amnio acid
                if (not turnflag) and (col<3): ColFlag = True	# Start collecting data before getting to first turn
                if col in ResPointIn: 
                    # ResWeightData[col] +=  [aminotable[AA].ResWeight] 
                    ResWeightData[col] +=  [aminotable[AA].SideChainVol] 	# Collecting total volume of residue pointing in
                    CountResFlag = turnflag or ColFlag
                    if CountResFlag:
                        # turnWeight +=  aminotable[AA].ResWeight
                        turnWeight +=  aminotable[AA].SideChainVol	# Collecting the volume of one side in the LBH
                if (col == base_length) and CountResFlag: 
                    TurnResWeight[iside] = turnWeight		# Storing side volume
                    turnWeight = 0
                    iside += 1					# Countiung number of sides
            else: row = 0 # note that all AA that are in a loop will be counted at row 0

            AA_data_collect[row,col] += 1	# add the amino acid data in the data collection matrix
            #if turnflag and NotInLoop:
                #tcol += 1
                #tcol = tcol%turn_length
                #tcol += turn_length*(tcol == 0)
                #if (AA in AAlist): turnWeight +=  aminotable[AA].ResWeight
                #if (tcol == turn_length):
                    #TurnResWeight[iturn] = turnWeight
                    #turnWeight = 0
                    #iturn += 1


        # Setting up loop status            
        if AA == '(': NotInLoop = False
        if NotInLoop: icol += 1   
        if AA == ')': NotInLoop = True
        if AA in ['|']: 
            icol -= 1
            turnflag = True

# Deleting row 0 and col 0
AA_data_collect= nm.delete(AA_data_collect,0,0)
AA_data_collect= nm.delete(AA_data_collect,0,1)
del(ResWeightData[0])
ResPosCounter = nm.sum(AA_data_collect,axis=0)		# Total Residues in each position      
TurnResWeight = nm.array([x for x in TurnResWeight if x !=0])	# Collecting only the non zero side residue weights
#nturns = float(len(TurnResWeight) - 2)		# number of data points 
#nturns = len(TurnResWeight)%3 + 3*(len(TurnResWeight)/3)	# number of data points
#TurnMeanWieght = nm.zeros(nturns)			# will be used to creat mean trun wieght
if iside != len(TurnResWeight):
    print 'Check side counter!!!   iside not equal TurnResWeight length'
nsides = float(iside)		# number of sides
print '\n---------------------------------------------'
print 'nsides',nsides
print 'iside',iside
print 'TurnResWeight - Side volume'
print TurnResWeight
print '---------------------------------------------\n'

#TurnResWeight = nm.array([x for x in TurnResWeight])
TurnResWeight = Convert_to_turn_volume(TurnResWeight)
print 'TurnResWeight - Full turn volume'
print TurnResWeight
print '#########################\n'


# Statistics of turn weight	
#MeanTurnWeight = sum(TurnMeanWieght)/nturns
nturns = float(len(TurnResWeight))	# numbers of full turns volume    
MeanTurnWeight = sum(TurnResWeight)/nturns
TurnWeightVariance = sum([(x - MeanTurnWeight)**2 for x in TurnResWeight])/(nturns-1.)
TurnWeightSigma = TurnWeightVariance**.5

print ttl	# Print the type of structure we are getting the information for
for i in range(len(AA_data_collect)): 	# print information for each amino acid   
    print NumtoAA[i+1],AA_data_collect[i],int(nm.sum(AA_data_collect[i]))

print '\n*** Results for LHB ***'    
print 'Statistical info on residue volume, for residues that are pointing into the structure \n'
for i in range(len(ResWeightData)):
    d = cp.copy(ResWeightData[i])
    d = [x for x in d if x > 0]		# not counting any zeros
    if len(d) > 1:
        n = float(len(d))		# number of data points
        dataMean = sum(d)/n
        dataVariance = sum([(x - dataMean)**2 for x in d])/(n-1.)
        Sigma = dataVariance**.5
        print '-------------------------------------'
        print 'For position %g' % (i+1)
        print 'dataMean : %.2f' % dataMean
        #print 'dataVariance : %.2f' % dataVariance
        print 'standard deviation : %.2f \n' % Sigma

print '----------------------------------------------'
print 'Mean Turn Volume %.2f' % MeanTurnWeight
print 'Turn Volume Sigma %.2f' % TurnWeightSigma
