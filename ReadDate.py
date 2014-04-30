# temporary list of file used in the threading programs

#from sys import *
import numpy as nm
import string as st
from AminoData import *
from TheaderCommonFunctions import *

def seq_file_info():
    # Read the protein sequence
    aminotable = AminoTable()	# Dictionary of amino acids
    file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//'
    #file_name = get_file_name(file_path)
    file_name = 'test//test'
    #file_name = 'test2'
    #file_name = 'test3'
    #file_name = 'protseq1.txt'
    #file_name = 'protseq2.txt'
    #file_name = '1m8n//1m8n'
    #file_name = '1TDT'
    #file_name = '1LXA'
    #file_name = '1m8n-part1'
    #file_name = '247-286.txt'
    #file_name = '248-286.txt'
    #file_name = 'SUP35_scap.txt'
    #file_name = 'ure2P_scap.txt'
    seq  = read_seq_file(file_name, file_path)
    return seq,file_path,file_name

def structure_file_info():
    #s_file_name = 'StructureFile-LeftBetaH-18.txt'
    s_file_name = 'StructureFile-LeftBetaH-15.txt'
    s_file_path = 'C://Users//youval//Documents//Education//work biophy//Structure and scoring rules//'    
    Structure_rules = read_structure_rules(s_file_name, s_file_path)
    return Structure_rules,s_file_path,s_file_name

def known_seq_info():
    # Read the protein sequence
    file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//'
    #file_name = '1g97//1g97_pdb'
    #file_name = '1m8n_pdb'
    #file_name = '1TDT_pdb'
    #file_name = '1m8n-part1_pdb'
    file_name = 'test//test_pdb'
    seq,seq1  = read_known_seq_file(file_name, file_path)
    return seq,seq1,file_path,file_name

def read_known_seq_file(file_name, file_path):
    '''
    Read an amino acid from a file
    '''
    AAlist = AminoList()
    aminotable = AminoTable()
    f = open(file_path + file_name)
    seq1 = []
    i = 0
    stot = []
    for line in f:			# read data
        if line[0] != '#':
            s = line.replace(',','')	# taking care of data seperated be ',' or by ';'
            s = s.replace(' ','')
            s = s.replace(';','')
            stot = stot + list(s)
            if stot[-1] == '\n': del(stot[-1])		# make sure there is no <CR> at the end of sequence
    f.close()

    seq1 = [x for x in stot]				# full thread
    if seq1[-1] == '\n': del(seq1[-1])			# make sure there is no <CR> at the end of sequence
    seq = [x for x in stot if x in AAlist]		# Only amino acid
    # seq_list = [aminotable[x] for x in seq]		# list of amino acid objects

    return st.join(seq,''),st.join(seq1,'')
# end of read sequence from file

def read_seq_file(file_name, file_path):
    '''
    Read an amino acid from a file
    '''
    aminotable = AminoTable()
    f = open(file_path + file_name)
    seq = []
    for line in f:			# read data
        if line[0]!= '#':
            s = line.replace(',',' ')	# taking care of data seperated be ',' or by ';'
            s = s.replace(';',' ')
            s = s.replace('',' ')
            s = s.split()
            seq = seq + s
    f.close()
    # seq_list = [aminotable[x] for x in seq]		# list of amino acid objects
    seq_s = st.join(seq,'')
    seq_s = seq_s.upper()
    return seq_s
# end of read sequence from file

def read_structure_rules(file_name, file_path):
    '''
	Read structure information, rules from a file
	'''

    f = open(file_path + file_name)
    tmp = []
    for line in f:					# read data	    
        s = line.split()
        if s[0]!='#':	
            a = ThreadStructure()
            a.PosNum = int(s[0])
            a.PointInOut = float(s[1])
            a.CornerFlag = int(s[2])
            a.CornerMarker = int(s[3])
            a.ResSize = float(s[4])
            a.ResSizeSigma = float(s[5])	
            a.TurnVolume = float(s[6])
            a.TurnVolSigma = float(s[7])
            a.Special = int(s[8])
            a.Turn = int(s[9])	
            a.DiSulDis = int(s[10])

            tmp = tmp + [a]
    f.close()	
    return tmp
# end of read sequence from file


def readProgParameters(Structure_rules):
    BaseData = BaseInputData(Structure_rules)
    # read shortest allowed loop length
    BaseData.shortl = read_var(1,1,'Enter shortest allowed loop length')
    # read longest allowed loop length  
    BaseData.longl = read_var(4,1,'Enter longest allowed loop length')
    # What range of score we would like to get information on   
    BaseData.dScore = read_var(5.0,2,'Enter results score range')   
    # Tolerance in scoring - scores that are within the tolerance will be consider to be a good path score   
    #BaseData.tolerance = read_var(0.0,2,'Enter scoring tolerance')      
    # minimum steps needed before next loop is allowed   
    BaseData.maxLoopFlag = read_var(17,1,'Enter minimum number of steps between loops')  
    # minimum steps needed before next corner cut is allowed   


    return BaseData


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


def Get_pdb_seq_info(file_name, file_path):
    '''
    read a sequance from pdb file
    '''
    namelist = AA_name_convert()
    f = open(file_path + file_name)
    seq = []
    AAnum = 0				# amino acid counter
    chain = ''				# chain name
    model_flag = True			# if more than 1 model, use only the first one
    stop_flag = False			# also for 1 model, use only the first one
    for line in f:			# read data
        s = line.split() 
        if (s[0] == 'ATOM') and model_flag:		# read only ATOM data from the first chain
            #stop_flag = True
            if AAnum == 0: chain = s[4]			# storing the name of the first chain
            if len(s[2]) < 4:
                if chain == s[4]:
                    if  AAnum < int(s[5]):	# if it is a new amino acid
                        AAnum = int(s[5])
                        seq += [namelist[s[3][-3:]]]
                else: 
                    break 
                

        #elif stop_flag and model_flag: model_flag = False	# stop recording amino acids


    f.close()
    seq_s = st.join(seq,'')
    return seq_s