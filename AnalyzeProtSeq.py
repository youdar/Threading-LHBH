'''
AnalyzeProtSeq.py is a program that run a series of tests on a protein sequence, produce plots files 
and text file that summerize the results


To run the test you need to setup the folder and job file in the Run_JobFile_ProtScan() function
and to make sure the job_file.txt is in this format:

line with "Folder::" locate the folder location. where all the proteins forders are
    Only one "Folder" line
    
line with "job::" for each protein data folder
    the proein file should have the same name as the folder

    Example:
    Folder:: C://Users//youval//Documents//Education//work biophy//Protein seq//Prion collection//
    job:: Primates_P04156_Human
    job:: Rodentia_P04925_Mouse

Note that at the moment we are only doing LBH-18
'''

import numpy as nm
import string as st
import pylab as pl
import copy as cp
import time 
from TheaderCommonFunctions import *
from ReadDate import *
from ThreadClasses import *
from AminoData import *

def Run_JobFile_ProtScan():
    job_file_path = 'C://Users//youval//Documents//Education/work biophy//Protein seq//Prion collection//'
    #job_file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//'
    #job_file_name = 'test.txt'
    job_file_name = 'job_file.txt'
    AnalyzeJob(job_file_name,job_file_path)


def AnalyzeJob(job_file_name='Job_File.txt',job_file_path=''):
    '''
    Job list file
    Doing only LBH-18

    line with "Folder::" locate the folder location. where all the proteins forders are
    Only one "Folder" line

    line with "job::" for each protein data folder
    the proein file should have the same name as the folder

    Example:
    Folder:: C://Users//youval//Documents//Education//work biophy//Protein seq//Prion collection//
    job:: Primates_P04156_Human
    job:: Rodentia_P04925_Mouse

    NO spaces in names
    line with "Done_job" represent completed jobs

    Before running the program, make sure you have a folder with the sequence file in it
    for each sequence you want to check

    '''
    # Structure rule file name and location
    s_file_name = 'StructureFile-LeftBetaH-18.txt'
    s_file_path = 'C://Users//youval//Documents//Education//work biophy//Structure and scoring rules//'       
    Structure_rules = read_structure_rules(s_file_name, s_file_path)
    ScoresP = ScoreParameters()		# Scoring parameters
    jobFileStr = job_file_path + job_file_name
    DataHeader,JobData = Get_file_data(jobFileStr)

    for i in range(len(JobData)):
        if JobData[i][0] in ['Folder','folder']:
            PathSTR = cp.copy(JobData[i][-1])
        if JobData[i][0] in ['job','Job']:
            ProtName = JobData[i][-1]
            JobData[i][0] = 'job_inprogress'
            tmpData = [':: '.join(s) + '\n' for s in JobData]
            jobFile = open(jobFileStr,'w')
            jobFile.writelines(DataHeader)
            jobFile.writelines(tmpData)
            jobFile.close()
            AnalyzeSeq(ProtName,PathSTR,Structure_rules,ScoresP)            
            DataHeader,JobData = Get_file_data(jobFileStr)
            JobData[i][0] = 'job_done'
            tmpData = [':: '.join(s) + '\n' for s in JobData]
            jobFile = open(jobFileStr,'w')
            jobFile.writelines(DataHeader)
            jobFile.writelines(tmpData)
            jobFile.close()
    print 'Job is Done'


def Get_file_data(jobFileStr):
    jobFile = open(jobFileStr)
    RawJobData = jobFile.readlines()		# All jobs information
    jobFile.close()
    JobData = []
    DataHeader = []
    for str in RawJobData:			# Maintain file header
        if str[0] == '#': DataHeader += [str]
        else:
            tmp = str[:-1].split('::')
            tmp = [x.strip() for x in tmp]
            JobData += [tmp]
    return DataHeader,JobData



def AnalyzeSeq(ProtName,PathSTR,Structure_rules,ScoresP):
    '''
    ProtName : protein name - to be printed on files and plots
    '''
    # Combine file destination path

    OutputFolderPath = PathSTR + ProtName + '//'

    # Reading the protien sequence
    seq  = read_seq_file(ProtName, OutputFolderPath)	# Get amino acid from pdb file
    # Extacting and calculating basic test info
    Hi_score = 60				# Score above which I'll pressent threads
    tol = 3					# The range of score from the run top score for which treads will be collected
    plot_filename = ProtName + '_plot_'
    seq_length = len(seq)
    fig_num = 1
    yLabelStr = 'Scores' 
    test_segment_loop = [[72,2],[72,8],[72,15],[42,2],[42,8]]
    #test_segment_loop = [[42,2]]
    #open file for sequence writing
    filestr = OutputFolderPath + ProtName + '_Threads.txt'
    LineSTR = 'Start at amino acid: '
    LoopSTR = '   Allowed loop length:'
    SegmentSTR = '    Using segment length: '

    f = open(filestr, 'w')
    f.write('Threading results collection for: ' + ProtName + '\n')
    f.write('\n')
    f.write('Score tolerance for threads collection: %.1f \n' % tol)
    f.write('\n')
    f.write('The Sequence: \n')
    seqCounter = 0
    for seqAA in seq:
        f.write(seqAA)
        seqCounter += 1
        if seqCounter%80 == 0:
            f.write('\n')
            seqCounter = 0
    f.write('\n')
    # Creating a series of plots for the whole sequence according the 
    # the test_segment_loop information
    # Note that Structure_rules can be for different structures
    for data in test_segment_loop:
        test_segment = data[0]
        tesl_loop = data[1]
        test_filename = plot_filename + str(test_segment) + '_' + str(tesl_loop)
        Output_file = OutputFolderPath + '//' + test_filename
        x,y,ScoreThreadList = Seq_file_score_LBH18(seq,Structure_rules,ScoresP,1,tesl_loop,test_segment,17,tol,Hi_score)
        # ScoreThreadList contains all the threads for which the score was more than Hi_score
        pl.figure(fig_num)
        TitleString = 'Scores for %s - steps of %i amino acids' % (ProtName,6)
        xLabelStr = 'Amino acids, Segments length: %i, Max loop: %i ' % (test_segment,tesl_loop)
        pl.title(TitleString)
        pl.xlabel(xLabelStr)
        pl.ylabel(yLabelStr)
        pl.plot(x,y,'r--o')  	# r fro red, --o both dash line and dots
        pl.ylim((0,120))
        pl.xlim((0,(len(seq)+10)))
        pl.savefig(Output_file, dpi=600, format='tif')
        fig_num += 1
        # Writing to file

        if ScoreThreadList != []:
            for tempData in ScoreThreadList:
                OutStr = LineSTR + str(tempData[0]) + SegmentSTR + str(test_segment) + LoopSTR + str(tesl_loop) 
                OutStr += '   Score: %.1f' % tempData[2]
                f.write('\n')
                f.write(OutStr + '\n')
                for th in tempData[1]:
                    f.write(th + '\n' )

    f.write('\n' + 'End')
    f.close()

def AnalyzeLBH18prot(segLength,LoopLength,seq):
    '''
    You need to be in the same folder as the protein file to run this program   
    This program is for producing protein scan plot only   
    Custom for LBH-18 only
    '''
    
    s_file_name = 'StructureFile-LeftBetaH-18.txt'
    s_file_path = 'C://Users//youval//Documents//Education//work biophy//Structure and scoring rules//'       
    Structure_rules = read_structure_rules(s_file_name, s_file_path)
    ScoresP = ScoreParameters()		# Scoring parameters
    

    # Reading the protien sequence
    OutFileName = 'plotdata'
    file_path = 'C://Users//youval//Documents//Education//Research//Paper 2 - Threading program//Plots//'
    # file_name = 'proteinSequence.txt'
    OutputFile = file_path + OutFileName +'.txt'
    # seq  = read_seq_file(file_name, file_path)
    print seq
    # Extacting and calculating basic test info
    Hi_score = 60				# Score above which I'll pressent threads
    tol = 3					# The range of score from the run top score for which treads will be collected
    plot_filename = 'test_plot_'
    seq_length = len(seq)
    yLabelStr = 'Scores' 
    test_segment_loop = [[int(segLength),int(LoopLength)]]
  

    #f.write('Threading results collection for: ' + ProtName + '\n')
    #f.write('\n')
    #f.write('Score tolerance for threads collection: %.1f \n' % tol)
    #f.write('\n')
    #f.write('The Sequence: \n')
    #seqCounter = 0
    #for seqAA in seq:
        #f.write(seqAA)
        #seqCounter += 1
        #if seqCounter%80 == 0:
            #f.write('\n')
            #seqCounter = 0
    #f.write('\n')
    # Creating a file with plot data points
    #open file for sequence writing
    data = [segLength,LoopLength]
    f = open(OutputFile, 'w')
    test_segment = data[0]
    tesl_loop = data[1]
    x,y,ScoreThreadList = Seq_file_score_LBH18(seq,Structure_rules,ScoresP,1,tesl_loop,test_segment,17,tol,Hi_score)

    
    for i in range(len(x)):
        f.write(str(x[i]))
        f.write(' ')
        f.write(str(y[i]))
        f.write('\n')
    f.close()    
    
    
# To run the job run  Run_JobFile_ProtScan()   
#Run_JobFile_ProtScan()




