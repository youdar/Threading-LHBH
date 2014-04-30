import numpy as nm
import string as st
from TheaderCommonFunctions import *
from ReadDate import *
from ThreadClasses import *
from AminoData import *
import pylab as pl
import copy as cp

#file_name = '1TUP.pdb'
#file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//pdb files//'
#file_path = 'C://Users//youval//Documents//Education//work biophy//Protein seq//Diseas related proteins//p53//'


#seq = Get_pdb_seq_info(file_name, file_path)
#print 'seq : ',seq
#print 'seq length = ',len(seq)

#x1,x2 = read_seq_file('1LXA', 'C://Documents and Settings//youval//My Documents//Education//work biophy//Protein seq//Original seq//')
#print 'x1',x1

filePath = 'C:\Users\youval\Documents\Education\work biophy\Ure2 Paper\URE2-1'
newfilePath = ''

for i in filePath:
    if i == '\\':
        newfilePath += '//'
    else:
        newfilePath += i
        
        
print new
    
    
    