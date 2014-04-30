import numpy as nm
import string as st
from TheaderCommonFunctions import *
from ReadDate import *
from ThreadClasses import *
from AminoData import *
import pylab as pl
import copy as cp

file_name = 'URE2-1_9-10000.pdb'
file_path = 'C:\Users\youval\Documents\Education\work biophy\Ure2 Paper\URE2-1'

file_path = convertPathNames(file_path)
print file_path + file_name

print Get_pdb_seq_info(file_name, file_path)

