from Protein_properties_v3 import *
from ThreadClasses import *

def AminoTable():
    '''
    Dictionary of amino acids
    '''
    aminotable = {'A':Ala(),'V':Val(),'L':Leu(),'I':Ile(),'F':Phe(),'W':Trp(),'G':Gly(),'S':Ser(),
                  'T':Thr(),'C':Cys(),'Y':Tyr(),'N':Asn(),'M':Met(),'P':Pro(),
                  'D':Asp(),'E':Glu(),'Q':Gln(),'K':Lys(),'R':Arg(),'H':His()}
    return aminotable

def AminoList():
    AAlist = ['A','V','F','I','L','P','M','D','E','K','R',
              'S','T','Y','C','N','Q','H','W','G']
    return AAlist

def AA_name_convert():
    Namelist = {'ALA':'A','VAL':'V','PHE':'F','ILE':'I','LEU':'L','PRO':'P','MET':'M','ASP':'D',
                'GLU':'E','LYS':'K','ARG':'R','SER':'S','THR':'T','TYR':'Y','CYS':'C',
                'ASN':'N','GLN':'Q','HIS':'H','TRP':'W','GLY':'G'}

    return Namelist

def AminoAcidtoNum():
    '''
    maps amino acid to numbers. there is no significans to the order
    '''
    AAtoNum = {'A':1,'V':2,'L':3,'I':4,'F':5,'W':6,'G':7,'S':8,
               'T':9,'C':10,'Y':11,'N':12,'M':13,'P':14,
               'D':15,'E':16,'Q':17,'K':18,'R':19,'H':20}
    return AAtoNum

def NumtoAminoAcid():
    '''
    maps amino acid to numbers. there is no significans to the order
    '''
    NumtoAA = {1:'A',2:'V',3:'L',4:'I',5:'F',6:'W',7:'G',8:'S',
               9:'T',10:'C',11:'Y',12:'N',13:'M',14:'P',
               15:'D',16:'E',17:'Q',18:'K',19:'R',20:'H'}
    return NumtoAA