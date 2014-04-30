import numpy as nm
import string as st
from TheaderCommonFunctions import *
from ReadDate import *
from ThreadClasses import *
from AminoData import *
import pylab as pl
import copy as cp

aminotable = AminoTable()	# Dictionary of amino acids
AAlist = AminoList()
# Read structure rules
Structure_rules,s_file_path,s_file_name = structure_file_info()

# read some program ristrictions
BaseData = BaseInputData(Structure_rules)
build_Reverse_in_sidechain_list(BaseData,Structure_rules)

# scoring parameters ratios
ScoresP = ScoreParameters()

#The threads we want to compare
#seqA1 = 'RNFYAN|FNLTIV|D(DY)TVTIG|DNVLIA|PNVTLS|V(TGHPVHHELRKNGEMYSF)PITIG|NNVWIG|SHVVIN'
#seqB1 = 'RNFYAN|FNLTIV|D(DY)TVTIG|DNVLIA|PNVTLS|V(TGHPVHH)ELRKN|GEMYSF|-PITIG|NNVWIG|SHVVIN'

#seqA1 = 'SHVVVN|GHTKIG|RDNEIY|QFASIG|(EVNQDLKYAGEP)TRVEIG|DRNRIR|ESVTIH|R(GTVQGGG)LTKVG|SDNLLM|INAHIA|HDCTVG|NRCILA|NNATLA|GHVSVD|DFAIIG|GMTAVH'
#seqB1 = 'SHVVVN|GHTKIG|RDNEIY|QFASIG|-EVNQD|(LKYAGEP)TRVEIG|DRNRIR|ESVTIH|R(GTVQGGG)LTKVG|SDNLLM|INAHIA|HDCTVG|NRCILA|NNATLA|GHVSVD|DFAIIG|GMTAVH'

#seqA1 = 'PYAHIR|PNSSLG|AQVHIG|NFVEVK|-GSSIG|ENTKAG|HLTYIG|-NCEVG|SNVNFG|AGTITV|(NYDGKNK)YKTVIG|NNVFVG|SNSTII|APVELG'
#seqB1 = 'PYAHIR|PNSSLG|AQVHIG|NFVEVK|-GSSIG|ENTKAG|HLTYIG|-NCEVG|SNVNFG|(A)GTITVN|YDGKNK|YKTVIG|NNVFVG|SNSTII|APVELG'

seqA1 = 'NGVSFV|(NP)EATYID|IDVEIA|SEVQIE|ANVTLK|GQTKIG|AETVLT|NGTYVV|-DSTIG|AGAVIT|-NSMIE|-ESSVA|DGVIVG|PYAHIR|PNSSLG|AQVHIG|NFVEVK|-GSSIG|ENTKAG|HLTYIG|-NCEVG|SNVNFG'
seqB1 = 'NGVSFV|(NP)EATYID|IDVEIA|SEVQIE|ANVTLK|GQTKIG|AETVLT|NGTYVV|D(STI)GAGAV|ITNSMI|EESSVA|DGVIVG|PYAHIR|PNSSLG|AQVHIG|NFVEVK|-GSSIG|ENTKAG|HLTYIG|-NCEVG|SNVNFG'

# Collect clean amino acid sequence and combine to 1 string
seqA = [x for x in seqA1 if x in AAlist]	
seqA = st.join(seqA,'')
seqB = [x for x in seqB1 if x in AAlist]
seqB = st.join(seqB,'')

# Get the correct thread score
FixScoreA = GetThreadScores(seqA,seqA1,Structure_rules,ScoresP,BaseData)
print '======== Start scoring second sequence  ============'
FixScoreB = GetThreadScores(seqB,seqB1,Structure_rules,ScoresP,BaseData)


# Printing info
print 'Total score results'
print '%.1f :: %s' % (FixScoreA,seqA1)
print '%.1f :: %s' % (FixScoreB,seqB1)
