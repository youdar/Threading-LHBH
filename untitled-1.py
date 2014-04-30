'''
threadQualityPenalty(t.thread)
Giving penalty for three or more consecutive abnormalities
Giving penalty if the ratio between the total number of corner cuts and loops divided 
by the total corners is larger than 0.4
'''

PenaltyCount = 0
nConsecutive = 0
turnFlag = 0
nCornerCut = 0
nLoop = 0 
nCorners = 0
nAbnormality = 0
QratioFlag = 0

thread = 'KLLDMA|-AQIAE|GMAFIE|ERNYIH|RDLRAA|-NILVS|DTLSCK|IADFG|-LARLI|-NEYTAR|-E'

for i in thread:
    if i=='|':
        nCorners += 1
        if turnFlag == 0:
            # reset structure abnormality counter
            print 'nConsecutive',nConsecutive
            if nConsecutive > 2:
                # collect abnormality count
                nAbnormality += (nConsecutive - 2)
            nConsecutive = 0
        turnFlag = 0
    if i=='-':
        nCornerCut += 1
        nConsecutive += 1
        turnFlag = 1
    if i=='(':
        nLoop += 1
        nConsecutive +=1
        turnFlag = 1
        
if nConsecutive > 2:
    nAbnormality += (nConsecutive - 2)      
# Quality measurments
if nCorners > 0:
    Qratio = float(nCornerCut + nLoop)/float(nCorners)
if Qratio > 0.4:
    # to many abnormalities
    QratioFlag = 1

print nCorners,(nCornerCut + nLoop)
print QratioFlag,nAbnormality