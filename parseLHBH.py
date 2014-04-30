import numpy as nm

def ParseLhbh(seq1):
    '''
    Divide a thread to overlaping thread segments
    
    Use for producing LHBH thread segments for scoring parameter evaluation
    
    At each step jump 3 corners forward
    12 corner marks in each segment
    '''
    CornerMarkNum = 12
    iStep = CornerMarkNum + 1
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    print 'List of thread segments'
    print seq1
    cornerCounterList = [-1]
    for i in range(len(seq1)):       
        if seq1[i] == '|':
            cornerCounterList.append(i)           
        i += 1
    cornerCounter = len(cornerCounterList) - 1
    i = 0
    while cornerCounter > iStep:
        print seq1[0 + cornerCounterList[i] + 1: cornerCounterList[i+iStep]]
        i += 3
        cornerCounter -= 3
        if cornerCounter < iStep:
            print seq1[0 + cornerCounterList[i] + 1::]
    print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
# end of ParseLhbh
    
    
    
    
ParseLhbh('sdfs|dfsdf|sdfsd|fsdf|sdfs|dfsdf|sdfs|dfsd|sdfgsdfg|sdfsdf|sdfsdf|sdfsdf|2342|cvbcxb|23r32|dfgd|dgdg|hgf')
    
    