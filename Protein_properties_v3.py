class AminoAcid:
    def __init__(self,name='AA'):
        self.name = name
        self.name3L = ''
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.ResWeight = 0 		# residue weight (weight - 56) backbone weight is 56
        self.ResVol = 0			# Residue volume from http://prowl.rockefeller.edu/aainfo/volume.htm
        self.SideChainVol = 0		# Side Chain volume is evaluated as ResVol - 0.9 Gly.ResVol
        self.Hydropathy = 0		# Hydropathy index
        self.n1 = 0			
        self.n2 = 0 
        # n values
        # -1 when the amino acid (AA) residue has an N donor, short residue.	
        # -2 when the AA residue has an O acceptor, short residue.
        # -3 when the AA residue has an N donor, long residue that able to bond across two turns.
        # -5 when the AA residue has an O acceptor, long residue that able to bond across two turns.
        # -7 when it is a Cystine(C) 
        # 0 when bond is not possible.
        # 1 when this N or O participating in a bond.   
        # A residu can only participate in one side-chain bond. So when a bond is created
        # for example with n1, n1 get the bond value and n2 will be assigned 0


    def __mul__ (self,other):
        # Evaluating side chain interaction
        Prod = self.donor * other.acceptor
        return Prod


# ############   Non Polar, HydroPhobic   ###########
class Ala(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'A')
        # Alanine	 
        # ###
        # CH3-CH(NH2)-COOH
        #
        #
        # Molecular weight 89.09 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.616 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 1.8 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point (when protonation accure) pH 6.01
        # pKa( alpha-COOH) 2.35 
        # pKa( alpha-NH2) 9.87
        # CAS # 56-41-7
        # PubChem ID 5950
        #
        self.name3L = 'ALA'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0	
        self.Hydropathy = 1.8
        self.ResWeight = 33
        self.ResVol = 88.6
        self.SideChainVol = 88.6-54.1


class Val(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'V')
        # Valine
        # #########
        # (CH3)2-CH-CH(NH2)-COOH
        #
        #
        # Essential AA (cannot be synthesized by humans)
        # Molecular weight 117.15 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.825 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 4.2 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.00
        # pKa( alpha-COOH) 2.39
        # pKa( alpha-NH2) 9.74
        # CAS # 72-18-4
        # PubChem ID 1182
        #
        self.name3L = 'VAL'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.Hydropathy = 4.2
        self.ResWeight = 61
        self.ResVol = 140.0
        self.SideChainVol = 140-54.1


class Leu(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'L')
        # Leucine	 
        # #############
        # (CH3)2-CH-CH2-CH(NH2)-COOH
        #
        #
        # Essential AA
        # Molecular weight 131.18 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.943 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 3.8 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.01
        # pKa( alpha-COOH) 2.33
        # pKa( alpha-NH2) 9.74
        # CAS # 61-90-5
        # PubChem ID 6106
        #
        self.name3L = 'LEU'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.Hydropathy = 3.8
        self.ResWeight = 75
        self.ResVol = 166.7
        self.SideChainVol = 166.7-54.1


class Ile(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'I')
        # Isoleucine
        # ###############
        # CH3-CH2-CH(CH3)-CH(NH2)-COOH
        #
        #
        # Essential AA
        # Molecular weight 131.18 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.943 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 4.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.05
        # pKa( alpha-COOH) 2.33
        # pKa( alpha-NH2) 9.74
        # CAS # 61-90-5
        # PubChem ID 6106
        #
        self.Hydropathy = 4.5
        self.ResWeight = 75
        self.name3L = 'ILE'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 166.7
        self.SideChainVol = 166.7-54.1

class Phe(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'F')
        # Phenylalanine
        # ######
        # Ph-CH2-CH(NH2)-COOH
        # The residue Ph-CH2 : C6H5-CH2  benzyl
        #
        # Essential AA
        # Molecular weight 165.19 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 1 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 2.8 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.49
        # pKa( alpha-COOH) 2.20
        # pKa( alpha-NH2) 9.31
        # CAS # 63-91-2
        # PubChem ID 994
        #
        self.Hydropathy = 2.8
        self.ResWeight = 109
        self.name3L = 'PHE'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 189.9
        self.SideChainVol = 189.9-54.1


class Trp(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'W')
        # Tryptophan 
        # ##############
        # Ph-NH-CH=C-CH2-CH(NH2)-COOH
        # |________|
        #
        # contains an indole functional group.
        # aromatic heterocyclic organic compound
        # It has a bicyclic structure, consisting of a six-membered benzene ring fused to 
        # a five-membered nitrogen-containing pyrrole ring
        #
        # Essential AA
        # Molecular weight 204.23 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.878 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -0.9 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.89
        # pKa( alpha-COOH) 2.46
        # pKa( alpha-NH2) 9.41
        # CAS # 73-22-3
        # PubChem ID 6305
        #
        self.Hydropathy = -0.9
        self.ResWeight = 148
        self.name3L = 'TRP'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 227.8
        self.SideChainVol = 227.8-54.1

class Met(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'M')
        # Methionine 
        # ############
        # CH3-S-(CH2)2-CH(NH2)-COOH
        # sulfur-containing residue
        # methyl donor  R-CH3
        # methionine is incorporated into the N-terminal position of all proteins 
        # in eukaryotes and archaea during translation, although it is usually removed 
        # by post-translational modification
        #
        # Essential AA
        # Molecular weight 149.21 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.738 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 1.9 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.74
        # pKa( alpha-COOH) 2.13
        # pKa( alpha-NH2) 9.28
        # CAS # 63-68-3
        # PubChem ID 876
        #
        self.Hydropathy = 1.9
        self.ResWeight = 93
        self.name3L = 'MET'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 162.9
        self.SideChainVol = 162.9-54.1


class Pro(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'P')
        # Proline
        # **********
        # NH-(CH2)3-CH-COOH
        # |_________|
        # Side chain bond to C alpha
        # exceptional conformational rigidity 
        # usually solvent-exposed.
        # lacks a hydrogen on the amide group, it cannot act as a hydrogen bond donor, 
        # only as a hydrogen bond acceptor.
        #
        # Molecular weight 115.13 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.711 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -1.6 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.30
        # pKa( alpha-COOH) 1.95
        # pKa( alpha-NH2) 10.64
        # CAS # 147-85-3
        # PubChem ID 614
        #
        self.Hydropathy = -1.6
        self.ResWeight = 59
        self.name3L = 'PRO'
        self.Hydrophobic = 1  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 0
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0,3:-3,4:-3,5:-3,6:-2}	# special value scores
        self.n1 = 0			
        self.n2 = 0 
        self.ResVol = 112.7
        self.SideChainVol = 112.7-54.1

# ############ Non Polar Uncharged	###########

class Gly(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'G')	
        # 
        # NH2-CH2-COOH
        # 
        # Molecular weight 75.07 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.501 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -0.4 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.06
        # pKa( alpha-COOH) 2.35
        # pKa( alpha-NH2) 9.78
        # CAS # 56-40-6
        # PubChem ID 750
        self.Hydrophobic = 0		# 1: Hydrophobic, 0: Hydrophilic
        self.Hydropathy = -0.4
        self.ResWeight = 19
        self.name3L = 'GLY'
        self.SpecialRes = {0:0,3:-3,5:-3}	# special value scores
        self.ResVol = 60.1
        self.SideChainVol = 60.1-54.1

# ############ Polar Uncharged	###########

class Ser(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'S')
        # Serine
        # ######
        # HO-CH2-CH(NH2)-COOH
        # 
        # Molecular weight 105.09 Da
        # Ploar
        # Acidity - Natural
        # Hydrophobicity 0.359 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -0.8 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.68
        # pKa( alpha-COOH) 2.19
        # pKa( alpha-NH2) 9.21
        # CAS # 56-45-1
        # PubChem ID 617
        #
        self.Hydropathy = -0.8
        self.ResWeight = 49
        self.name3L = 'SER'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 89.0
        self.SideChainVol = 89-54.1


class Thr(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'T')
        # Threonine
        # ##########
        # CH3-CH(OH)-CH(NH2)-COOH
        # bearing an alcohol group
        #
        # Essential AA
        # Molecular weight 119.12 Da
        # Ploar
        # Acidity - Natural
        # Hydrophobicity 0.450 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -0.7 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.60
        # pKa( alpha-COOH) 2.09
        # pKa( alpha-NH2) 9.10
        # CAS # 72-19-5
        # PubChem ID 6288
        #
        self.Hydropathy = -0.7
        self.ResWeight = 63
        self.name3L = 'THR'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 116.1
        self.SideChainVol = 116.1-54.1

class Cys(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'C')
        # Cysteine
        # ######
        # HS-CH2-CH(NH2)-COOH
        # thiol (R-S-H) side chain 
        # Has Sulfur in side chain
        #
        # Molecular weight 121.16 Da
        # Ploar
        # Acidity - Natural
        # Hydrophobicity 0.680 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index 2.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.05
        # pKa( alpha-COOH) 1.92
        # pKa( alpha-NH2) 10.70
        # CAS # 59-90-4
        # PubChem ID 5862
        #
        self.Hydropathy = 2.5
        self.ResWeight = 65
        self.name3L = 'CYS'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.n1 = -7			
        self.n2 = 0
        self.SpecialRes = {0:0}	# special value scores
        self.ResVol = 108.5
        self.SideChainVol = 108.5-54.1



class Tyr(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'Y')
        # Tyrosine  
        # ###########
        # HO-p-Ph-CH2-CH(NH2)-COOH
        # 
        # Molecular weight 181.19 Da
        # Non ploar
        # Acidity - Natural
        # Hydrophobicity 0.880 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -1.3 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.64
        # pKa( alpha-COOH) 2.20
        # pKa( alpha-NH2) 9.21
        # CAS # 60-18-4
        # PubChem ID 1153
        #
        self.Hydropathy = -1.3
        self.ResWeight = 125
        self.name3L = 'TYR'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = 0			
        self.n2 = 0
        self.ResVol = 193.6
        self.SideChainVol = 193.6-54.1

class Asn(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'N')
        # Asparagine
        # ##########
        # H2N-CO-CH2-CH(NH2)-COOH
        # N Donor - NH2
        # 
        # has carboxamide as the side chain's functional group(R-CO-NH2)
        # side chain can form hydrogen bond interactions with the peptide backbone
        # often found near the beginning and the end of alpha-helices, 
        # and in turn motifs in beta sheets.
        # 
        # Molecular weight 132.12 Da
        # Ploar
        # Acidity - Natural
        # Hydrophobicity 0.236 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.41
        # pKa( alpha-COOH) 2.14
        # pKa( alpha-NH2) 8.72
        # CAS # 70-47-3
        # PubChem ID 236
        #
        self.Hydropathy = -3.5
        self.ResWeight = 76
        self.name3L = 'ASN'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = -1			
        self.n2 = -2
        self.ResVol = 114.1
        self.SideChainVol = 114.1-54.1


class Gln(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'Q')
        # Glutamine  
        # #############
        # H2N-CO-(CH2)2-CH(NH2)-COOH
        # N Donor - NH2 
        #
        # Molecular weight 146.14 Da
        # Ploar
        # Acidity - Natural
        # Hydrophobicity 0.251 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 5.65
        # pKa( alpha-COOH) 2.17
        # pKa( alpha-NH2) 9.13
        # CAS # 56-85-9
        # PubChem ID 5950
        #
        self.Hydropathy = -3.5 
        self.ResWeight = 90
        self.name3L = 'GLN'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.n1 = -1			
        self.n2 = -2
        self.SpecialRes = {0:0}	# special value scores
        self.ResVol = 143.8
        self.SideChainVol = 143.8-54.1

#  ##########   Polar Acidic   ###########

class Asp(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'D')
        # Aspartic acid
        # ########
        # HOOC-CH2-CH(NH2)-COOH
        # 
        # Molecular weight 133.10 Da
        # Ploar
        # Acidity - Acidic
        # Hydrophobicity 0.028 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 2.85
        # pKa( alpha-COOH) 1.99
        # pKa( alpha-NH2) 9.90
        # CAS # 56-84-8
        # PubChem ID 5960
        #
        self.Hydropathy = -3.5 
        self.ResWeight = 77
        self.name3L = 'ASP'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 1
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop self.loop = 0
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = -2			
        self.n2 = 0
        self.ResVol = 111.1
        self.SideChainVol = 111.1-54.1

class Glu(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'E')
        # ###########
        # HOOC-(CH2)2-CH(NH2)-COOH
        # 
        # Molecular weight 147.13 Da
        # Ploar
        # Acidity - Acidic
        # Hydrophobicity 0.043 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 3.15
        # pKa( alpha-COOH) 2.10
        # pKa( alpha-NH2) 9.47
        # CAS # 56-86-0
        # PubChem ID 611
        #
        self.Hydropathy = -3.5 
        self.ResWeight = 91
        self.name3L = 'GLU'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 1
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop self.loop = 0
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = -2			
        self.n2 = 0
        self.ResVol = 138.4
        self.SideChainVol = 138.4-54.1


# ##############  Polar Basic   #############

class Lys(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'K')
        # Lysine
        # ##########
        # H2N-(CH2)4-CH(NH2)-COOH
        # often participates in hydrogen bonding 
        # N Donor - NH2
        #
        # Essential AA
        # Molecular weight 146.19 Da
        # Ploar
        # Acidity - Basic
        # Hydrophobicity 0.283 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.9 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 6.90
        # pKa( alpha-COOH) 2.16
        # pKa( alpha-NH2) 9.06
        # CAS # 56-87-1
        # PubChem ID 866
        #
        self.Hydropathy = -3.9
        self.ResWeight = 90
        self.Hydropathy = -3.9
        self.ResWeight = 90
        self.name3L = 'LYS'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 1
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop self.loop = 0
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = -1			
        self.n2 = 0
        self.ResVol = 168.6
        self.SideChainVol = 168.6-54.1


class Arg(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'R')
        # Arginine
        # ###################
        # HN=C(NH2)-NH-(CH2)3-CH(NH2)-COOH
        # N Donor - NH2
        # 
        # Molecular weight 174.20 Da
        # Ploar
        # Acidity - Basic (strong)
        # Hydrophobicity 0.000 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -4.5 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point (when protonation accure) pH 10.76
        # pKa( alpha-COOH) 1.82
        # pKa( alpha-NH2) 8.99
        # CAS # 74-79-3
        # PubChem ID 5950
        #
        self.Hydropathy = -4.5 
        self.ResWeight = 118
        self.name3L = 'ARG'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 1
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop self.loop = 1
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = -1			
        self.n2 = -3
        self.ResVol = 173.4
        self.SideChainVol = 173.4-54.1

class His(AminoAcid):
    def __init__(self):   
        AminoAcid.__init__(self,'H')
        # Histidine 
        # ################
        # NH-CH=N-CH=C-CH2-CH(NH2)-COOH
        # |__________|
        # N Donor - NH
        # The imidazole side chain has two nitrogens with different properties
        # 
        # Molecular weight 155.15 Da
        # Ploar
        # Acidity - Basic (week)
        # Hydrophobicity 0.165 (Analytical Bio chemistry 193:11,72-82 Elsevier 1991)
        # Hydrophathy index -3.2 (J.Mol.Bio(1982) 157, 105-132)
        # Isoelectric point 7.60
        # pKa( alpha-COOH) 1.80
        # pKa( alpha-NH2) 9.33
        # CAS # 71-00-1
        # PubChem ID 773
        #
        self.Hydropathy = -3.2
        self.ResWeight = 99
        self.name3L = 'HIS'
        self.Hydrophobic = 0  		# 1: Hydrophobic, 0: Hydrophilic
        self.charge = 0.5
        self.polar = 1
        self.corner = 0			# Would prefer to be at a corner : give positive value
        self.loop = 0			# cost/benefit when on a loop
        self.size = 0			# Residue size   (0:1)   0:ignor size, 1:Large residue
        self.SpecialRes = {0:0}		# Special characteristic of residue
        self.n1 = -1			
        self.n2 = 0
        self.ResVol = 153.2
        self.SideChainVol = 153.2-54.1