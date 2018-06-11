#!/usr/bin/env python3
# sequenceAnalysis.py
# Name: Stephanie Gardner (sggardne)
# Group Members: none

import sys


class FastAreader :
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

class NucParams:

    '''
    Creates nucleotide, codon, and amino acid dictionaries from
    a given genome/DNA string

    instantiation:
    myNuc = NucParams()
    usage:
    seq is the desired DNA string,
    myNuc.addSequence(seq)
    '''

    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }

    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    allowedBases = 'ATCGUN'

    aaComp = {}
    nucComp = {}
    codonComp = {}
    codon = []

    def __init__ (self):
        '''init method builds dictionaries from rnaCodonTable dictionary.'''

        self.nucComp = {nuc:0 for nuc in self.allowedBases.upper()} #nucleotide composition dictionary

        for codon,aminos in self.rnaCodonTable.items():
            self.codonComp[codon] = 0 #codon composition dictionary
            self.aaComp[aminos] = 0 #amino acid dictionary


    def addSequence(self, addedSeq):
        '''add sequence adds values for each base, codon, and amino acid in the data.'''

        for aminos in addedSeq.upper(): #uppercases the sequence string
            if aminos in self.allowedBases: #makes sure all bases are valid
                self.nucComp[aminos] += 1 #adds 1 to value for each base


        for pos in range(0,len(addedSeq.upper()),3):
            self.codon = addedSeq.upper()[pos:pos + 3] #parses the sequence string by 3 chars to create codons
            rnaCodon = self.codon.replace('T','U') #makes sure all codons are rna codons
            if rnaCodon in self.codonComp:
                self.codonComp[rnaCodon] += 1 #adds 1 to the value of each codon iterated
                self.aaComp[self.rnaCodonTable[rnaCodon]] += 1 #adds 1 to the value of each amino acid iterated


    '''the rest of these functions simply return
        the dicitonaries that can be used in other programs'''

    def aaComposition(self):

        return self.aaComp

    def nucComposition(self):

        return self.nucComp

    def codonComposition(self):

        return self.codonComp

    def nucCount(self):
        return (sum(self.nucComp.values()))
