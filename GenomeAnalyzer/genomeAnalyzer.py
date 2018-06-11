#!/usr/bin/env python3
# Name: Stephanie Gardner (sggardne)
# Group Members: none

from sequenceAnalysis import NucParams,FastAreader
import sys



def genomeAnalyzer():

    '''
    This genomeAnalyzer takes FastA files and prints out
    sorted relative codon information, sequence length
    information in Mb, and GC % content of the given genome.
    '''

    myReader = FastAreader('/Users/stephaniegardner/Desktop/BME160/Lab04/HomoSapiensMitochondrion.fa') #instantiation of FasAreader class
    myNuc = NucParams() #instantiation of NucParams class
    for head, seq in myReader.readFasta() :#usage of FastAreader class
        myNuc.addSequence(seq) #usage of NucParams class

    '''Sequence length: takes total nucleotide counts and converts to Mb'''
    length = myNuc.nucCount()/1000000
    print('sequence length = {:.2f}Mb'.format(length),"\n")

    '''GC Content: adds all G's and C's found in nucComp dictionary and divides by total nucleotides'''
    nucComp = myNuc.nucComposition()
    gc = nucComp['G'] + nucComp['C']
    gc = gc/myNuc.nucCount()
    print('GC Content = {:.2%}'.format(gc),"\n")

    #Individual Codon Analysis
    codonComp = myNuc.codonComposition()
    aaComp = myNuc.aaComposition()
    for codon,aa in sorted(myNuc.rnaCodonTable.items(),key=lambda t: t[1]+t[0]):
    #^^crazy lambda function Logan helped me with. It dictates what is a codon and what is an amino acid in rnaCodonTable
        total = aaComp[aa] #total number of the certain amino acid being iterated
        codonCount = codonComp[codon] #total number of codons associated with amino being iterated

        if total != 0:
                val = (codonCount/total)
                #^^takes the amount of codon for a certain amino and divdes it by total times the amino occured in the sequence
        else:
                val = (codonCount/1)
                #^^ if the amino acid was not found in the genome, its value is 0%

        print('{} : {} {:5.1f}% ({:6d})'.format(codon,aa,val*100, codonCount))


genomeAnalyzer()
