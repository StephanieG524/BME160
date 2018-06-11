#!/usr/bin/env python3
# findUnique.py
# Name: Stephanie Gardner (sggardne)
# Group Members: Lisa (Ldefeo)

from sequenceAnalysis import FastAreader, NucParams


class tRNA:

'''
this class creates a power set of only essential and
unique tRNAvsequences from a file input

instantiation:
for eachtRNA in tRNA.tRNAobjects:
    eachtRNA.buildUniques()
    eachtRNA.findEssentials()

usage:
in command line -
python findUnique.py < fileName
'''
    tRNAobjects = list() #holds each tRNA object that has its own head, seq and power set


    def __init__(self,head,seq):

        ''''when u inizialize me you give me a sequence and i make the
        power set for myself and append it to the tRNA objets list to
        call onto later'''

        self.powerSet = set()
        for i in range(len(seq)): #makes the powerSet for each tRNA
            for j in range(i+1,len(seq)+1):
                self.powerSet.add(seq[i:j])

        tRNA.tRNAobjects.append(self)

        self.head = head
        self.seq = seq
        self.uniqueSubs = set()
        self.essentialSubs = set()



    def buildUniques(self):

        '''
        finds all the sub sequences that are not found in any other tRNA set
        '''
        allButMe = set()

        for eachtRNA in tRNA.tRNAobjects:
            if eachtRNA is not self:
                allButMe = allButMe.union(eachtRNA.powerSet)

        self.uniqueSubs = self.powerSet - allButMe


    def findEssentials(self):
        '''
        finds and saves the smallest posisble unique sequence
        '''

        nonEssential = set()

        for element in self.uniqueSubs:
            i = self.seq.find(element)
            stop = i + len(element)
            bigL = i - 1
            bigR = stop + 1
            '''Lisa helped me with identifying what should be added to the non essential list'''
            if not bigL < 0:
                nonEssential.add(self.seq[bigL:stop])
            if not bigR < 0:
                nonEssential.add(self.seq[i:bigR])

        self.essentialSubs = self.uniqueSubs - (nonEssential)



def main():

    tRnaFinder = FastAreader('')

    for head,seq in tRnaFinder.readFasta():
        allPowerSets = tRNA(head,seq) #makes powerSet for each trna sequence

    sortTrnas = sorted(tRNA.tRNAobjects, key = lambda t:t.head) #sorts the tRNAs so they can be iterated through alphabetically

    for eachtRNA in tRNA.tRNAobjects:
        eachtRNA.buildUniques()
        eachtRNA.findEssentials()

    for everyTrna in sortTrnas:
        print(everyTrna.head)
        print(everyTrna.seq)
        sortedtRna = sorted(everyTrna.essentialSubs, key=lambda x:everyTrna.seq.find(x)) #Lisa also here
        for e in sortedtRna:
            position = everyTrna.seq.find(e)
            print('{}{}'.format('.' * position, e))



if __name__=="__main__":
    main()
