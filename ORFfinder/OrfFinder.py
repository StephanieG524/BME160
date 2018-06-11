# #!/usr/bin/env python3
# OrfFinder.py
# Stephanie Gardner (sggardne)
# Group Members: none

from sequenceAnalysis import NucParams,FastAreader

'''
OrfFinder.py takes fa files and collectes open reading frame data.
This includes: Frame(+ denotes top strand/ - denotes bottom), start
postion, end position, and length of ORF's.
Hanging start and stop ORF's are also included in this program.
Defaults are implemented but can be changed by the user through the
command line.
'''

class CommandLine() :
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - this program finds ORFs in a DNA fasta file' ,
                        epilog = 'Program epilog - options for codons and ORF length are included', add_help = True, prefix_chars = '-',
                        usage = '%(prog)s input output -option1[default]' )
        self.parser.add_argument('inFile', action = 'store',
                                  help='input file name')
        self.parser.add_argument('outFile', action = 'store',
                                  help='output file name')
        self.parser.add_argument('-lG', '--longestGene', action = 'store',
                                 nargs='?', const=True, default=False,
                                 help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int,
                                 choices= (100,200,300,500,1000),
                                 default=100, action = 'store',
                                 help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append',
                                 choices = ['ATG','TTG','GTG'],
                                 default = ['ATG'],nargs='?',
                                 help='start Codon') #allows multiple options
        self.parser.add_argument('-t', '--stop', action = 'append',
                                 default = ['TAG','TGA','TAA'],
                                 nargs='?', help='stop Codon')
        self.parser.add_argument('-v', '--version', action='version',
                                 version='%(prog)s 0.1')
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

class OrfFinder():

    '''
    OrfFinder class finds and inputs Orf's into a dictionary
    and writes them to a file (type and name of file input by user in
    main class but has defaults)

    instantiation:
    orfReader = FastAreader(myCommandLine.args.inFile)
    with open(myCommandLine.args.outFile, 'w') as opFile:
        for head, seq in orfReader.readFasta():

    usage:
    finding = OrfFinder(bases,head,myCommandLine.args.minGene,myCommandLine.args.longestGene,myCommandLine.args.start,reverseStrand)
    '''


    stopCodons = ['TAG','TAA','TGA'] #stop codons for top DNA strand (+)


    def __init__(self, seq, header,minGene=0,start=['ATG'],longestGene=False,revSeq=False):
        self.seq = seq
        self.header = header
        self.minGene = minGene
        self.longestGene = longestGene
        self.revSeq = revSeq

        self.start = start
        self.start = ['ATG']

        self.minGene = minGene


        self.frameList = [] # has all the frames to print

    def findOrf (self,frameFile,reverse=False):
        ''' finds and labels where the ORFs are in the DNA on both the
            top and bottom strand and passes the info to a method
            that saves it into a dictionary'''

        startList = [] #make a list to save starts in
        if reverse:
            dataSeq = self.revSeq #use revSeq if reverse == True
        else:
            dataSeq = self.seq

        for frame in (0,1,2): #there are three frames to iterate through in DNA
            startAppended = False #no start found yet
            for pos in range(frame, len(dataSeq), 3): #searches through the dataSeq by 3 positions starting at the frame
                codon = ''.join(dataSeq[ pos : pos+3]) #makes each codon into a string
                if codon in self.start:
                    startList.append(pos) #if its a start codon save the information
                    startAppended = True
                if codon in self.stopCodons and startAppended: #if we have passed a start:
                    startAppended = False
                    if reverse: #pass index information to saveOrfDict
                        self.saveOrfDict(len(dataSeq) - (pos+3) + 1, len(dataSeq) - (startList[0]+1) + 1,pos-(startList[0])+3,-(frame+1),frameFile)
                    else: #pass index information to saveOrfDict
                        self.saveOrfDict(startList[0]+1,pos+3,pos-(startList[0])+3,frame+1,frameFile)
                    startList = [] #reset startList


    def findReverse (self,frameFile):
        ''' passes the reverse compliment strand into findOrf method'''

        self.findOrf(frameFile,reverse=True)


    def saveOrfDict(self,start, stop, length, frame,frameFile):
        ''' saves data from findOrf method into a dictionary with appropriate labeling'''

        self.frameList.append(dict(frame=frame,start=start,stop=stop,length=length))

    def writeFrameFile(self,frameFile):
        ''' this method sorts the frameList and writes it to a file
            from largest to smallest ORF length with the header'''

        frameFile.write('{} \n'.format(self.header))

        for frame in sorted(self.frameList, key=lambda k: k['length'], reverse=True):
            formatString = '{:+d} {:>5d}..{:>5d} {:>5d}\n'.format(frame['frame'],frame['start'],frame['stop'],frame['length'])
            frameFile.write(formatString)


def main(inCL=None):
    '''
    Find some genes.
    '''
    if inCL is None:
        myCommandLine = CommandLine(['tass2.fa',
                                    'tass2ORFdata-ATG-100.txt',
                                    '--longestGene'])

                                     #create default values
    else :
        myCommandLine = CommandLine(inCL)


    #open file created
    orfReader = FastAreader(myCommandLine.args.inFile)
    with open(myCommandLine.args.outFile, 'w') as opFile:

        #loop through the sequences within file and calculate ORFs then STDOUT them
        for head, seq in orfReader.readFasta():

            #reverses strand to be read
            bases = list(seq)
            reverseStrand = reverseDNA(bases)

            #Puts seq through program and outputs in frame file
            finding = OrfFinder(bases,head,myCommandLine.args.minGene,myCommandLine.args.longestGene,myCommandLine.args.start,reverseStrand)
            finding.findOrf(opFile)
            finding.findReverse(opFile)
            finding.writeFrameFile(opFile)



def reverseDNA(bases):
    '''
    returns the reverse compliment/bottom strand to be read by OrfFinder class.
    '''
    revBases = []
    for base in bases:
        if base == 'A':
            revBases.append('T')
        if base == 'T':
            revBases.append('A')
        if base == 'G':
            revBases.append('C')
        if base == 'C':
            revBases.append('G')

    return revBases[::-1] #after all the bases have been switched this reverses the whole list



if __name__ == "__main__":
    main()
