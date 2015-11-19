#Chris Dale
#Bioinformatics HW2
#Must be opened in Python3

import string

#generic function for producing string and list of book
def openGeneFile(nameFile):
    #returns a list of all of the characters in the input file
    theSequence=open(nameFile,"r",encoding="utf-8")
    theSeqString=theSequence.read()  # read the entire file
    theSequence.close()
    theSeqList=list(theSeqString)
    for i in range(len(theSeqList)-1,0,-1):
        if theSeqList[i].isalpha()==False:
            del theSeqList[i]
    return(theSeqList)

def giveReverseCompliment(inSequence):
    #returns the reverse compliment of the entire sequence
    reverseComp=[]
    for base in range(len(inSequence)):
        if inSequence[base]=="A":
            reverseComp+="T"
        elif inSequence[base]=="T":
            reverseComp+="A"
        elif inSequence[base]=="C":
            reverseComp+="G"
        elif inSequence[base]=="G":
            reverseComp+="C"
    return(reverseComp[::-1]) #anyList[::-1] will give a reversed version of anyList


def giveReadFrames(inSequence):
    #returns 3 forward reading frames
    leng=len(inSequence)
    seq1=inSequence[0:leng-(leng%3)]
    #seq2=inSequence[1:leng-1-(leng%3)] #what it was initally
    #seq3=inSequence[2:leng-2-(leng%3)] #what it was initally
    #seq2=inSequence[1:((leng-1)%3+1)] #broken
    #seq3=inSequence[2:((leng-2)%3+2)] #broken
    seq2=inSequence[1:((leng-1)-(leng-1)%3+1)]
    seq3=inSequence[2:((leng-2)-(leng-2)%3+2)]
    return(seq1,seq2,seq3)

def isStop(inSequence):
    #end is TGA TAA or TAG
    if inSequence==["T","A","A"] or inSequence==["T","A","G"] or inSequence==["T","G","A"]:
        return True
    else:
        return False

def countORFs(inSequence):
    #start is ATG
    #end is TGA TAA or TAG
    outCount=0
    for trip in range(len(inSequence)):
        if trip%3==0 and inSequence[trip:trip+3]==["A","T","G"]:
            counter=0
            iterator=trip
            while iterator<len(inSequence):
                if iterator%3==0 and counter>500 and isStop(inSequence[iterator:iterator+3]):
                    outCount+=1
                iterator+=1
                counter+=1
    return outCount

def main():
    #seqF=openGeneFile("theFileame.txt")
    seqF=openGeneFile("HPV11alpha_GUMC-AJ.txt")
    seqR=giveReverseCompliment(seqF)
    seqF1,seqF2,seqF3=giveReadFrames(seqF)
    seqR1,seqR2,seqR3=giveReadFrames(seqR)
    print("An open reading frame is defined as an area within a genome that is downstrean of the 'start' signal of ATG, and upstream of a stop signal of TAA, TAG or TGA.")
    print("These areas are able to code for the mRNA, or messenger RNA, that gets translated into a protien. The larger the number of ORFs")
    print("a sequence contains, the larger the number of different protiens it could possibly code for. Each triplet of bases codes for a protien, therefore the given")
    print("sequence is analyzed 6 times: once for each forward reading frame and once for each reading frame on the other strand of the DNA, which is the reverse compliment of the first sequence.")
    #print(seqF1[:10]," ~~ ",seqF1[len(seqF1)-10:])
    print("Forward 1 ORFs: ",countORFs(seqF1))
    #print(seqF2[:10]," ~~ ",seqF2[len(seqF2)-10:])
    print("Forward 2 ORFs: ",countORFs(seqF2))
    #print(seqF3[:10]," ~~ ",seqF3[len(seqF3)-10:])
    print("Forward 3 ORFs: ",countORFs(seqF3))
    #print(seqR1[:10]," ~~ ",seqR1[len(seqR1)-10:])
    print("Reverse 1 ORFs: ",countORFs(seqR1))
    #print(seqR2[:10]," ~~ ",seqR2[len(seqR2)-10:])
    print("Reverse 2 ORFs: ",countORFs(seqR2))
    #print(seqR3[:10]," ~~ ",seqR3[len(seqR3)-10:])
    print("Reverse 3 ORFs: ",countORFs(seqR3))

main()
