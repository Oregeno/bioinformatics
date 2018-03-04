#!/usr/bin/python

import sys
import numpy as np
import os
import subprocess
import re
import pexpect
import datetime
#import matplotlib as mpl

#COPY THIS CODE OR SOMETHING THERABOUTS - replace A56U with the name of sequence of intrest
#rawStrucs = getStrucs()
#strucs = rawToStruc(rawStrucs)
#sequenceDict = getSequencesDict()
#structureSet1, structureSet2 = partitionStructures(strucs,"A56U","U22G")
#getEnt(sequenceDict['A56U'],structureSet1)






#get the list of dot comma structures for wild type and others
def getStrucs():
	current = os.path.dirname(os.path.abspath(__file__))
	stringlist = []

	for filename in os.listdir(current):
	    if filename.endswith(".txt"): 
		target = current + "//" + filename
		f = open(target,'r')
		structuresRaw = f.read()
	

		continue
	    else:
		continue

	return structuresRaw #tested

def getGeneratedStrucs():
	current = os.path.dirname(os.path.abspath(__file__))
	stringlist = []

	for filename in os.listdir(current):
	    if filename.endswith(".fas"): 
		target = current + "//" + filename
		f = open(target,'r')
		structuresRaw = f.read()
	

		continue
	    else:
		continue

	return structuresRaw #tested

#input z values for the files - maybe pull them from RNAfold
#Calulate Z values from FE's and probabilities

def rawToStruc(structuresRaw):

        strucData = re.findall('>(.*)_(.*)_(.*)\n(.*)',structuresRaw)

        return strucData  #tested


def getLogZ(sequence):

        child = pexpect.spawn('RNAfold -p')
        child.logfile = open("foldLog",'w+')
        child.expect('....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8')
        child.sendline(sequence)
        child.expect('')
        child.sendline('a')
        child.expect('....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8')
        string = child.before
        child.terminate()

        
        energy = re.findall('energy = ([-\s]?\S*) ',string)
        probability = re.findall(' in ensemble ([-\s]?\S*);',string)
        #print(energy[0])
        #print(probability[0])
        #print(-(float(energy[0]))-np.log(float(probability[0])))
        return float(energy[0]), float(probability[0]), -(float(energy[0]))-np.log(float(probability[0]))

#Feed non-wt sequences and wt dot comma structures into FE to get energies

def getEnergy(sequences,structures):
	
	current = os.path.dirname(os.path.abspath(__file__))
	args = sequences + " " + structures


	child = pexpect.spawn('RNAeval')
	child.setecho(True)
	child.logfile = open("mylog", "w+")
	child.expect('....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8')
	child.sendline(sequences)
	child.expect("")
	child.sendline(structures)
	child.expect("")

	child.sendline("a")
	


	#child.expect("use")
	child.expect('....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8')
	child.sendline(sequences)
	child.sendline(structures)
	#child.sendline("")


	child.terminate()

	#print("before fixing" + child.before)

        needsFixing, basesToFix = checkWarnings(child.before)

        if needsFixing == 1:
                newStructure = fixBases(basesToFix,structures)

                child = pexpect.spawn('RNAeval')
                child.setecho(True)
                child.logfile = open("mylog", "w+")
                child.expect('....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8')
                child.sendline(sequences)
                child.expect("")
                #print(newStructure)
                child.sendline(newStructure)
                child.expect("")

                child.sendline("a")
                


                #child.expect("use")
                child.expect('....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8')
                child.sendline(sequences)
                child.sendline(structures)
                #child.sendline("")


                child.terminate()

                #print('after fixing ' + child.before)
                
                res = re.findall("energy = ([-\s]?\S*) ", child.before)
                return res
                
	else:
                res = re.findall("energy = ([-\s]?\S*) ", child.before)
                return res

#convert energies to probabilities

def convertToLogP(energy,logz):

        logp = -energy - logz
        #print("logp" + str(logp))
        #print(energy)
        #print(logz)
        return logp
        
        
        


#calculate div pairwise between wt and altered sequences

#getEnergies("
#((.(((.(((((((((.(((((((((.(((...(((((......))))).))))))))))))))))))........)))...))).)).................(((.(((....((((((......)))))).....((((......))))..))).)))..........(((((.((........)).)))))...","GCAGTTCGGCGGTCCCGCGGGTCTGTCTCTTGCTTCAACAGTGTTTGGACGGAACAGATCCGGGGACTCTCTTCCAGCCTCCGACCGCCCTCCGATTTCCTCTCCGCTTGCAACCTCCGGGACCATCTTCTCGGCCATCTCCTGCTTCTGGGACCTGCCAGCACCGTTTTTGTGGTTAGCTCCTTCTTGCCAACCAACC")
#((.(((.(((((((((.((((((.((.(((...(((((......))))).))))).))))))))))))........)))...))).)).................(((.(((....((((((......)))))).....((((......))))..))).)))..........(((((.((........)).)))))...
def getSequencesDict():
        sequences = {
                'wt':  'GCAGTTCGGCGGTCCCGCGGGTCTGTCTCTTGCTTCAACAGTGTTTGGACGGAACAGATCCGGGGACTCTCTTCCAGCCTCCGACCGCCCTCCGATTTCCTCTCCGCTTGCAACCTCCGGGACCATCTTCTCGGCCATCTCCTGCTTCTGGGACCTGCCAGCACCGTTTTTGTGGTTAGCTCCTTCTTGCCAACCAACC',
                'wt2': 'GCAGTTCGGCGGTCCCGCGGGTCTGTCTCTTGCTTCAACAGTGTTTGGACGGAACAGATCCGGGGACTCTCTTCCAGCCTCCGACCGCCCTCCGATTTCCTCTCCGCTTGCAACCTCCGGGACCATCTTCTCGGCCATCTCCTGCTTCTGGGACCTGCCAGCACCGTTTTTGTGGTTAGCTCCTTCTTGCCAACCAACC',
                'A56U':'GCAGTTCGGCGGTCCCGCGGGTCTGTCTCTTGCTTCAACAGTGTTTGGACGGAACTGATCCGGGGACTCTCTTCCAGCCTCCGACCGCCCTCCGATTTCCTCTCCGCTTGCAACCTCCGGGACCATCTTCTCGGCCATCTCCTGCTTCTGGGACCTGCCAGCACCGTTTTTGTGGTTAGCTCCTTCTTGCCAACCAACC',
                'C10U':'GCAGTTCGGUGGTCCCGCGGGTCTGTCTCTTGCTTCAACAGTGTTTGGACGGAACAGATCCGGGGACTCTCTTCCAGCCTCCGACCGCCCTCCGATTTCCTCTCCGCTTGCAACCTCCGGGACCATCTTCTCGGCCATCTCCTGCTTCTGGGACCTGCCAGCACCGTTTTTGTGGTTAGCTCCTTCTTGCCAACCAACC',
                'C14G':'GCAGTTCGGCGGTGCCGCGGGTCTGTCTCTTGCTTCAACAGTGTTTGGACGGAACAGATCCGGGGACTCTCTTCCAGCCTCCGACCGCCCTCCGATTTCCTCTCCGCTTGCAACCTCCGGGACCATCTTCTCGGCCATCTCCTGCTTCTGGGACCTGCCAGCACCGTTTTTGTGGTTAGCTCCTTCTTGCCAACCAACC',
                'U22G':'GCAGTTCGGCGGTCCCGCGGGGCTGTCTCTTGCTTCAACAGTGTTTGGACGGAACAGATCCGGGGACTCTCTTCCAGCCTCCGACCGCCCTCCGATTTCCTCTCCGCTTGCAACCTCCGGGACCATCTTCTCGGCCATCTCCTGCTTCTGGGACCTGCCAGCACCGTTTTTGTGGTTAGCTCCTTCTTGCCAACCAACC'
                }
        return sequences




def getEnt(sequence,structureSet):

        sampleNumber = len(structureSet)

        probabilityList = []

        energy, probability, logZ = getLogZ(sequence)

        for structure in structureSet:
                energyRaw = getEnergy(sequence,structure[3])
                energy = float(re.findall("([-\s]?\S*)",energyRaw[0])[0])
                logProbability = convertToLogP(energy,logZ)

                probabilityList.append(logProbability)

        ent = np.sum(probabilityList)/sampleNumber
                

        return ent

def getDiv(sequence1,structureSet1,sequence2,structureSet2,sequenceDict):

        #print(sequence1)
        #print(sequence2)
        #print(structureSet1[:10])
        #print(structureSet2[:10])
        #initialise function parameters
        #print(structureSet1)
        div = 0

        #print(structureSet1[0][1])
        label1 = structureSet1[0][1]
        label2 = structureSet2[0][1]
        s2e_f_s1s = []
        s1e_f_s2s = []
        #print(label1)
        #print(label2)

        #get the log z's for turning energies into probabilities

        energy1, probability1, logZ1 = getLogZ(sequence1)
        energy2, probability2, logZ2 = getLogZ(sequence2)

        #print(energy1, probability1, logZ1)
        #print(energy2, probability2, logZ2)
        
        #get probabilities for sequence 2 on strucure set 1

        for data in structureSet1:
                #print(data)
                energyRaw = getEnergy(sequence2,data[3])
                energy = float(re.findall("([-\s]?\S*)",energyRaw[0])[0])
                logProbability = convertToLogP(energy,logZ2)

                #print("Sequence 2 for struc 1 E "+str(energy))
                #print("Sequence 2 for struc 1 P "+str(np.exp(logProbability)))
                

                se2st1 = logProbability

                energyRaw = getEnergy(sequence1,data[3])
                energy = float(re.findall("([-\s]?\S*)",energyRaw[0])[0])
                se1st1 = convertToLogP(energy,logZ1)
                #print(se1st1)

                #print("Sequence 1 for struc 1 E "+str(energy))
                #print("Sequence 1 for struc 1 P "+str(np.exp(se1st1)))
                
                s2e_f_s1s.append((se1st1-se2st1))

                
                #print(se1st1)
                #print(se2st1)

                #print(se2st1,se1st1,s2e_f_s1s)
        #get probabilities for sequnce 1 on structure set 2

        for data in structureSet2:


                #logprob of P(x)
                energyRaw = getEnergy(sequence1,data[3])
                energy = float(re.findall("([-\s]?\S*)",energyRaw[0])[0])
                logProbability = convertToLogP(energy,logZ1)
                se1st2 = logProbability

                #print("Sequence 1 for struc 2 P "+str(np.exp(se1st2)))


                energyRaw = getEnergy(sequence2,data[3])
                energy = float(re.findall("([-\s]?\S*)",energyRaw[0])[0])
                se2st2 = convertToLogP(energy,logZ2)

                #print("Sequence 2 for struc 2 P "+str(np.exp(se2st2)))
                

                s1e_f_s2s.append((np.exp(se1st2-se2st2))*(se1st2-se2st2))
                #print(se1st2,se2st2)
                #print("Importance Sampled"+str(s1e_f_s2s))
                #print(np.logaddexp(se1st2,-se2st2))
                #print(np.exp(se1st2-se2st2))
                #print(se1st2)
                #print(se2st2)
        
        #get mean of probabilities over probabilities

        #combine probabilities


        div = [np.sum(s2e_f_s1s), np.sum(s1e_f_s2s),np.sum(s2e_f_s1s)+np.sum(s1e_f_s2s)]
        #print(div)

        


        return div

def partitionStructures(strucs,label1,label2):

        structureSet1 = []
        structureSet2 = []

        for data in strucs:

                if data[1] == label1:

                        structureSet1.append(data)

                elif data[1] == label2:

                        structureSet2.append(data)
                
        return structureSet1, structureSet2



def generatePermutations(sequence):

        permutationList = [[],[],[],[]]
        numDict = {'T':0,
                   'G':1,
                   'C':2,
                   'A':3}

        sequenceSplit = list(sequence)

        for baseIndex in range(len(sequenceSplit)):
                #print(baseIndex)
                temp = sequenceSplit
                #print(temp)

                for base in ['A','T','G','C']:
                        
                        permutedBaseSequence = list(temp)
                        
                        if permutedBaseSequence[baseIndex] != base:
                                permutedBaseSequence[baseIndex] = base
                                #print(permutedBaseSequence)

                                permutationList[numDict[base]].append(permutedBaseSequence)
                        else:
                                permutationList[numDict[base]].append(permutedBaseSequence)

        #print(permutationList)

        #print(permutationList[0][0],permutationList[1][0],permutationList[2][0],permutationList[3][0])



        return permutationList
        


                                

        #print(permutationList)
#.......(((((((((.(((((((((.((.((.(((((......))))))))))))))))))))((......))........)))))))................(((.(((.(((((((((......))))))...............)))...))).)))..........(((((.((........)).)))))...
#.......(((((((((.((((((.((.((.((.(((((......))))))))))).))))))))((......))........)))))))................(((.(((.(((((((((......))))))...............)))...))).)))..........(((((.((........)).)))))...

#multiple broken
#.......(((((((((.(((((((((.((.((.(((((......))))))))))))))))))))((......))........)))))))................(((.(((.(((((((((......))))))...............)))...))).)))..........(((((.((........)).)))))...
#.......(((((((((.((((((.((.((....(((.(......).)))..)))).))))))))((......))........)))))))................(((.(((....((((((......)))))).....((((......))))..))).)))..........(((((.((........)).)))))...



#str1
#.......(((((((((.(((((((((.((.((.(((((......))))))))))))))))))))((......))........)))))))................(((.(((.(((((((((......))))))...............)))...))).)))..........(((((.((........)).)))))...
#seq2str1
#.......(((((((((.(..(.(.(......(.(..(........)..))....).).)..)))((......))........)))))))................(((.(((.(((((((((......))))))...............)))...))).)))..........(((((.((........)).)))))...

def checkWarnings(string):
        #print("bases To Fix")
        #print(string)
        basesToFixRaw = re.findall('bases (.*) and (.*) \(', string)

        #print(basesToFixRaw)
        if not basesToFixRaw:
                needsFixing = 0
                basesToFix = 0
        else:
                needsFixing = 1
                basesToFix = basesToFixRaw
        
        return needsFixing, basesToFix

def fixBases(basesToFix,structure):

        #print('made it to fix bases')

        structureArray = list(structure)

        for basePairs in basesToFix:

                for baseIndex in basePairs:

                        structureArray[int(baseIndex)-1] = "."

        newStructure = ''.join(structureArray)
        #print(newStructure)

        return newStructure


def saveSequences(sequenceArray):

        with open('permutedSequences.fas','w') as sequenceFile:

                for typeIndex, baseType in enumerate(sequenceArray):

                        for indexIndex, baseIndex in enumerate(baseType):
                                indexToBase = {'0':'A',
                                               '1':'T',
                                               '2':'G',
                                               '3':'C'}
                                
                                sequenceFile.write('>'+'wt'+'_'+str(indexIndex)+ '_' + str(indexToBase[str(typeIndex)]) + '\n')
                                sequenceFile.write(''.join(baseIndex) + '\n')

                sequenceFile.close()
                #print(indexIndex)

def makeSequencesDict(sequenceList):

        sequenceDict = {}

        for sequenceData in sequenceList:

                sequenceDict[sequenceData[1]] = sequenceData[3]
        
        return sequenceDict

def compareAll(sequenceListPath,wildTypePath,startingAt,endingAt):


        #pull the sequences out of list into array

        sequenceListRaw = getSequences('permutedSequences.fas')
        sequenceListUnsorted = fixSequences(sequenceListRaw)
        sequenceList = sorted(sequenceListUnsorted,key = lambda baseIndex: int(baseIndex[2]))
        #make a dictionary matching identities to sequences
        sequenceDict = makeSequencesDict(sequenceList)
        #print(sequenceList)
        #get ensemble for wild type

 
        filePath = str(sequenceListPath) + '//' + str(wildTypePath)
        wildTypeEnsemble = getEnsemble(filePath)
#issue is that the sequences not fixed properly
#GCAGTTCGGCGGTCCCGCGGGTCTGTCTCTTGCTTCAACAGTGTTTGGACTGAACAGATCCGGGGACTCTCTTCCAGCCTCCGACCGCCCTCCGATTTCCTCTCCGCTTGCAACCTCCGGGACCATCTTCTCGGCCATCTCCTGCTTCTGGGACCTGCCAGCACCGTTTTTGTGGTTAGCTCCTTCTTGCCAACCAACC
        #sequenceFile = open('permutedSequences.fas','r')

        #sequences = sequenceFile.read()

        #structureArray = rawToStruc(sequences)

        #indexSequenceDict = {}

        #for indexSequencePair in structureArray:

                #indexSequenceDict[indexSequencePair[2]] = indexSequencePair[3]

                #print(indexSequenceDict)
        startTime = str(datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y"))
        
        for alternateSequence in sequenceList[startingAt:endingAt]:

                with open('calculatedDivs3//CalculatedEnergies_'+str(startingAt)+'_'+startTime,'a+') as energyFile:

                
#GCAGTTCGGCGGTCCCGCGGGTCTGTCTCTTGCTTCAACAGTGTTTGGACTGAACAGATCCGGGGACTCTCTTCCAGCCTCCGACCGCCCTCCGATTTCCTCTCCGCTTGCAACCTCCGGGACCATCTTCTCGGCCATCTCCTGCTTCTGGGACCTGCCAGCACCGTTTTTGTGGTTAGCTCCTTCTTGCCAACCAACC
                        #print(alternateSequence)

                        filePath = str(sequenceListPath) + '//' + str('>' + alternateSequence[0]+'_'+alternateSequence[1] + '.txt')
                        alternateEnsemble = getEnsemble(filePath)
                        filePath = str(sequenceListPath) + '//' + str('>' + alternateSequence[0]+'_'+alternateSequence[1] + '.fas')
                        aSequence = getAlternateSequence(filePath)
                        #print(aSequence[0])
                        #print(filePath)
                        #print(alternateSequence)
			#print(sequenceDict[alternateSequence[1]])
			#print(alternateEnsemble)
                        div = getDiv(sequenceDict['0_T'],wildTypeEnsemble,aSequence[0],alternateEnsemble,sequenceDict)
                        print('#################################' + str(div))
                        energyFile.write('wt'+'_'+ alternateSequence[1] + str(div) + '\n')

                        
                        
                        

def fixSequences(sequenceArray):

        fixedSequences = []

        for sequenceData in sequenceArray:

                #print(sequenceData)

                fixedSequences.append([sequenceData[0],sequenceData[1]+'_'+sequenceData[2],sequenceData[1],sequenceData[3]])

        return fixedSequences

def getSequences(sequencesRaw):

        with open(sequencesRaw,'r') as sequenceFile:

                sequenceRead = sequenceFile.read()

                sequences = rawToStruc(sequenceRead)

                #print(sequences)

        return sequences

def getEnsemble(filePath):
        #print(filePath)
        ensembleFile = open(filePath,'r')
        
        ensembleRead = ensembleFile.read()

        #print(ensembleRead)

        ensembleArray= rawToStruc(ensembleRead)

        return ensembleArray

def getAlternateSequence(filePath):

        ensembleFile = open(filePath,'r')
        
        ensembleRead = ensembleFile.read()

        #print(ensembleRead)

        ensembleArray = re.findall('>.*\n(.*)',ensembleRead)

        return ensembleArray

def plotMap(directory):

        #collect all the files and raed them to strings

        mapStrings = getMapStrings(directory)

        #read strings to list

        mapList = makeMapList(mapStrings)

        #use id's to place values in array

        mapArray = makeMapArray(mapArray)

        #show image

        maplotlib.imshow(mapArray)

        

        #plot with imshow
        

def getMapStrings(directory):

	mapStringList = []

	for filename in os.listdir(directory):
	    if filename.endswith(".txt"): 
		target = directory + "//" + filename
		f = open(target,'r')
		mapStringList.append(f.read())
	

		continue
	    else:
		continue

	return mapStringList #tested

 #       with open(filePath,'r') as 

#actual code starts

rawStrucs = getStrucs()
strucs = rawToStruc(rawStrucs)
sequenceDict = getSequencesDict()
structureSet1, structureSet2 = partitionStructures(strucs,"A56U","U22G")
#print(getLogZ(sequenceDict['wt']),getLogZ(sequenceDict['A56U']),getLogZ(sequenceDict['C14G']),getLogZ(sequenceDict['C10U']),getLogZ(sequenceDict['U22G']))


#print(np.sum(getDiv(structureSet1, structureSet2, sequenceDict)))

#sequenceArray = generatePermutations(sequenceDict['wt'])

#saveSequences(sequenceArray)

#generatedStrucs = getGeneratedStrucs()

compareAll('all_ensembles','>wt_0_T.txt',int(sys.argv[1]),int(sys.argv[2]))

#A>U/T
#G>C
#A56U C10U C14G U22G
#3.45795307979
#34.5527948305

#34.6651583289
#348.135595478
