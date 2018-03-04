# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 14:17:35 2015

@author: Michael
"""
import numpy as np
import os

def dotBracketToPairedSites(dotBracket):
    """ convert a dot bracket string to a paired sites list """    
    
    pairedSites = np.zeros(len(dotBracket), dtype='int32')    
    stack = []
    for j in xrange(len(dotBracket)):
        if dotBracket[j] == "(":
            stack.append(j)
        elif dotBracket[j] == ")":
            x = stack.pop()
            y = j
            pairedSites[x] = y + 1
            pairedSites[y] = x + 1
    
    return pairedSites

def nChars(c, length):
    """ string of c repeated length times"""
    return c*length
    
def getStructureStar(length):
    """ required for computing the mountain metric diameter """
    if length % 2 == 0:
        return dotBracketToPairedSites(nChars('(', (length - 2) / 2) + ".." + nChars(')', (length - 2) / 2))
    else:
        return dotBracketToPairedSites(nChars('(', (length - 1) / 2) + "." + nChars(')', (length - 1) / 2))

def getStructureZero(length):
    """ required for computing the mountain metric diameter """
    return dotBracketToPairedSites("."*length)    


def getMountainVector(pairedSites, weighted):
    if weighted:
        v = np.zeros(len(pairedSites))
        for i in xrange(len(pairedSites)):
            v[i] = v[i-1]
            if pairedSites[i] != 0:
                if pairedSites[i] > i:
                    v[i] += 1.0/(pairedSites[i]-1.0-i)
                else:
                    v[i] -= 1.0/(i-(pairedSites[i]-1.0))
    else:        
        v = np.zeros(len(pairedSites), 'int32')
        for i in xrange(len(pairedSites)):
            v[i] = v[i-1]
            if pairedSites[i] != 0:
                if pairedSites[i] > i:
                    v[i] += 1
                else:
                    v[i] -= 1
    return v

def getMountainDistance(mountainVector1, mountainVector2, power):
    d = 0.0
    for i in xrange(len(mountainVector1)):
        d += np.power(np.abs(mountainVector1[i]-mountainVector2[i]),power)
    return np.power(d, 1.0/power)

def getMountainDiameter(length, power, weighted):
    """ required to compute the normalized mountain distance, the mountain diameter is the maximum possible distance between two structures of a given length """
    return getMountainDistance(getMountainVector(getStructureStar(length), weighted),getMountainVector(getStructureZero(length), weighted),power)

def calculateNormalizedMountainDistance(pairedSites1, pairedSites2, power, weighted):
    mv1 = getMountainVector(pairedSites1,weighted)
    mv2 = getMountainVector(pairedSites2,weighted)
    return  getMountainDistance(mv1, mv2, power) / getMountainDiameter(len(pairedSites1), power, weighted)
    
import subprocess
def callRNAsuboptLinux(sequence_name, sequence, samples, temperature):
    """ shells to RNAsubopt to stochastically sample structures using the given input sequence and temperature """
    
    infile = "%s.fas" % sequence_name
    
    writer = open(infile, "w")
    writer.write(">seq\n")
    writer.write(sequence+"\n")
    writer.close()    
    
    binaryPath = "/auto/dtchome/casablancasantrasv/Biostatistics/bin/RNAsubopt"
    args = "%s --temp=%f --stochBT=%d < \"%s\"" % (binaryPath, temperature, samples, infile)
    print args
    #args = "%s --temp=%f --stochBT=%d < \"%s\" > %s" % (binaryPath, temperature, samples, infile, outfile)
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
    out, err = proc.communicate()

    print err
    lines = out.split("\n")

    sample = []

    for line in lines:
        line2 = line.strip()
        #print line2
        if not line2.startswith(">") and len(line2) == len(sequence.strip()):
            sample.append(line2)
            #sample.append(dotBracketToPairedSites(line))
    proc.wait()
    #print sample
    return sample # returns a list of structures in dotbracket format
 

s1 = dotBracketToPairedSites("...(..).......")
s2 = dotBracketToPairedSites("..............")
#print(calculateNormalizedMountainDistance(s1, s2, 1.0, True))


#samples = callRNAsuboptWindows("ACGAGCGCGATGATGCTAGCTGATCGTAACGCTT", length(num_samples), 37.0)
#print samples
#print(len(samples))



    
def read_fasta(infile):
    
    seq_names = []
    sequences = []
    with open(infile, 'r') as fastafile:
	lines = fastafile.readlines()
        for line in lines: 
            if line.startswith('>'):
	        
                seq_names.append(line.strip().strip('\n'))
            elif len(line.strip()) > 0:
		sequences.append(line.strip())
 
    print seq_names

    return seq_names, sequences


#generate text files and mountain vector matrix for RNA sequences
def generate_ensembles(filename, num_samples): 
    
    
    names, sequences = read_fasta(filename) 
    
    i = 0
    for name in names: 
        #generates ensemble in dot bracket notation
        samples = callRNAsuboptLinux(name, sequences[i], num_samples, 37.0) 
        # appends ensemble in dot bracket notation to a list
	
	with open(name + '.txt', 'w') as output_file:
	    j = 0 
	    for sample in samples: 
		output_file.write(name.split('_')[0] + '_' + name.split('_')[1] + '_%d\n' % j)
                output_file.write(sample + '\n')
                j += 1 
        i += 1
    
	


    

file_names = generate_ensembles("permutedSequences.fas", 100)
            



