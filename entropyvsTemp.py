#!/usr/bin/python

import sys
import numpy as np
import os
import subprocess
import re
import pexpect
import datetime
import divcalc9 as dc

#reads fasta file (>seq followed by newline followed by sequence followed by newline, start loop again)
#and outputs a list of names of structures and list of structures 

def read_fasta(infile):
    

    sequences = []
    with open(infile, 'r') as fastafile:
	lines = fastafile.readlines()
        for line in lines: 
            if line.startswith('>'):
	        seq_name = line.strip().strip('\n')
            elif len(line.strip()) > 0:
		sequences.append(line.strip())
 
    print seq_name

    return seq_name, sequences


import subprocess
def callRNAsuboptLinux(sequence, samples, temperature):
    """ shells to RNAsubopt to stochastically sample structures using the given input sequence and temperature """
    
    infile = "input.fas"
    
    writer = open(infile, "w")
    writer.write(">seq\n")
    writer.write(sequence+"\n")
    writer.close()    
    
    binaryPath = "//auto//dtchome//casablancasantrasv//Biostatistics//bin//RNAsubopt"
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


#T_list = [22,27,32,37,42,47,52]
num_sim = 10
num_samples = 100


t = float(sys.argv[1])*5+22
seq_index = 0

seqs = ['wt', 'A56U', 'C10U', 'C14G', 'U22G']
temp_ent = np.zeros((len(seqs),num_sim))

for seq in seqs:
    
    for sim in range(num_sim):
        print(sim)
        sequence = dc.getSequencesDict()[seq]
        print(sequence)
        samples = callRNAsuboptLinux(sequence, num_samples, t)
        structureSet = []
        for sample in samples:
            structureSet.append(["0","0","0",sample])
        entropy = dc.getEnt(sequence,structureSet)
        print("this is entropy:" + str(entropy))
        temp_ent[seq_index][sim] = entropy
        print(temp_ent[seq_index][sim])
    seq_index += 1
    print(seq_index/9*100)
np.save("%s" % t, temp_ent)


