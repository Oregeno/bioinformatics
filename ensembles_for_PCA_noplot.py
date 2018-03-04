#!/usr/bin/python

"""
Created on Wed Mar 11 14:17:35 2015

@author: Michael
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import csv


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
    #print args
    #args = "%s --temp=%f --stochBT=%d < \"%s\" > %s" % (binaryPath, temperature, samples, infile, outfile)
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
    out, err = proc.communicate()

    #print err
    lines = out.split("\n")

    sample = []

    for line in lines:
        line2 = line.strip()
        #print line2
        if not line2.startswith(">") and len(line2) == len(sequence.strip()):
            sample.append(line2)
            
    proc.wait()
    #print sample
    return sample # returns a list of structures in dotbracket format
 




#generate text files and mountain vector matrix for RNA sequences
def generate_ensembles(folder, num_samples): 
    
    files = os.listdir(folder)

    samples_list = []
    file_names = []
    with open("dot_bracket.txt", 'w') as dotbracket_textfile:

        for f in files:     
	    if f.endswith(".fas"): 
                file_names.append(f.split('.')[0])       
                input_file = open(folder+f, 'r')
                lines = input_file.readlines()
	        for l in lines:
                    if (not l.startswith(">")):
	                sequence = l.strip("\n")
                #generates ensemble in dot bracket notation
                samples = callRNAsuboptLinux(f.split('.')[0], sequence, num_samples, 37.0) 
    	        # appends ensemble in dot bracket notation to a list
		samples_list.append(samples)

	rna_length = len(samples_list[0][0])

	mountains = np.zeros((num_samples*len(file_names), rna_length))
	pairedsites = np.zeros((num_samples*len(file_names), rna_length))
	pairedsites_original = np.zeros((num_samples*len(file_names), rna_length))
	
	for i in range(0, len(samples_list)):
			
		for j in range(0, len(samples_list[0])):
			dotbracket_textfile.write('>' + str(file_names[i]) + '_' + str(j) + '\n')
			dotbracket_textfile.write(str(samples_list[i][j] + '\n'))
		        v = getMountainVector(dotBracketToPairedSites(samples_list[i][j]), True)
                        #print len(v)
			
			mountains[i*num_samples+j] = v

                        ps = dotBracketToPairedSites(samples_list[i][j])
			pairedsites_original[i*num_samples+j] = np.array(ps)
			pairedsites[i*num_samples+j] = np.array([min(e,1) for e in ps])


      

        np.save("pairedSites", pairedsites)
        np.save("mountains", mountains)
        #print pairedsites_original[0]
	


    return file_names

file_names = generate_ensembles("./Data/", 1000)
            

num_samples = 1000

def PCA(file_names, metric, plot=False):
 

    # Load in matrix
    if metric=='mountain':
        matrix = np.load('mountains.npy')
    elif metric=='paired':
        matrix = np.load('pairedSites.npy')
        
    #print matrix[0]

    # calculate the covariance matrix and from there the eigenvalues/vectorsd
    cov_mat = np.cov(matrix.T)
    eig_val_cov, eig_vec_cov = np.linalg.eig(cov_mat)

    # make some tuples of eigenvalues and vectors
    eig_pairs = [(np.abs(eig_val_cov[i]), eig_vec_cov[:,i]) for i in range(len(eig_val_cov))]

    # Sort the tuples from high to low
    eig_pairs.sort(key=lambda x: x[0], reverse=True)

    #calculate the contribution of the two top vectors and store in variable top2
    eigvals = [v[0] for v in eig_pairs]
    normeigvals = eigvals / np.sum(eigvals) * 100
    top2 = normeigvals[0:2]

    # choose the two vectors with the highest eigenvalues
    eigen_matrix = np.hstack((eig_pairs[0][1].reshape(matrix.shape[1],1), eig_pairs[1][1].reshape(matrix.shape[1],1)))

    #transform the original data
    transformed = eigen_matrix.T.dot(matrix.T)


    if plot==True: 

        color_list = ['blue','red', 'green', 'black', 'orange', 'magenta']

        figure1 = plt.figure()
        for i in range(0, 1000):

            if i == 1:
                plt.plot(transformed[1,i],transformed[0,i], 'o', markersize=7, color='%s' % color_list[0], alpha=0.3, label = "%s" % file_names[0])
            else:
                plt.plot(transformed[1,i],transformed[0,i], 'o', markersize=7, color='%s' % color_list[0], alpha=0.3)
                
        for i in range(1000, 2000):
            if i == 1001:
                plt.plot(transformed[1,i],transformed[0,i], 'x', markersize=7, color='%s' % color_list[1], alpha=0.3, label = "%s" % file_names[1])
            else:
                plt.plot(transformed[1,i],transformed[0,i], 'x', markersize=7, color='%s' % color_list[1], alpha=0.3)

        for i in range(2000, 3000):
            if i == 2001:
                plt.plot(transformed[1,i],transformed[0,i], 'o', markersize=7, color='%s' % color_list[2], alpha=0.3, label = "%s" % file_names[2])
            else:
                plt.plot(transformed[1,i],transformed[0,i], 'o', markersize=7, color='%s' % color_list[2], alpha=0.3)
                
        for i in range(3000, 4000):
            if i == 3001:
                plt.plot(transformed[1,i],transformed[0,i], 'o', markersize=7, color='%s' % color_list[3], alpha=0.3, label = "%s" % file_names[3])
            else:
                plt.plot(transformed[1,i],transformed[0,i], 'o', markersize=7, color='%s' % color_list[3], alpha=0.3)

        for i in range(4000, 5000):
            if i == 4001:
                plt.plot(transformed[1,i],transformed[0,i], 'o', markersize=7, color='%s' % color_list[4], alpha=0.3, label = "%s" % file_names[4])
            else:
                plt.plot(transformed[1,i],transformed[0,i], 'o', markersize=7, color='%s' % color_list[4], alpha=0.3)


        for i in range(5000, 6000):
            if i == 5001:
                plt.plot(transformed[1,i],transformed[0,i], 'o', markersize=7, color='%s' % color_list[5], alpha=0.3, label = "%s" % file_names[5])
            else:
                plt.plot(transformed[1,i],transformed[0,i], 'o', markersize=7, color='%s' % color_list[5], alpha=0.3)


        plt.xlabel('eigenvector 1')
        plt.ylabel('eigenvector 2')
        plt.legend(fontsize = 'small')

        if metric=='mountain':
            plt.xlim(-3,3)
            plt.ylim(-3,3)
            plt.title('RNA ensemble representation in 2D space (mountain metric)')
            plt.savefig("mountain_all.png")
        elif metric=='paired':
            plt.xlim(-4,5)
            plt.xlim(-4,5)
            plt.title('RNA ensemble representation in 2D space (paired sites metric)')
            plt.savefig("paired_all.png")
            
        n=0
        for name in file_names: 
            fig = plt.figure()
            plt.xlabel('eigenvector 1')
            plt.ylabel('eigenvector 2')
            plt.legend(fontsize = 'small')
            #plt.xlim(-3,3)
            #plt.ylim(-3,3)
            for i in range(n*1000, (n+1)*1000):
                plt.plot(transformed[1,i],transformed[0,i], 'o', markersize=7, color= '%s' % color_list[n], alpha=0.3)

            if metric=='mountain':
                plt.title('RNA ensemble representation of ' + name + ' in 2D space (mountain metric)')
                plt.savefig("mountain_%s.png" % name)
            
            elif metric=='paired':
                plt.title('RNA ensemble representation of ' + name + ' in 2D space (paired sites metric)')
                plt.savefig("paired_%s.png" % name)

            n+=1
            
    return top2

    
# actual code starts 


num_samples = 1000
num_simulations = 100

with open("topvectors.csv",'w') as vectorfile:
    vectorf = csv.writer(vectorfile)
    paired_top2 = np.zeros((1,num_simulations))
    mountain_top2 = np.zeros((1,num_simulations))

    for sim in range(num_simulations):
        print("Sim = %d" % sim)
        file_names = generate_ensembles("./Data/", 1000)
        topvectors = []
        if sim==0: #generate only one representative set of plots
            pairedtop2 = PCA(file_names, plot=True,metric='paired')
            print(pairedtop2)
            topvectors.extend(pairedtop2)
            paired_top2[0][sim] = sum(pairedtop2)

            mountaintop2 = PCA(file_names,plot=True,metric='mountain')
            print(mountaintop2)
            topvectors.extend(mountaintop2)
            mountain_top2[0][sim] = sum(mountaintop2)

            #in the csv file the first two columns are for the paired metric, the 3rd and 4th are for the mountain metric           
            vectorf.writerow(topvectors)
        else: 
            pairedtop2 = PCA(file_names, plot=False,metric='paired')
            topvectors.extend(pairedtop2)
            paired_top2[0][sim] = sum(pairedtop2)

            mountaintop2 = PCA(file_names,plot=False,metric='mountain')
            topvectors.extend(mountaintop2)
            mountain_top2[0][sim] = sum(mountaintop2)

            #in the csv file the first two columns are for the paired metric, the 3rd and 4th are for the mountain metric           
            vectorf.writerow(topvectors)

    np.save("mountain_top2", mountain_top2)
    np.save("paired_top2", paired_top2)




