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
