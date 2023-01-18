import dna

def load(filename) :
    """Takes the filename of a .fna file and returns the DNA sequence and description line"""
    with open(filename, 'r') as file_read:
        if filename.endswith(".fna") :
            for line in file_read:
                line = line.strip("\n")
                if line.startswith(">"):
                    description = str(line.strip(">"))
                else:
                    sequence = [acid.upper() for acid in line]
        else :
            sequence = list()
            description = str("file was not loaded")
    return(sequence, description)

def stats(sequence) :
    """Take a sequence and return a table that includes the number of times a nucleic acid code occurs."""
    table_k =('A', 'C', 'G', 'T', 'N', 'U', 'K', 'S', 'Y', 'M', 'W', 'R', 'B', 'D', 'H', 'V', '-')
    other = list(set(sequence).difference(set(table_k)))
    table = dict.fromkeys(table_k)
    for k_values in table.keys() :
        table[k_values] = sequence.count(k_values)
    counter = 0
    for val in other:
        counter += sequence.count(val)
    table["Other"] = counter   
    return (table)

def format_sequence(sequence, first_index, last_index) :
    """Take a sequence along with two indices and return the subsequence with a particular format"""
    extract = sequence [first_index : last_index + 1]
    new_extract = list()
    for i in range (0, len(extract), 80) :
        up_extract = extract [i : i + 80]
        new_extract.append(''.join(up_extract))
    return (new_extract)

def write(filename, description, sequence, first_index, last_index):
    """Take a description, sequence, and sequence range, and write to a .fna file"""
    file_w = open (filename, "w")
    file_w.write(">" + str(description) + "\n")
    new_file = [file_w.write(i) for i in sequence [first_index : last_index + 1] ] 
        

def dna2protein(dna_sequence):
    """Convert the DNA sequence to its corresponding protein sequence"""
    with open("dna2protein.csv", 'r') as file:
        line = file.readline()
        data = list()
        while line:
            data.append(line.strip().split(','))
            line = file.readline()
    acid = list()
    for i in range(0,len(data)):
        codes = (data[i][0] + data[i][1] + data[i][2])
        acid.append(codes)    
    table = dict.fromkeys(acid)
    j = 0
    for k_values in table.keys() :
        table[k_values] = data[j][4]
        j += 1
    table.update ({"???" : "?"})
    default = table.setdefault("???", "?")
    sequence = list()
    for i in range (0, len(dna_sequence), 3) :
        extract = dna_sequence [i : i + 3]
        sequence.append(''.join(extract))
    protein_sequence = list()
    for i in range(0, len(sequence)) :
        if sequence[i] in table.keys() :
            protein_sequence.append(table.get(sequence[i]))
        else :
            protein_sequence.append(default)
    return (protein_sequence, table)


import time



#Part 1 - [6.9]
seq, info = dna.load('chr22.fna')
##print(info)
##print(seq[0:80])

t1 = time.time()

#Part 2 - [15.8]
table = dna.stats(seq)
print(table)
##
###Part 3 - [0.049]
##data = dna.format_sequence(seq,100,300)
##print(data)

###Part 4 - [0.015]
##dna.write('BRCA1_subseq.fna', 'fragment sequence of BRCA1.fna', seq,100,300)

###Part 5 - [0.78]
##seq = list('AAAGTTAAATAATAAATAGGTGAA')
##seq_to_find = list('AAA')
##matches = dna.find(seq,seq_to_find)
##print(matches)
##
###Part 6, 7, 8 - [0.05]
##line = ''
##seq = list('AAAGTTAAATAATAAATAGGTGAA')
##print(line.join(seq))
##seq = dna.add(seq,list('NNNNN'),5)
##print(line.join(seq))
##seq = dna.delete(seq,6,4)
##print(line.join(seq))
##seq = dna.replace(seq,list('HHH'),5, 1)
##print(line.join(seq))
##
###Part 9 - [0.039]
##dna_seq = list('AAAGTTAAATAATAAATAGGTGAA')
##pro_seq, table = dna.dna2protein(dna_seq)
##print(pro_seq)
##print(table)

t2 = time.time()

t3 = time.time()

seq2, info2 = load('chr22.fna')
##print(info2)
##print(seq2[0:80])

table = stats(seq2)
print(table)

t4 = time.time()

print(t2-t1)
print(t4-t3)
