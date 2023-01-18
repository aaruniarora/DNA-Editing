def load(filename) :
    """Takes the filename of a .fna file and returns the DNA sequence and description line"""
    with open(filename, 'r') as file_read:
        sequence = list()
        if filename.endswith(".fna") :
            for line in file_read:
                line.rstrip()
                line = line.strip("\n")
                if line.startswith(">"):
                    description = str(line.strip(">"))
                else:
                    for g in line:
                        g = g.upper()
                        sequence.append(g)
        else :
            sequence = sequence
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
    for p in other:
        counter += sequence.count(p)
    table["Other"] = counter
    return (table)

def format_sequence(sequence, first_index, last_index) :
    """Take a sequence along with two indices and return the subsequence with a particular format"""
    extract = sequence [first_index : last_index+1]
    new_extract = list()
    for i in range (0, len(extract), 80) :
        up_extract = extract [i:i+80]
        new_extract.append(''.join(up_extract))
    return (new_extract)

def write(filename, description, sequence, first_index, last_index):
    """Take a description, sequence, and sequence range, and write to a .fna file"""
    file_w = open (filename, "w")
    file_w.write(">" + str(description)+"\n")
    for i in sequence [first_index:last_index+1] : 
        file_w.write(i)

def find(sequence, sequence_to_find) : 
    """Find a sequence within another sequence and record the indices where they occurred."""
    matches = list()
    for i in range (0, len(sequence)-1) :
        if sequence[i] == sequence_to_find [0] :
            if sequence [i:i+len(sequence_to_find)] == sequence_to_find [0:len(sequence_to_find)] :
                matches.append(i)
    return (matches)            

def add(sequence, sequence_to_add, index): 
    """Add a sequence into an existing sequence at a specified index"""
    new_seq = str()
    if 0 < int(index) < len(sequence) :
        new_seq = sequence [0:index] + sequence_to_add + sequence [index:len(sequence)]
    else :
        new_seq = sequence + sequence_to_add
    return (new_seq)

def delete(sequence, index, number_of_codes):
    """Delete a subsequence from a sequence as specified by a starting index and the number of codes to delete"""
    new_seq = str()
    if 0 < int(index + number_of_codes) < len(sequence) :
        new_seq = sequence [0:index] + sequence [index + number_of_codes:len(sequence)]
    else :
        new_seq = sequence [0:index]
    return (new_seq)

def replace(seq, seq_to_add, ind, num_codes): 
    """Replace a section of a sequence with a new subsequence"""
    if ((int(ind)+int(num_codes)) < len(seq)):
        seq[ind:ind+num_codes] = seq_to_add
        new_seq = ''.join(seq)
    else :
        new_seq = str()
        new_seq = seq [:] + seq_to_join [:]
    return(new_seq)
    
def dna2protein(dna_sequence): 
    """Convert the DNA sequence to its corresponding protein sequence"""
    with open("dna2protein.csv", 'r') as file:
        line = file.readline()
        data = list()
        while line:
            data.append(line.strip().split(','))
            line = file.readline()
    var1 = list()
    for i in range(0,len(data)):
        codes = (data[i][0] + data[i][1] + data[i][2])
        var1.append(codes)    
    table = dict.fromkeys(var1)
    j = 0
    for k_values in table.keys() :
        table[k_values] = data[j][4]
        j = j + 1
    table.update ({"???" : "?"})  
    sequence = list()
    for i in range (0, len(dna_sequence), 3) :
        extract = dna_sequence [i:i+3]
        sequence.append(''.join(extract))
    protein_sequence = list()
    for i in range(0, len(sequence)) :
        if sequence[i] in table.keys() :
            protein_sequence.append(table.get(sequence[i]))
    return (protein_sequence, table)
