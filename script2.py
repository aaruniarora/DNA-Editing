import dna

import time

t1 = time.time()

#Part 1 - [0.0467]
print ("Part 1")
seq, info = dna.load('BRCA1.fna')
seq2, info2 = dna.load('Trial.txt')
seq3, info3 = dna.load('BRCA2.fna')
seq4, info4 = dna.load('BRCA4.fna')
print(info)
print(seq[0:80])
print("\n")
print(info2)
print(seq2[0:80])
print("\n")
print(info3)
print(seq3[0:80])
print("\n")
print(info4)
print(seq4[0:80])
print("\n")

#Part 2 - [0.0842]
print ("Part 2")
seq10 = list("atgcatgcatgc!")
table = dna.stats(seq)
table3 = dna.stats(seq3)
table4 = dna.stats(seq4)
table10 = dna.stats(seq10)
print(table)
print(table3)
print(table4)
print(table10)
print("\n")

#Part 3
print ("Part 3")
data = dna.format_sequence(seq,80950,81070)
data = dna.format_sequence(seq10,0,10)
print(data)
print("\n")

#Part 4
print ("Part 4")
dna.write('BRCA1_subseq.fna', 'fragment sequence of BRCA1.fna', seq, 0, 1000)
print("\n")

#Part 5
print ("Part 5")
seq_to_find = list('ATTCG')
matches = dna.find(seq,seq_to_find)
print(matches)
print("\n")

#Part 6, 7, 8
print ("Part 6, 7, 8")
line = '';
seq = list('AAAGTTAAATAATAAATAGGTGAA')
print(line.join(seq))
seq = dna.add(seq,list('NNNNN'),5) #should be placed before the index, otherwise end
print(line.join(seq))
seq = dna.delete(seq,28,4) #delete from index to number/until the end/out of bounds > nothing
print(line.join(seq))
seq = dna.replace(seq,list('HHH'),28, 0) #index > sequence to delete from, seq_to_add before the index/ same as add&del
print(line.join(seq))
print("\n")

#Part 9
print ("Part 9")
dna_seq = list('aaa-TTAAATAATAAATAGGTGAAK')
pro_seq, table = dna.dna2protein(dna_seq)
print(pro_seq)
print(table)
print("\n")

t2 = time.time()
print(t2-t1)
