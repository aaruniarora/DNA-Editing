import dna

import time

start1 = time.time()

#Part 1
print ("Part 1")
seq, info = dna.load('chr22.fna')
print(info)
print(seq[0:80])
print("\n")

end1 = time.time()
start2 = time.time()

#Part 2
print ("Part 2")
table = dna.stats(seq)
print(table)
print("\n")

end2 = time.time()
start3 = time.time()

#Part 3
print ("Part 3")
data = dna.format_sequence(seq,100,300)
print(data)
print("\n")

end3 = time.time()
start4 = time.time()

#Part 4
print ("Part 4")
dna.write('BRCA1_subseq.fna', 'fragment sequence of BRCA1.fna', seq,100,300)
print("\n")
end4 = time.time()

start5 = time.time()

#Part 5
print ("Part 5")
seq = list('AAAGTTAAATAATAAATAGGTGAA')
seq_to_find = list('AAA')
matches = dna.find(seq,seq_to_find)
print(matches)
print("\n")

end5 = time.time()

#Part 6, 7, 8
print ("Part 6, 7, 8")
line = '';
seq = list('AAAGTTAAATAATAAATAGGTGAA')
print(line.join(seq))

#Part 6
start6 = time.time()
seq = dna.add(seq,list('NNNNN'),5)
print(line.join(seq))
end6 = time.time()

#Part 7
start7 = time.time()
seq = dna.delete(seq,6,4)
print(line.join(seq))
end7 = time.time()

#Part 8
start8 = time.time()
seq = dna.replace(seq,list('HHH'),5, 1)
print(line.join(seq))
end8 = time.time()

print("\n")

start9 = time.time()

#Part 9
print ("Part 9")
dna_seq = list('AAAGTTAAATAATAAATAGGTGAA')
pro_seq, table = dna.dna2protein(dna_seq)
print(pro_seq)
print(table)
print("\n")

end9 = time.time()

print("1 : " + str(end1-start1))
print("2 : " + str(end2-start2))
print("3 : " + str(end3-start3))
print("4 : " + str(end4-start4))
print("5 : " + str(end5-start5))
print("6 : " + str(end6-start6))
print("7 : " + str(end7-start7))
print("8 : " + str(end8-start8))
print("9 : " + str(end9-start9))
