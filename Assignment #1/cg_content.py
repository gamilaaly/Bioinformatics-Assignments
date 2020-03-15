sequenceFile = open('DNAsequence.txt', 'r')
sequence = sequenceFile.read()

G_bases_count = sequence.count('G')
C_bases_count = sequence.count('C')

CG_content_percentage = (G_bases_count + C_bases_count) / len(sequence) * 100

print("The number of cytosine bases = " + str(C_bases_count) + "\n" + "The number of guanine bases =  " + str(G_bases_count) + "\n" +"The CG content percentage = " + str(CG_content_percentage) + " %")

