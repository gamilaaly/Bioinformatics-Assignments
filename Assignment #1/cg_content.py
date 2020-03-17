# Open the file which contains the sequence, and store the sequence as a string in a variable to be possessed
sequenceFile = open('DNAsequence.txt', 'r')
sequence = sequenceFile.read()

# The count function takes the char needed to be counted as an argument and returns its count in a given string
# These functions gets the number of cytosine and guanine bases in a given DNA
G_bases_count = sequence.count('G')
C_bases_count = sequence.count('C')

# 
CG_content_percentage = (G_bases_count + C_bases_count) / len(sequence) * 100

print("The number of cytosine bases = " + str(C_bases_count) + "\n"
      "The number of guanine bases =  " + str(G_bases_count) + "\n"
      "The CG content percentage = " + str(CG_content_percentage) + " %")
