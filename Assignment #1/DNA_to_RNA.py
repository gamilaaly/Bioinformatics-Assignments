# Open the file which contains the sequence, and store the sequence as a string in a variable to be possessed
sequenceFile = open('DNAsequence.txt', 'r')
DNA_sequence = sequenceFile.read()

# Conversion function from RNA to DNA that replace each Thymine base to Uracil base
# replace is a built in function in python which takes the char needed to be replaced as first argument and the char that will be replaced with as a second argument
RNA_sequence = DNA_sequence.replace('T', 'U')

print("DNA Sequence : " + DNA_sequence + "\n"
      "RNA Sequence : " + RNA_sequence)
