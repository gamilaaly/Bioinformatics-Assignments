sequenceFile = open('DNAsequence.txt', 'r')
DNA_sequence = sequenceFile.read()

RNA_sequence = DNA_sequence.replace('T','U')

print("DNA Sequence : " + DNA_sequence + "\n"
      "RNA Sequence : " + RNA_sequence)