from math import floor

# Open the file which contains the sequence, and store the sequence as a string in a variable to be possessed
sequenceFile = open('DNAsequence.txt', 'r')
sequence = sequenceFile.read()

# Each 3 bases in the sequence can make an amino acid
possible_amino_acid_number = floor(len(sequence) / 3)

# So, the spare bases in the sequence will be the remaining after all possible amino acids are constructed
spare_bases = len(sequence) % 3

print("DNA sequence : " + sequence + "\n" +
      "Possible number of amino acids resulting from this sequence = " + str(possible_amino_acid_number) + "\n" +
      "Number of bases that do not make a whole amino acid = " + str(spare_bases)
      )
