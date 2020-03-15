from math import  floor
sequenceFile = open('DNAsequence.txt', 'r')
sequence = sequenceFile.read()

possible_amino_acid_number = floor(len(sequence)/3)
spare_bases = len(sequence) % 3

print("DNA sequence : " + sequence + "\n" +
      "Possible number of amino acids resulting from this sequence = " + str(possible_amino_acid_number) + "\n" +
      "Number of bases that do not make a whole amino acid = " + str(spare_bases)
      )
