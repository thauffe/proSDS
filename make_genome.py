import sys
from Bio import SeqIO

fasta_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))

output = open(sys.argv[2], 'w')
id = sys.argv[1].split('/')
output.write('>' + id[-1] + '\n')

for key in fasta_dict:
    output.write(str(fasta_dict[key].seq))

output.close()