from Bio import SeqIO
from sys import argv

if __name__ == "__main__":
    in_recs = SeqIO.parse(argv[1], 'clustal')
    SeqIO.write(in_recs, argv[2], 'fasta')

