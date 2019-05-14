from Bio import SeqIO
from sys import argv

if __name__ == "__main__":
    default_score = 1000
    best_recs = {}
    for rec in SeqIO.parse(argv[1], 'fasta'):
        score = float(rec.description.split()[-1])
        rec.description = rec.id
        if best_recs.get(rec.id, (default_score, None))[0] > score:
            best_recs[rec.id] = (score, rec)
    SeqIO.write([r[1] for r in best_recs.values()], argv[2], 'fasta')
