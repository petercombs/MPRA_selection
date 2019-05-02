from sys import argv, stderr, stdout
from Bio import SeqIO, AlignIO
from Bio.Seq import MutableSeq, Seq
from json import dump

if __name__ == "__main__":
    aligned_seqs = [rec for rec in SeqIO.parse(argv[1], "fasta")]
    masked_seqs = {rec.id: rec for rec in SeqIO.parse(argv[2], "fasta")}

    mutseqs = {rec.id: MutableSeq(str(rec.seq)) for rec in aligned_seqs}
    SeqIO.write(masked_seqs.values(), stderr, 'fasta')

    for recid, rec in masked_seqs.items():
        aligncol = 0
        for i, base in enumerate(rec.seq):
            while aligncol < len(mutseqs[recid]) and  mutseqs[recid][aligncol] == "-":
                aligncol += 1
            if base == "N":
                print(f"Masking col {aligncol:03d} due to mask in {recid}", file=stderr)
                for aid, arec in mutseqs.items():
                    if arec[aligncol] != "-" and aid != "original":
                        arec[aligncol] = "N"
            aligncol += 1

    mutseqs2 = [
        SeqIO.SeqRecord(id=key, seq=value.toseq(), description="")
        for key, value in mutseqs.items()
    ]

    dump(
        {key: str(value) for key, value in mutseqs.items()},
        open(argv[3].replace("fasta", "json"), "w"),
    )
    SeqIO.write(mutseqs2, argv[3], "fasta")
