from Bio import SeqIO, AlignIO
from Bio.SeqIO import SeqRecord
from Bio.Seq import Seq
from sys import argv, stderr
from re import compile
from collections import defaultdict
from os import path, makedirs
from subprocess import Popen
from collections import Counter

aligner = compile("FastML-[^/]*")
tcoffee = (
    "t_coffee -seq {dir}/{node}.fasta -outfile {dir}/{node}.out.aln -newtree {dir}/{node}.tree"
)


def get_consensus(aln_file):
    aln = AlignIO.read(aln_file, "clustal")
    out_bases = []
    default_seq = [rec for rec in aln if rec.id.endswith("tcoffee")][0]
    for base in range(len(aln[0])):
        bases = Counter(aln[:, base])
        if len(bases) == 1:
            out_bases.append(bases.most_common(1)[0][0])
            continue
        best, second = bases.most_common(2)
        if best[1] > second[1]:
            out_bases.append(best[0])
        else:
            print("Ambiguous base: {} {} {}".format(aln_file, base, bases), file=stderr)
            out_bases.append(default_seq[base])
    return "".join(out_bases)


if __name__ == "__main__":
    out_fname = argv[1]
    outdir = path.dirname(out_fname)
    in_fnames = argv[2:]

    all_seqs = defaultdict(list)

    for fname in in_fnames:
        prog = aligner.findall(fname)[0]
        for rec in SeqIO.parse(fname, "fasta"):
            old_id = rec.id
            rec.id = rec.id + "_" + prog
            all_seqs[old_id].append(rec)

    tempout = path.join(outdir, "tempseqs")
    makedirs(tempout, exist_ok=True)
    jobs = []
    for node, seqs in all_seqs.items():
        SeqIO.write(seqs, path.join(tempout, node + ".fasta"), "fasta")
        jobs.append(Popen(tcoffee.format(node=node, dir=tempout).split()))

    for job in jobs:
        job.wait()

    out_recs = []
    for node in all_seqs:
        seq = get_consensus(path.join(tempout, node + ".out.aln"))
        out_recs.append(SeqRecord(seq=Seq(seq), id=node))

    SeqIO.write(out_recs, out_fname, "fasta")
