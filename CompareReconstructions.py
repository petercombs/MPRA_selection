from Bio import SeqIO
from subprocess import Popen

comparisons = {
    "mammals": ("A1", "N1"),
    "primates": ("A2", "N4"),
    "haplorhini": ("A3", "N5"),  # (exclude Microcebus)
    "simiformi": ("A4", "N6"),  # (exclude Tarsius)
    "catarrhini": ("A5", "N7"),  # (Macaques, Hylobates, and Hominidae)
    "hominoidae": ("A6", "N19"),  # (Hylobates and Hominidae)
    "hominidae": ("A7", "N20"),  # (Hominidae)
    "homoninae": ("A8", "N21"),  # (Homininae)
    "hominini": ("A8", "N22"),  # (Hominini)
    "pan": ("A9", "N23"),  # (Pan)
    "newworld": ("A10", "N24"),  # (Callithrix and Saimiri)
}

smith = {
    rec.id: rec
    for rec in SeqIO.parse("Reference/SmithData/Data/ALDOB_Sequences.fasta", "fasta")
}
combs = {
    rec.id: rec
    for rec in SeqIO.parse("enhancers/ALDOB/merged/merged_seq.fasta", "fasta")
}

jobs = []
for comparison, (smithname, combsname) in comparisons.items():
    smith_seq = smith[smithname]
    combs_seq = combs[combsname]
    SeqIO.write([smith_seq, combs_seq], "tmp/{}.fasta".format(comparison), "fasta")
    jobs.append(
        Popen(
            [
                "clustalo",
                "--dealign",
                "--force",
                "-i",
                f"tmp/{comparison}.fasta",
                "-o",
                f"tmp/{comparison}.clustal",
                "--outfmt",
                "clustal",
            ]
        )
    )

for job in jobs:
    job.wait()
