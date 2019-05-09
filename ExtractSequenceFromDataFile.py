from sys import argv, stderr
from collections import defaultdict
from textwrap import wrap

if __name__ == "__main__":
    filetype = "?"
    chrom = "?"
    for line in open(argv[2]):
        data = line.strip().split()
        try:
            pos = int(data[1])
        except:
            continue

        if filetype == "?":
            if len(data) == 10:
                filetype = "kircher"
                seq = defaultdict(lambda: set("N"))
            elif len(data) == 5:
                filetype = "patwardhan"
                seq = defaultdict(lambda: set("ATGC"))

        if filetype == "kircher":
            if data[9] == argv[1]:
                chrom = data[0]
                seq[pos] = {data[2]}
        elif filetype == "patwardhan":
            if data[0].upper() == argv[1]:
                seq[pos].discard(data[2])
            #else:
                #print("Skipping row for ", data[0].upper())
        else:
            raise ValueError("We cannot infer the correct filetype")

    min = min(seq.keys())
    max = max(seq.keys())
    # We know that TCF7L2 has a missing base in the data file, so I'm going to
    # use range instead of just looping over the bases.
    seq = wrap("".join(seq[i].pop() for i in range(min, max+1)))
    print(f">{argv[3]} chr{chrom}:{min}-{max} 100.0", *seq, sep="\n")
