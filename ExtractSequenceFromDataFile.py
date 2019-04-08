from sys import argv, stderr
from collections import defaultdict
from textwrap import wrap

if __name__ == "__main__":
    filetype = "?"
    chrom = "?"
    seq = defaultdict(lambda: set("ATGC"))
    for line in open(argv[2]):
        data = line.strip().split()
        try:
            pos = int(data[1])
        except:
            continue

        if filetype == "?":
            if len(data) == 10:
                filetype = "nadav"
            elif len(data) == 5:
                filetype = "patwardhan"

        if filetype == "nadav":
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
    seq = wrap("".join(seq[i].pop() for i in sorted(seq.keys())))
    print(f">Homo_sapiens chr{chrom}:{min}-{max} 100.0", *seq, sep="\n")
