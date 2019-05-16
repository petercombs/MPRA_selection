import pandas as pd
import numpy as np
from sys import argv
from re import compile
from matplotlib.pyplot import scatter, xlabel, ylabel, savefig, text, gca

numerator_direction = compile("; ?\(([0-9]*)/")
numerator_neutral = compile("\)/\(([0-9]*)/")

if __name__ == "__main__":
    kukns = {}
    kukn_ps = {}
    kukn_empirical = {}
    kdkns = {}
    kdkn_ps = {}
    kdkn_empirical = {}

    dus = {}
    dns = {}
    dds = {}

    pus = {}
    pns = {}
    pds = {}

    for fname in argv[1:]:
        enhancer = fname.split("/")[1]
        for line in open(fname):
            data = line.split()
            try:
                if line.startswith("Overall Ku/Kn"):
                    kukns[enhancer] = float(data[2])
                    kukn_ps[enhancer] = float(data[3][3:-1])
                    dus = int(numerator_direction.findall(line)[0])
                elif line.startswith("Overall Kd/Kn"):
                    kdkns[enhancer] = float(data[2])
                    kdkn_ps[enhancer] = float(data[3][3:].split(";")[0])
                    dds = int(numerator_direction.findall(line)[0])
                elif line.startswith("Possible up"):
                    pus[enhancer] = int(data[-1])
                elif line.startswith("Possible down"):
                    pds[enhancer] = int(data[-1])
                elif line.startswith("Possible neutral"):
                    pns[enhancer] = int(data[-1])
                elif line.startswith("Empirical KuKn"):
                    val = float(data[-1])
                    kukn_empirical[enhancer] = min(val, 1 - (val)) * 2
                elif line.startswith("Empirical KdKn"):
                    val = float(data[-1])
                    kdkn_empirical[enhancer] = min(val, 1 - (val)) * 2
            except ValueError:
                continue
    out = pd.DataFrame(
        {
            "KuKn": kukns,
            "KdKn": kdkns,
            "KuKn_fisher_p": kukn_ps,
            "KdKn_fisher_p": kdkn_ps,
            "KuKn_empirical_p": kukn_empirical,
            "KdKn_empirical_p": kdkn_empirical,
            "ObservedUp": dus,
            "PossibleUp": pus,
            "ObservedDown": dds,
            "PossibleDown": pds,
            "ObservedNeutral": dns,
            "PossibleNeutral": pns,
        }
    )

    out.to_csv("figures/data.tsv", sep="\t")
    scatter(
        -np.sign(np.log(out.KdKn)) * np.log10(out.KdKn_empirical_p),
        -np.sign(np.log(out.KuKn)) * np.log10(out.KuKn_empirical_p),
    )

    for ix in out.index:
        text(
            -np.sign(np.log(out.KdKn.loc[ix])) * np.log10(out.KdKn_empirical_p.loc[ix]),
            -np.sign(np.log(out.KuKn.loc[ix])) * np.log10(out.KuKn_empirical_p.loc[ix]),
            ix,
        )

    ax = gca()
    ax.spines["bottom"].set_position(("data", 0))
    ax.spines["left"].set_position(("data", 0))
    ax.spines["top"].set_alpha(0)
    ax.spines["right"].set_alpha(0)

    savefig("figures/data.png")
