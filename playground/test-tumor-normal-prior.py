import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math


HETEROZYGOSITY = 0.001
PLOIDY = 2
EFFMUT = 350000
GENOME_N = 3.5e9


def prior_normal(f):
    m = f * PLOIDY
    if f == 0.0:
        return HETEROZYGOSITY * sum(1.0 / m for m in range(1, PLOIDY + 1))
    else:
        return HETEROZYGOSITY / m


def prior(f_tumor, f_normal):
    f_somatic = abs(f_tumor - f_normal)
    p = EFFMUT * (1 / (f_somatic ** 2)) / GENOME_N
    return p * prior_normal(f_normal)
    #p = math.log(EFFMUT) - (2 * math.log(f_somatic) + math.log(GENOME_N))
    #return math.exp(p) * prior_normal(f_normal)


for f_normal in [0.0, 0.5, 1.0]:
    x = np.linspace(0.0, 1.0, 100)
    y = [prior(f, f_normal) for f in x]
    print(math.fsum(y))
    plt.semilogy(x, y, "-", label="normal={}".format(f_normal))

sns.despine()
plt.legend()
plt.xlabel("tumor allele freq")
plt.ylabel("density")

plt.savefig("tumor-normal-prior.svg", bbox_inches="tight")
