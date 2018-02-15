from scipy.stats import poisson
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import math


def prob(max_softclip, observed_softclips, coverage, allelefreq, readlen):
    varcov = coverage * allelefreq
    p = 1.0
    for s in range(readlen):
        softclip_count = observed_softclips[s]
        c = varcov if s <= max_softclip else 0
        p *= poisson.pmf(softclip_count, c)
    return p


observed_softclips = defaultdict(int)
#observed_softclips.update({14: 1, 34: 1, 4: 1})
coverage = 0.3
allelefreq = 0.5
readlen = 100

x = np.arange(0, 100)
y = np.array([prob(s, observed_softclips, coverage, allelefreq, readlen) for s in x])
marginal = math.fsum(y)
y /= marginal

plt.plot(x, y)
plt.show()
