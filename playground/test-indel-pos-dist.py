from scipy.stats import poisson
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import math
from functools import partial

def prob(start, end, observed_indels, readlen):
    if start == end:
        varcov = 0
    elif start > end:
        return 0.0
    else:
        varcov = sum(c for p, c in observed_indels.items() if p >= start and p < end) / (end - start)
    p = 1.0
    for s in range(readlen):
        count = observed_indels[s]
        c = varcov if s >= start and s < end else 0
        p *= poisson.pmf(count, c)
    return p


observed = defaultdict(int)
observed.update({25: 1, 52: 1, 46: 1, 42: 1, 37: 1})
readlen = 100

y = np.array([prob(start, end, observed, readlen) for start in range(readlen) for end in range(readlen)])
y = y.reshape((readlen, readlen))
marginal = math.fsum(y.reshape(10000))
y /= marginal

plt.imshow(y)
plt.show()
