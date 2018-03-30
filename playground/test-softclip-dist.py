from scipy.stats import poisson
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import math


def prob(max_softclip, observed_softclips, readlen):
    if max_softclip == 0:
        varcov = 0
    else:
        varcov = sum(c for l, c in observed_softclips.items() if l <= max_softclip) / max_softclip
    print(max_softclip, varcov)
    p = 1.0
    for s in range(readlen):
        softclip_count = observed_softclips[s]
        c = varcov if s <= max_softclip else 0
        p_ = poisson.pmf(softclip_count, c)
        if s <= max_softclip and softclip_count > c:
            p_ = 1.0
        p *= p_
    return p


observed_softclips = defaultdict(int)
#observed_softclips.update({2: 2000, 14: 1, 34: 1, 4: 1})
observed_softclips.update({1: 1073, 2: 1089, 3: 1042, 4: 1226, 5: 1360, 6: 531, 7: 536, 8: 539, 9: 601, 10: 650, 11: 323, 12: 369, 13: 389, 14: 395, 15: 409, 16: 241, 17: 245, 18: 251, 19: 246, 20: 234, 21: 188, 22: 184, 23: 181, 24: 166, 25: 175, 26: 129, 27: 140, 28: 121, 29: 130, 30: 126, 31: 95, 32: 100, 33: 84, 34: 95, 35: 91, 36: 59, 37: 63, 38: 52, 39: 76, 40: 51, 41: 59, 42: 51, 43: 50, 44: 40, 45: 46, 46: 24, 47: 30, 48: 19, 49: 28, 50: 39, 51: 21, 52: 20, 53: 9, 54: 10, 55: 7, 56: 12, 57: 2})
#observed_softclips.update({14: 1, 15: 4})
readlen = 75


plt.subplot(211)
x = np.arange(0, readlen + 1)
y = [observed_softclips.get(i, 0) for i in x]

plt.plot(x, y)

plt.subplot(212)
x = np.arange(0, readlen + 1)
y = np.array([prob(s, observed_softclips, readlen) for s in x])
marginal = math.fsum(y)
if marginal:
    y /= marginal

plt.plot(x, y)

plt.show()
