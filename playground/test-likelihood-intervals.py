import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

likelihoods = pd.read_table("test-likelihood-intervals-somatic.txt", squeeze=True, names=["lh"])
afs = np.linspace(0.0, 1.0, len(likelihoods))

plt.plot(afs, likelihoods)
plt.show()
