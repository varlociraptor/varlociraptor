import matplotlib.pyplot as plt
import numpy as np

def prob(obs, af, prob_sample_alt):
    p = 1
    for o in obs:
        p *= af * prob_sample_alt * o[1] + (1 - af * prob_sample_alt) * o[0]
    return p


ref_obs = [(1, 0) for _ in range(50)]
alt_obs = [(0, 1) for _ in range(50)]
obs = ref_obs + alt_obs
obs = ref_obs

def plot_prob(prob_sample_alt):
    plt.figure()
    x = np.linspace(0.0, 1.0, 10)
    y = [prob(obs, af, prob_sample_alt) for af in x]
    plt.plot(x, y)
    plt.show()


plot_prob(0.0001)
