import numpy as np
import matplotlib.pyplot as plt
from scipy . linalg import null_space

# Import the MarkovChain class from markovchain.py
from markovchain import MarkovChain
P = (1/15)*np.array([[8, 1, 6], [3, 5, 7], [4, 9, 2]]) # Transition matrix
mc = MarkovChain(P, ['0', '1', '2'],percentages=True)
mc.draw()

v = null_space(P - np.eye(3))
v = v/sum(v)
print("Limiting distribution is:", v)