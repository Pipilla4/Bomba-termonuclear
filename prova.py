print('hello')
import numpy as np
import matplotlib.pyplot as plt

x = np.array([1,2,3,3,1,1,3,1,2,3,1,1,2,3,3,1,2,2,3,1])
valors, comptes = np.unique(x, return_counts=True)
print(comptes)