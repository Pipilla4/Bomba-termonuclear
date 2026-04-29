import numpy as np
import matplotlib.pyplot as plt

passos = 10000

T = 10 # kT/E

N = 1000 # Nombre de partícules
estats = [0,1,10] # Estats possibles d'una partícula
x = np.random.choice([0,1,2], size=N) # Microstat

def P_acc(y,x):
    return min(1,np.exp((-estats[y]+estats[x])/T))

for i in range(passos):
    particula = np.random.randint(N) # Triem una partícula a l'atzar
    canvi = np.random.choice([s for s in [0,1,2] if s != x[particula]]) # Triem un nou estat per la partícula
    if np.random.rand() < P_acc(canvi,x[particula]): # Acceptem o rebutgem el microstat
        x[particula] = canvi

plt.hist(x, bins=3)
plt.xlabel('Estat')
plt.ylabel('Partícules')
plt.show()