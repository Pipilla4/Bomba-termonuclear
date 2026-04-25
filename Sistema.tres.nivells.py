import numpy as np
import matplotlib.pyplot as plt

passos = 1000

T = 300
beta = 1/(1.38e-23*T)

N = 1000 # Nombre de partícules
estats = np.array([0,1,10]) # Estats possibles d'una partícula
x = np.array([int(N/3)+1,int(N/3),int(N/3)]) # Microstat (Nombre de partícules en cada estat)
def H(x):
    return np.dot(estats,x) #Energia del microstat
def P(x):
    return np.exp(-beta*H(x)) #Probabilitat del microstat
def P_acc(y,x):
    return min(1,P(y)/P(x))

for i in range(passos):
    y = [0,0,0]
    for j in range(N):
        y[np.random.randint(3)] += 1 # Generem un nou microstat aleatori
    if np.random.rand() < P_acc(y,x): # Acceptem o rebutgem el microstat
        x=y

print(x)