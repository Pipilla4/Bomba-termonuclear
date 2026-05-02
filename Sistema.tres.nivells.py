import numpy as np
import matplotlib.pyplot as plt

N = 1000  # nombre de partícules
precisio  = 0.01  # precisió de la simulació

T = np.logspace(-1, 3, num=70)  # temperatures entre 0.01 i 1000

estats = np.array([0, 1, 10])  # estats possibles d'una partícula
X = np.random.choice(estats, size=N)  # estats inicials de les partícules

def P_acc(a,b,temperatura):
    delta = (a - b)
    if delta <= 0:
        return 1.0
    else:
        return np.exp(-delta / temperatura)

ocupacio0 = []
ocupacio1 = []
ocupacio2 = []
for t in T:
    x = X.copy()
    error = precisio + 1
    while error > precisio:
        E_inicial = np.sum(x)
        for i in range(2*N):
            particula = np.random.randint(N) # Triem una partícula a l'atzar
            canvi = np.random.choice([s for s in estats if s != x[particula]]) # Triem un nou estat per la partícula
            if np.random.rand() < P_acc(canvi,x[particula],t): # Acceptem o rebutgem el microstat
                x[particula] = canvi
        E_final = np.sum(x)
        error = abs(E_final - E_inicial) / E_inicial
    ocupacio0.append(np.sum(x == 0))
    ocupacio1.append(np.sum(x == 1))
    ocupacio2.append(np.sum(x == 10)) 

def f(T,E):
    return np.exp(-E/T)/(1 + np.exp(-1/T) + np.exp(-10/T))

plt.plot(T, ocupacio0,'.', label='Estat 0')
plt.plot(T, ocupacio1, '.', label='Estat 1')
plt.plot(T, ocupacio2, '.', label='Estat 2')

plt.plot(T, N*f(T,0),'--', label='Teoria 0')
plt.plot(T, N*f(T,1),'--', label='Teoria 1')
plt.plot(T, N*f(T,10),'--', label='Teoria 2')

plt.vlines(10*np.log(N),0,N,color='black',linestyle='dashed',label=r'$T_c$')
plt.xscale('log')
plt.xlabel(r'$kT/\varepsilon$') 
plt.ylabel('Ocupació')
plt.legend()
plt.show()