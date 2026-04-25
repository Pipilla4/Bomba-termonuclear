from vpython import *
import numpy as np

# --- PARÁMETROS BÁSICOS (Q1: Natoms = 500) ---
win = 500
Natoms = 500  
L = 1  
gray = color.gray(0.7)
mass = 4E-3/6E23  
Ratom = 0.03  # Valor a discutir en Q1 [cite: 11]
k = 1.4E-23  
T = 300  # Temperatura de la colectividad canónica (Q2) 
dt = 1E-5
sigma = (k*T/mass)**(1/2)

animation = canvas(width=win, height=win, align='left')
animation.range = L
animation.title = 'Gas de esferas duras - Termostato de Andersen (Q2)'
animation.caption = "Distribución de velocidades (Izq) y Energía Total (Der)."

# --- CONFIGURACIÓN DEL CONTENEDOR ---
d = L/2+Ratom
r = 0.005
box_points = [vector(-d,-d,-d), vector(-d,-d,d), vector(d,-d,d), vector(d,-d,-d), vector(-d,-d,-d)]
boxbottom = curve(color=gray, radius=r, pos=box_points)
boxtop = curve(color=gray, radius=r, pos=[p + vector(0,2*d,0) for p in box_points])
for i in range(4):
    c = curve(color=gray, radius=r)
    c.append([box_points[i], box_points[i] + vector(0,2*d,0)])

Atoms = []
p = []
apos = []
pavg = sqrt(2*mass*1.5*k*T) 

for i in range(Natoms):
    x, y, z = L*random()-L/2, L*random()-L/2, L*random()-L/2
    if i == 0:
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=color.cyan, 
                            make_trail=True, retain=100, trail_radius=0.3*Ratom))
    else: 
        Atoms.append(sphere(pos=vector(x,y,z), radius=Ratom, color=gray))
    apos.append(vec(x,y,z))
    theta, phi = pi*random(), 2*pi*random()
    px = pavg*sin(theta)*cos(phi)
    py = pavg*sin(theta)*sin(phi)
    pz = pavg*cos(theta)
    p.append(vector(px,py,pz))

# --- CONFIGURACIÓN DE GRÁFICOS Y ESTADÍSTICA (Q2) ---
deltav = 100 
def barx(v): return int(v/deltav)
nhisto_bins = int(4500/deltav)
histo = [0.0]*nhisto_bins

# Gráfico de velocidades (Maxwell-Boltzmann)
gg = graph(width=win, height=0.4*win, xmax=3000, align='left',
           xtitle='v (m/s)', ytitle='N. Átomos')
theory = gcurve(color=color.blue, width=2)
for v in range(0, 3001, 10):
    theory.plot(v, (deltav/10)*Natoms*4*pi*((mass/(2*pi*k*T))**1.5) * exp(-0.5*mass*(v**2)/(k*T))*(v**2)*10)
vdist = gvbars(color=color.red, delta=deltav)

# Gráfico de Energía (Estudio de la Colectividad Canónica Q2)
ge = graph(width=win, height=0.4*win, align='right', 
           xtitle='Pasos (t)', ytitle='Energía Total (J)')
energy_plot = gcurve(color=color.green)
E_teorica_line = gcurve(color=color.black, dot=True)
E_teorica = (3/2) * Natoms * k * T 

accum = [[deltav*(i+.5), 0] for i in range(int(3000/deltav))]

def update_histogram():
    """Recalcula el conteo real de partículas por rango de velocidad."""
    global histo
    histo = [0.0] * nhisto_bins
    for i in range(Natoms):
        v = p[i].mag / mass
        bin_idx = barx(v)
        if bin_idx < len(histo):
            histo[bin_idx] += 1

def checkCollisions():
    hitlist = []
    r2 = (2*Ratom)**2
    for i in range(Natoms):
        ai = apos[i]
        for j in range(i):
            if mag2(ai - apos[j]) < r2: hitlist.append([i,j])
    return hitlist

# --- BUCLE PRINCIPAL ---
nhisto = 0
while True:
    rate(300)
    
    # Sincronización del histograma
    update_histogram()
    for i in range(len(accum)): 
        accum[i][1] = (nhisto * accum[i][1] + histo[i]) / (nhisto + 1)
    
    if nhisto % 10 == 0: 
        vdist.data = accum
        # Estudio de energía (Q2) 
        K_total = sum([mag2(pi)/(2*mass) for pi in p])
        energy_plot.plot(nhisto, K_total)
        E_teorica_line.plot(nhisto, E_teorica)
    
    nhisto += 1

    # Evolución temporal
    for i in range(Natoms): 
        apos[i] = apos[i] + (p[i]/mass)*dt
        Atoms[i].pos = apos[i]

    # Colisiones entre esferas
    hitlist = checkCollisions()
    for i, j in hitlist:
        vrel = (p[j]-p[i])/mass
        rrel = apos[i]-apos[j]
        if vrel.mag == 0 or rrel.mag > 2*Ratom: continue
        ptot = p[i] + p[j]
        rrel_unit = norm(rrel)
        p[i] = p[i] - dot(p[i] - ptot*0.5, rrel_unit) * rrel_unit * 2
        p[j] = ptot - p[i]

    # Paredes: Termostato de Andersen (Q2) 
    for i in range(Natoms):
        loc = apos[i]
        if abs(loc.x) > L/2:
            p[i].x = (1 if loc.x < 0 else -1) * abs(mass * np.random.normal(0, sigma))
            apos[i].x = (L/2 - 0.001) * np.sign(loc.x)
        if abs(loc.y) > L/2:
            p[i].y = (1 if loc.y < 0 else -1) * abs(mass * np.random.normal(0, sigma))
            apos[i].y = (L/2 - 0.001) * np.sign(loc.y)
        if abs(loc.z) > L/2:
            p[i].z = (1 if loc.z < 0 else -1) * abs(mass * np.random.normal(0, sigma))
            apos[i].z = (L/2 - 0.001) * np.sign(loc.z)
