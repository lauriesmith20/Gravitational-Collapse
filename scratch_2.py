import matplotlib.pyplot as plt
import numpy as np
import random
from matplotlib.animation import FuncAnimation
from IPython import display







# Parameters
n = 10 # No. of Particles
steps = 100
G = 6.67408 * 10 ** -11
dt = 10
e = 10


posplt = np.zeros((steps+1,n,3))
for s in range(n):
    posplt[0,s,:] = (random.randint(1,10000),random.randint(1,10000),0)

#Initial Velocities
velocity = np.zeros((steps+1,n,3))
for s in range(n):
    velocity[1,s,:]= (0,0,0)


# (Inital)Positions
#posplt = np.zeros((steps + 1, n, 3))
#for s in range(n):
   # posplt[0, s, :] = (random.randint(1, 100000), random.randint(1, 100000), 0)

# Initial Velocities
#velocity = np.zeros((steps + 1, n, 3))
#for s in range(n):
    #velocity[1, s, :] = (0, 0, 0)

# Masses
mass = np.zeros((n))
mass[:]=1e10

# Energies for each Particle
GPE = np.zeros((steps, n))
KE = np.zeros((steps, n))
Energy = np.zeros((steps, n))

# Value Stores
radii = np.zeros((steps, n, n))  # radius between 2 particles
F = np.zeros((steps + 1, n, n, 3))  # force between 2 particles
TotalGPE = []
TotalKE = []
TotalE = []
Eperc = []
names = ['1', '2', '3', '4', '5']
# Calculations
for x in range(steps):
    if x%1000 ==0:
        print(x)
    for p in range(n):
        for z in range(p, n):
            dx = (posplt[x, z, 0] - posplt[x, p, 0])
            dy = (posplt[x, z, 1] - posplt[x, p, 1])
            dz = (posplt[x, z, 2] - posplt[x, p, 2])
            r = ((dx ** 2) + (dy ** 2) + (dz ** 2)) ** (1 / 2)
            radii[x, p, z] = r
            if r != 0:
                F[x, p, z, :] = (G * mass[p] * mass[z] / abs((r ** 2 + e ** 2) ** (3 / 2))) * np.array((dx, dy, dz))
                F[x, z, p, :] = -F[x, p, z, :]
                GPE[x, p] += -G * mass[p] * ((mass[z] / r))
                # GPE[x,z] += -G*mass[p]*((mass[z]/r))

        Ftot = np.array((sum(F[x, p, :, 0]), sum(F[x, p, :, 1]), sum(F[x, p, :, 2])))
        # print(Ftot)

        dv = Ftot * (dt / mass[p])

        newv = (velocity[x, p, :] + dv)
        velocity[x + 1, p, :] = newv

        newpos = posplt[x, p, :] + (newv * dt)
        posplt[x + 1, p, :] = newpos

        # Energy
        magv = (velocity[x, p, 0] ** 2 + velocity[x, p, 1] ** 2 + velocity[x, p, 2] ** 2) ** (1 / 2)
        KE[x, p] = 0.5 * mass[p] * (magv ** 2)
        Energy[x, p] = KE[x, p] + GPE[x, p]

    # Energy Totals
    TotalKE.append(sum(KE[x, :]))
    TotalGPE.append(sum(GPE[x, :]))
    TotalE.append(sum(Energy[x, :]))
    Eperc.append(100 * (TotalE[x] - TotalE[0]) / TotalE[0])


# GRAPHS
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg'


fig = plt.figure()
lines=plt.plot([], "o")
line = lines[0]

def animate(frame):
    x = posplt[frame, 0, 0]
    y = posplt[frame, 0, 1]
    line.set_data((x,y))


anim = FuncAnimation(fig, animate, frames=100, interval=20)
video = anim.to_html5_video()
html =  display.HTML(video)
display.display(html)
