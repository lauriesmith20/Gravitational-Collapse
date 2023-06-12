import matplotlib.pyplot as plt
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D

# Parameters
n = 20 # No. of Particles
random.seed(100)
steps =20000
G = 6.67408 * 10 ** -11
dt = 25
soft = 10
vth = 1000
rho_g = 1e-6
rho_p = 5000
S = 10
size = np.zeros(n)
for s in range(n):
    size[s] = S
planetr = 1e8

combinationr = 2 * S

posplt = np.zeros((steps + 1, n, 3))
plot = np.empty((steps + 1, n, 3))

for x in range(n):
    posplt[0, x, :] = (
        (random.randint(0, 10000)/10), (random.randint(0, 10000)/10), (random.randint(0, 10000)/10))

radius = 1.52e11
circlex = np.linspace(0, 1e10, num=(n))
circley = []
negcircley = []
for z in circlex:
    y = ((radius ** 2) - (z ** 2)) ** (1 / 2)
    negy = -y
    circley.append(y)
    negcircley.append(negy)

#
#
# Masses
mass = np.zeros(n)
mass[:] = 4 / 3 * np.pi * (size ** 3) * rho_p

# Initial Velocities
velocity = np.zeros((steps + 1, n, 3))

# Energies for each Particle
GPE = np.zeros((steps, n))
KE = np.zeros((steps, n))
Energy = np.zeros((steps, n))

# Value Stores
momentum = np.zeros((steps+1, n, 3))
momentumSum = []
absmomentumSum =[]
permomentum =[]

radii = np.zeros((steps, n, n))  # radius between 2 particles
Fcollapse = np.zeros((steps + 1, n, n, 3))  # force between 2 particles
F = np.zeros((steps + 1, n, n, 3))
TotalGPE = []
TotalKE = []
TotalE = []
Eperc = []
names = ['1', '2', '3', '4', '5']
# Calculations
cup = np.zeros((n, n))
for x in range(steps):
    if x % 1000 == 0:
        print(x)

    for p in range(n):
        for z in range(p, n):
            dx = (posplt[x, z, 0] - posplt[x, p, 0])
            dy = (posplt[x, z, 1] - posplt[x, p, 1])
            dz = (posplt[x, z, 2] - posplt[x, p, 2])
            r = ((dx ** 2) + (dy ** 2) + (dz ** 2)) ** (1 / 2)
            radii[x, p, z] = r

            if r != 0:
                if r < combinationr:
                    if cup[p, z] == 0:
                        velocity[x, p, :] = ((mass[p] * velocity[x, p, :]) + (mass[z] * velocity[x, z, :])) / (
                                mass[p] + mass[z])
                        velocity[x, z, :] = (mass[p] * velocity[x, p, :] + mass[z] * velocity[x, z, :]) / (
                                mass[p] + mass[z])

                        # velocity[x + 1, z, :] = velocity[x + 1, p, :]
                        mass[p] += mass[z]
                        mass[z] = 0
                        posplt[x, z, :] = posplt[x, p, :]
                        cup[p, z] = 1
                        newsize = (size[p]**3 + size[z]**3)**(1/3)
                        size[p] = newsize
                        size[z] = 0
                        print(x, mass[p], mass[z])
                    else:
                        posplt[x, z, :] = posplt[x, p, :]

                if cup[p, z] == 1:
                    posplt[x, z, :] = posplt[x, p, :]

                Fcollapse[x, p, z, :] = (G * mass[p] * mass[z] / abs((r ** 2 + soft ** 2) ** (3 / 2))) * np.array(
                    (dx, dy, dz))
                Fcollapse[x, z, p, :] = -Fcollapse[x, p, z, :]
                GPE[x, p] += -G * mass[p] * (mass[z] / r)
                # GPE[x,z] += -G*mass[p]*((mass[z]/r))
        Fdx = -(4 / 3) * np.pi * (size[p] ** 2) * vth * rho_g * velocity[x, p, 0]
        Fdy = -(4 / 3) * np.pi * (size[p] ** 2) * vth * rho_g * velocity[x, p, 1]
        Fdz = -(4 / 3) * np.pi * (size[p] ** 2) * vth * rho_g * velocity[x, p, 2]

        Ftot = np.array(
            (sum(Fcollapse[x, p, :, 0]) + Fdx, sum(Fcollapse[x, p, :, 1]) + Fdy, sum(Fcollapse[x, p, :, 2]) + Fdz))

        # print(mass[p])
        # print(Fdx, Fdy, Fdz)
        # print(Ftot)
        if mass[p] != 0:
            dv = Ftot * (dt / mass[p])
        if mass[p] == 0:
            dv = 0

        if sum(velocity[x + 1, p, :]) == 0:
            newv = (velocity[x, p, :] + dv)
            velocity[x + 1, p, :] = newv

        newpos = posplt[x, p, :] + (velocity[x + 1, p, :] * dt)
        posplt[x + 1, p, :] = newpos

        # totalV = (velocity[x, p, 0] ** 2 + velocity[x, p, 1] ** 2 + velocity[x, p, 2] ** 2) ** (1 / 2)

        # Energy
        momentum[x, p, :] = (mass[p]*velocity[x,p,0],mass[p]*velocity[x,p,1],mass[p]*velocity[x,p,2])

        magv = (velocity[x, p, 0] ** 2 + velocity[x, p, 1] ** 2 + velocity[x, p, 2] ** 2) ** (1 / 2)
        KE[x, p] = 0.5 * mass[p] * (magv ** 2)
        Energy[x, p] = KE[x, p] + GPE[x, p]

    ptot = ((sum(momentum[x,:,0]))**2 + (sum(momentum[x,:,1]))**2 + (sum(momentum[x,:,2]))**2)**(1/2)
    absptot = ((sum(abs(momentum[x,:,0])))**2 + (sum(abs(momentum[x,:,1])))**2 + (sum(abs(momentum[x,:,2])))**2)**(1/2)
    if absptot != 0:
        perp = (ptot / absptot) *100
    else:
        perp = 0
    if abs(perp) ==100:
        perp = 0

    momentumSum.append(ptot)
    absmomentumSum.append(absptot)
    permomentum.append(perp)
    # Energy Totals
    TotalKE.append(sum(KE[x, :]))
    TotalGPE.append(sum(GPE[x, :]))
    TotalE.append(sum(Energy[x, :]))
    # Eperc.append(100 * (TotalE[x] - TotalE[0]) / TotalE[0])

# GRAPHS
print(mass[:])
print(size[:])

colors = ['blue', 'red', 'orange', 'green']

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for x in range(n):
    ax.scatter(posplt[0, x, 0], posplt[0, x, 1], posplt[0, x, 2])
    ax.plot(posplt[:, x, 0], posplt[:, x, 1], posplt[:, x, 2])

plt.figure()
plt.axis('equal')
for x in range(n):
    # plt.plot(circlex, circley, color='black')
    # plt.plot(circlex, negcircley, color='black')
    plt.plot(posplt[:, x, 0], posplt[:, x, 1])
    plt.scatter(posplt[0, x, 0], posplt[0, x, 1])

# plt.figure()
# plt.axis('equal')
# for n in range(n):
# plt.plot(posplt[-9700:steps-1, n, 0], posplt[-9700:steps-1, n, 1])

#t = np.linspace(0, (steps * dt), num=steps)  # Time Values
plt.show()

plt.figure()
t = np.arange(1,steps+1)
plt.plot(t, permomentum)
plt.xlabel("Number of Steps")
plt.ylabel("Percentage Error in Momentum (%)")
plt.show()

# plt.figure()
# plt.plot(t, TotalGPE)
# plt.show()

# plt.figure()
# plt.plot(t, TotalE)
# plt.show()
# print(Eperc[-1:]) #Change in Energy
#
