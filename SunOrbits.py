import matplotlib.pyplot as plt
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D

# Parameters
n = 2  # No. of Particles
pn = 1
steps = 63000
G = 6.67408 * 10 ** -11
dt = 5e4
e = 1e4
vth = 100
rho_g = 1e-6
rho_p = 5000
size = 0.1
planetr = 1e8
drag_v = 30

posplt = np.zeros((steps + 1, n, 3))
posplt[0, 0, :] = (0, 0, 0)
for x in range(1, pn):
    posplt[0, x, :] = (0, 1.5e11 + (1e9 * (x-1)), 0)

radius = 1.52e11
circlex = np.linspace(0, 1e10, num= (n - pn))
circley = []
negcircley = []
for z in circlex:
    y = ((radius ** 2) - (z ** 2)) ** (1 / 2)
    negy = -y
    circley.append(y)
    negcircley.append(negy)

for s in range(pn, n):
    posplt[0, s, :] = (circlex[s - pn], circley[s - pn] , 0)
#
#
# Masses
mass = np.zeros(n)
mass[0] = 2e30
mass[1:pn] = 1e22

mass[pn:n] = 4 / 3 * np.pi * (size ** 3) * rho_p


def Orbitals(orbitalposition):
    orbitalradius = (orbitalposition[0] ** 2 + orbitalposition[1] ** 2 + orbitalposition[2] ** 2) ** (1 / 2)
    if orbitalposition[0] != 0 and orbitalposition[1] != 0:

        if orbitalposition[0] > 0 and orbitalposition[1] > 0:
            signx = 1
            signy = -1
        elif orbitalposition[0] > 0 > orbitalposition[1]:
            signx = -1
            signy = 1
        elif orbitalposition[0] < 0 and orbitalposition[1] < 0:
            signx = -1
            signy = 1
        elif orbitalposition[0] < 0 < orbitalposition[1]:
            signx = 1
            signy = -1

    elif orbitalposition[0] == 0 and orbitalposition[1] != 0:
        signx = (orbitalposition[1] / abs(orbitalposition[1]))
        signy = 1

    elif orbitalposition[1] == 0 and orbitalposition[0] != 0:
        signx = 1
        signy = -(orbitalposition[0] / abs(orbitalposition[0]))

    if orbitalposition[1] != 0:
        answer = ((signx * ((G * mass[0] / orbitalradius) ** (1 / 2)) * np.cos(
            ((((G * mass[0]) / orbitalradius ** 3) ** (1 / 2)) * dt)
            + np.arctan(orbitalposition[0] / orbitalposition[1]))),

                  (signy * ((G * mass[0] / orbitalradius) ** (1 / 2)) * np.sin(
                      ((((G * mass[0]) / orbitalradius ** 3) ** (1 / 2)) * dt)
                      + np.arctan(orbitalposition[0] / orbitalposition[1]))),

                  0)
    else:
        answer = ((signx * ((G * mass[0] / orbitalradius) ** (1 / 2)) * np.cos(
            ((((G * mass[0]) / orbitalradius ** 3) ** (1 / 2)) * dt)
            + np.pi / 2)),

                  (signy * ((G * mass[0] / orbitalradius) ** (1 / 2)) * np.sin(
                      ((((G * mass[0]) / orbitalradius ** 3) ** (1 / 2)) * dt)
                      + np.pi / 2)),

                  0)
    return answer


# Initial Velocities
velocity = np.zeros((steps + 1, n, 3))
velocity[0, 0, :] = (0, 0, 0)

for v in range(1, pn):
    velocity[0, v, :] = Orbitals((posplt[0, v, :]))
    velocity[0, v, :] = Orbitals((posplt[0, v, :]))

for s in range(pn, n):
    velocity[0, s, :] = Orbitals(posplt[0, s, :])

# Energies for each Particle
GPE = np.zeros((steps, n))
KE = np.zeros((steps, n))
Energy = np.zeros((steps, n))

# Value Stores
radii = np.zeros((steps, n, n))  # radius between 2 particles
Fcollapse = np.zeros((steps + 1, pn, pn, 3))  # force between 2 particles
F = np.zeros((steps + 1, (n - pn), pn, 3))
TotalGPE = []
TotalKE = []
TotalE = []
Eperc = []
names = ['1', '2', '3', '4', '5']
# Calculations
for x in range(steps):
    if x % 1000 == 0:
        print(x)
    for p in range(pn):
        for z in range(p, pn):
            dx = (posplt[x, z, 0] - posplt[x, p, 0])
            dy = (posplt[x, z, 1] - posplt[x, p, 1])
            dz = (posplt[x, z, 2] - posplt[x, p, 2])
            r = ((dx ** 2) + (dy ** 2) + (dz ** 2)) ** (1 / 2)
            radii[x, p, z] = r
            if r != 0:
                Fcollapse[x, p, z, :] = (G * mass[p] * mass[z] / abs((r ** 2 + e ** 2) ** (3 / 2))) * np.array(
                    (dx, dy, dz))
                Fcollapse[x, z, p, :] = -Fcollapse[x, p, z, :]
                GPE[x, p] += -G * mass[p] * (mass[z] / r)
                # GPE[x,z] += -G*mass[p]*((mass[z]/r))
        Ftot = np.array((sum(Fcollapse[x, p, :, 0]), sum(Fcollapse[x, p, :, 1]), sum(Fcollapse[x, p, :, 2])))

        # print(Fdx, Fdy, Fdz)
        # print(Ftot)

        dv = Ftot * (dt / mass[p])

        newv = (velocity[x, p, :] + dv)
        velocity[x + 1, p, :] = newv

        newpos = posplt[x, p, :] + (newv * dt)
        posplt[x + 1, p, :] = newpos

    for p in range(pn, n):
        for z in range(pn):
            dx = (posplt[x, z, 0] - posplt[x, p, 0])
            dy = (posplt[x, z, 1] - posplt[x, p, 1])
            dz = (posplt[x, z, 2] - posplt[x, p, 2])
            r = ((dx ** 2) + (dy ** 2) + (dz ** 2)) ** (1 / 2)
            radii[x, z, p] = r
            if r != 0:
                F[x, p - pn, z, :] = (G * mass[p] * mass[z] / abs((r ** 2 + e ** 2) ** (3 / 2))) * np.array(
                    (dx, dy, dz))

            if r < planetr:
                mass[z] += mass[p]
                mass[p] = 0
                posplt[x, p, :] = (0, 0, 0)
                velocity[x, p, :] = (0, 0, 0)

        totalV = (velocity[x, p, 0] ** 2 + velocity[x, p, 1] ** 2 + velocity[x, p, 2] ** 2) ** (1 / 2)
        Fdx = -(4 / 3) * np.pi * (size ** 2) * vth * rho_g * (velocity[x, p, 0] / totalV) * drag_v
        Fdy = -(4 / 3) * np.pi * (size ** 2) * vth * rho_g * (velocity[x, p, 1] / totalV) * drag_v
        Fdz = -(4 / 3) * np.pi * (size ** 2) * vth * rho_g * (velocity[x, p, 2] / totalV) * drag_v

        Ftot = np.array((sum(F[x, p - pn, :, 0]) + Fdx, sum(F[x, p - pn, :, 1]) + Fdy, sum(F[x, p - pn, :, 2]) + Fdz))

        # print(Fdx, Fdy, Fdz)
        # print(Ftot)
        if mass[p] != 0:
            dv = Ftot * (dt / mass[p])
        if mass[p] == 0:
            dv = 0

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
    # Eperc.append(100 * (TotalE[x] - TotalE[0]) / TotalE[0])

# GRAPHS

colors = ['blue','red','orange','green']
for i in range(pn,n):
    colors.append('black')

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
    plt.plot(posplt[:, x, 0], posplt[:, x, 1], color= colors[x])
    plt.scatter(posplt[0, x, 0], posplt[0, x, 1])

t = np.linspace(0, (steps * dt), num=steps)  # Time Values
plt.show()

# plt.figure(figsize=(13,10))
##plt.figure()
# plt.plot(t, TotalKE)
# plt.show()

# plt.figure()
# plt.plot(t, TotalGPE)
# plt.show()

# plt.figure()
# plt.plot(t, TotalE)
# plt.show()
# print(Eperc[-1:]) #Change in Energy
#
