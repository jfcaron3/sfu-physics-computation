#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@title: PHYS 126 Electrodynamics
@author: Jean-François Caron
"""

import os
import datetime
import matplotlib.pyplot
import matplotlib.animation
import numpy

# Set colourblind-friendly colours.
matplotlib.pyplot.style.use('tableau-colorblind10')

# Set up timestamp for files and plot directory.
timestamp = datetime.datetime.now().replace(microsecond=0).isoformat()
plotdir = "plots"
os.makedirs(plotdir, exist_ok=True)


def savefig(fig, name, suffix):
    """Given a figure or other saveable object, a name, and
    and a suffix for file format, puts together the full path
    and saves the object."""
    prefix = "PHYS126_Electrodynamics"
    fullname = f"{prefix}_{name}_{timestamp}.{suffix}"
    fullpath = os.path.join(plotdir, fullname)
    if hasattr(fig, "savefig"):
        fig.savefig(fullpath)
    elif hasattr(fig, "save"):
        fig.save(fullpath)
    else:
        raise AttributeError(f"{fig} has no savefig or save method.")
    return


c = 299792458  # m/s speed of light
mu0 = 1.25663706212e-6  # NA^2 magnetic constant
e0 = 1/(mu0*c*c)  # electric constant
k = 1/(4*numpy.pi*e0)  # Coulomb constant


class Particle:
    def __init__(self, mass, charge, position, velocity):
        self.mass = mass
        self.charge = charge
        self.position = numpy.array(position)
        self.velocity = numpy.array(velocity)

    def move(self, E_field, B_field):
        self.position += self.velocity*dt
        E_force = self.charge*E_field
        B_force = self.charge*numpy.cross(self.velocity, B_field)
        total_force = E_force + B_force
        acceleration = total_force/self.mass
        self.velocity += acceleration*dt

    def E_field(self, p):
        r = p - self.position  # Vector from particle to p
        r2 = numpy.dot(r, r)
        rhat = r/numpy.sqrt(r2)
        return k*self.charge/r2 * rhat


# This simulation can handle multiple charges and the
# electric fields they produce, external electric fields,
# and external magnetic fields.  Notably it does NOT include
# the magnetic fields produced by the moving charges.

p1_start = [1.0, 0.0]
p1_v0 = [0.0, 0.1]
p2_start = [-1.0, 0.0]
p2_v0 = [0.0, -0.0]
p1 = Particle(1, 1e-6, p1_start, p1_v0)
p2 = Particle(1, -1e-6, p2_start, p2_v0)

# Calculate potential and kinetic energies of system
U1 = p1.charge*p2.E_field(p1.position)
U2 = p2.charge*p1.E_field(p2.position)
Utot = U1+U2
K1 = 0.5*p1.mass*numpy.dot(p1.velocity, p1.velocity)
K2 = 0.5*p2.mass*numpy.dot(p2.velocity, p2.velocity)
Ktot = K1+K2
print(f"Total potential energy is {Utot} J")
print(f"Total kinetic energy is {Ktot} J")

dt = 1  # Time step to use for simulation, in seconds
t0, t1 = 0, 1000  # Total simulation time.
t = t0  # Elapsed simulation time
n_iters = int((t1-t0)/dt)
i = 0  # Iteration number
p1_coord = numpy.zeros((n_iters, len(p1_start)))
p2_coord = numpy.zeros((n_iters, len(p1_start)))

while t < t1:
    p1_coord[i] = p1.position
    p2_coord[i] = p2.position

    # r12 = p1.position - p2.position
    # d12 = numpy.sqrt(numpy.dot(r12, r12))
    # print(i, t, d12)
    E12 = p1.E_field(p2.position)  # Field at p2 from p1
    E21 = p2.E_field(p1.position)  # Field at p1 from p2
    p1.move(E21, [0, 0])
    p2.move(E12, [0, 0])
    t += dt
    i += 1
p1_coord = numpy.transpose(p1_coord)
p2_coord = numpy.transpose(p2_coord)

fig1 = matplotlib.pyplot.figure()
axe1 = fig1.add_subplot()
fig1.suptitle("Motion of two charged particles")
axe1.plot(p1_coord[0], p1_coord[1], label="p1")
axe1.plot(p2_coord[0], p2_coord[1], label="p2")
axe1.set_xlabel("x (m)")
axe1.set_ylabel("y (m)")
axe1.legend(loc="upper left")
savefig(fig1, "fig1", "png")
