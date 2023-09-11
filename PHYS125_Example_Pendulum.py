#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@title: PHYS 125 Example Code
@author: Jean-François Caron
"""

import numpy
Aθ = 10  # Swing amplitude in degrees
g = 9.81  # m/s2
l = 1  # Length of pendulum in metres
ω = numpy.sqrt(abs(g)/l)  # Angular frequency in radians/s

def θ(t):
    return Aθ*numpy.cos(ω*t)  # N.b.: cos takes radians
def x(t):
    return l*numpy.sin(numpy.radians(θ(t)))
def y(t):
    y = -l*numpy.cos(numpy.radians(θ(t)))
    #  Subtract l to set y = 0 at bottom of swing.
    return y+1


npoints = 100
T = 2*numpy.pi/ω  # Oscillation period in seconds
times = numpy.linspace(0, 2*T, npoints)
x_values = x(times)
y_values = y(times)


import os
import datetime
import matplotlib.pyplot
matplotlib.pyplot.style.use('tableau-colorblind10')

fig1 = matplotlib.pyplot.figure()
fig1.suptitle(f"{l}m pendulum with amplitude {Aθ}°")
axe1 = fig1.add_subplot()
axe1.plot(times, x_values, label="x")
axe1.plot(times, y_values, label="y")
axe1.set_xlabel("Time (s)")
axe1.set_ylabel("Deviation (m)")
axe1.legend(loc="lower right")

os.makedirs("plots", exist_ok=True)  # Make plots directory
timestamp = datetime.datetime.now().replace(microsecond=0).isoformat()
fig1.savefig(f"plots/PHYS125_Pendulum_Coordinates_{timestamp}.png")

fig2 = matplotlib.pyplot.figure()
fig2.suptitle(f"{l}m pendulum with amplitude {Aθ}°")
axe2 = fig2.add_subplot(aspect="equal")
# Draw 3 lines, one for path, and two for extremes of pendulum motion.
# Keep the plot return values for animating later.
y_values = y_values-1  # Restore the full pendulum length.
line2a = axe2.plot(x_values, y_values, label="xy")
line2b = axe2.plot([0, min(x_values)], [0, max(y_values)], label="left")
line2c = axe2.plot([0, max(x_values)], [0, max(y_values)], label="right")
axe2.set_xlabel("x (m)")
axe2.set_ylabel("y (m)")
axe2.set_xlim(-0.5, 0.5)
axe2.set_ylim(-1.0, 0)
fig2.savefig(f"plots/PHYS125_Pendulum_Path_{timestamp}.png")

import matplotlib.animation
def update(frame):
    la = line2a[0]
    lb = line2b[0]
    lc = line2c[0]
    x, y = x_values[frame], y_values[frame]
    la.set_data([x], [y])  # Just the pendulum bob
    lb.set_data([0, x], [0, y])  # The pendulum length
    # lc outlines the maximum movement.
    lc.set_data([min(x_values), 0, max(x_values)],
                [max(y_values), 0, max(y_values)])
    return

fa = matplotlib.animation.FuncAnimation
nframes = len(times)
dt = 1000*(times[-1]-times[0])/nframes
ani = fa(fig2, func=update, frames=nframes, interval=dt)
ani.save(f"plots/PHYS125_Pendulum_Animation_{timestamp}.gif")
