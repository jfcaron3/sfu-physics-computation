#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@title: PHYS 125 Computational Assignment
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
    prefix = "PHYS125_BallisticMotion"
    fullname = f"{prefix}_{name}_{timestamp}.{suffix}"
    fullpath = os.path.join(plotdir, fullname)
    if hasattr(fig, "savefig"):
        fig.savefig(fullpath)
    elif hasattr(fig, "save"):
        fig.save(fullpath)
    else:
        raise AttributeError(f"{fig} has no savefig or save method.")
    return


# Ballistic motion constants
g = -9.81  # m/s² gravity (note the sign)
v0 = 50  # Initial speed m/s
x0, y0 = 0, 0  # Initial coordinates
t0 = 0  # Initial time
npoints = 91  # Number of points to use in plots

# Calculate the trajectory for a given launch angle


def x(t, θ_d):
    """Return the x-coordinate of ballistic motion at
    given launch angle.

    The x-coordinate is just uniform motion with no
    acceleration x(t) = x_0 + v_0x*t"""
    θ_r = numpy.radians(θ_d)
    v0x = v0*numpy.cos(θ_r)
    return x0 + v0x*t


def y(t, θ_d):
    """Return the y-coordinate of ballistic motion at
    given launch angle.

    The xycoordinate is just motion with constant acceleration
    y(t) = y_0 + v_0y*t + ½at²"""
    θ_r = numpy.radians(θ_d)
    v0y = v0*numpy.sin(θ_r)
    return y0 + v0y*t + 0.5*g*t**2

# How much time should we plot?  Find t when y(t) = 0
# 0 = ½*g*t² + v0y*t + y0
# 0 = t² + (2*v0y/g)*t + 2*y0/g
# Use the quadratic formula with a = 1, b = 2*v0y/g, c = 2*y0
# t± = (-2*v0y/g ± √((2*v0y/g)² - 4*2*y0))/2
# Note that the result with - always gives 0, this is
# the other point at which y = 0, at launch.


def t_1(θ_d):
    """Return the projectile flight time as a function of
    launch angle."""
    θ_r = numpy.radians(θ_d)
    v0y = v0*numpy.sin(θ_r)
    tf = (-2*v0y/g + numpy.sqrt((2*v0y/g)**2 - 8*y0))/2
    return tf


θ_plot = 75  # Degrees
times = numpy.linspace(t0, t_1(θ_plot), npoints, endpoint=True)
xes = x(times, θ_plot)  # Array of x-coordinates
yes = y(times, θ_plot)  # Array of y-coordinates
fig1 = matplotlib.pyplot.figure()
axe1 = fig1.add_subplot()
fig1.suptitle(f"Coordinates at Fixed Projectile Angle {θ_plot}°")
axe1.plot(times, xes, label="x")  # Plot x vs t
axe1.plot(times, yes, label="y")  # Plot y vs t
axe1.set_xlabel("Time (s)")
axe1.set_ylabel("Distance (m)")
axe1.legend(loc="upper left")
savefig(fig1, "fig1", "png")

fig2 = matplotlib.pyplot.figure()
axe2 = fig2.add_subplot(aspect="equal")
fig2.suptitle(f"Path of Projectile at Angle {θ_plot}°")
# Note: retain the line2 object to make an animation later.
line2 = axe2.plot(xes, yes, label="xy")  # Plot y vs x
axe2.set_xlabel("x (m)")
axe2.set_ylabel("y (m)")
# Set the x and y axis limits equal to get a "square" plot
xmin, xmax = min(xes), max(xes)
ymin, ymax = min(yes), max(yes)
extra = 5  # A bit of extra space around the plot
themin = min(xmin, ymin)-extra
themax = max(xmax, ymax)+extra
axe2.set_xlim(themin, themax)
axe2.set_ylim(themin, themax)
savefig(fig2, "fig2", "png")

# Now we want to plot the maximum height and distance
# reached as a function of angle.
# Since we already calculated t_final, we can just
# plug that into x(t_final,θ) and y(t_final/2,θ).
θs_d = numpy.linspace(0, 90, npoints, endpoint=True)
x_maxes = x(t_1(θs_d), θs_d)  # Horizontal distance travelled
y_maxes = y(t_1(θs_d)/2, θs_d)  # Maximum heights reached
fig3 = matplotlib.pyplot.figure()
axe3 = fig3.add_subplot()
fig3.suptitle("Maximum x and y coordinates.")
axe3.plot(θs_d, x_maxes, label="x")
axe3.plot(θs_d, y_maxes, label="y")
axe3.set_xlabel("Angle (°)")
axe3.set_ylabel("Maximum distance (m)")
axe3.legend(loc="upper left")
savefig(fig3, "fig3", "png")


def update(frame):
    """Update the data in line2 based on animation frame number."""
    theline = line2[0]
    theline.set_data((xes[0:frame], yes[0:frame]))


fa = matplotlib.animation.FuncAnimation  # Shorten this long name
nframes = len(times)
dt = 1000*(times[-1]-times[0])/nframes  # Time between frames
ani = fa(fig=fig2, func=update, frames=nframes,
         interval=dt, repeat=False)
savefig(ani, "ani3", "gif")
