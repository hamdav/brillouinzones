import time
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.ticker import MaxNLocator

import Brillouin as B


############################
#   User variables
############################

# Set the number of brillouin zones to plot
# A value of 80 takes about 70 seconds on my machine,
# 20 take about 4 seconds
n = 20

# Set the second primitive basis vector.
# The first is always [1,0].
# If you want some other first primitive basis vector,
# the author suggests you tilt your head and zoom in

for x in np.linspace(-0.5,0.5,21):
    v = [x, math.sqrt(3)/2]


    ############################
    #   Calculation of points
    ############################

    start = time.time()

    bz = B.getBrillouinZonePoints(v, n)

    end = time.time()
    print(f"Time taken: {end - start}")

    ############################
    #   Plotting
    ############################

    plt.style.use('seaborn-darkgrid')

    patches = []

    # Plot "backwards" so the later zones do not cover the earlier
    for i in range(n, 0, -1):
        patches.append(Polygon(bz[i]))

    pc = PatchCollection(patches, alpha=1)

    # Construct the colors
    cmap = plt.get_cmap("tab20b")
    colors = cmap(np.tile(np.linspace(0,1,20), (n-1) // 20 + 1))
    pc.set_color(colors[::-1])

    # Create the figure and axis
    fig, ax = plt.subplots()

    # Add the collection of polygons
    ax.add_collection(pc)

    # Set limits
    xmax = max([p[0] for p in bz[n]])
    ymax = max([p[1] for p in bz[n]])
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-ymax, ymax)

    # Make aspect ratio 1:1
    ax.axis('equal')

    # Make axis fill figure
    ax.set_position([0, 0, 1, 1])

    # Set ticks to the integers
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))

    # Show the plot
    #plt.show()

    # Uncomment to save the figure to disk
    plt.savefig(f"brillouin_{v[0]:.2f}-{v[1]:.2f}.pdf")
