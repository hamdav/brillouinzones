import Brillouin as B

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

bz = B.BZ([0.4,0.8])

p = Polygon(bz.bzones[1])
pc = PatchCollection([p], alpha=1)

fig, ax = plt.subplots()
ax.add_collection(pc)
plt.show()

