import Brillouin as B

import time

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


start = time.time()

n = 20
bz = B.BZ([0.4,0.8], n)
#bz = B.BZ([0.5,np.sqrt(3)/2], n)

end = time.time()
print(f"Time taken: {end - start}")

patches = []

for i in range(n, 0, -1):
    patches.append(Polygon(bz.getBrillouinZone(i)))

pc = PatchCollection(patches, alpha=1)
colors = np.linspace(0,1,len(patches))
pc.set_array(colors)

fig, ax = plt.subplots()
ax.add_collection(pc)
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
plt.show()
