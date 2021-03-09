# brillouinzones

Includes python function to find the points marking the corners in any brillouin zone.

### Principle

The first brillouin zone is the area where the origin is the closest lattice point.
The n:th brillouin zone is the area where the origin is the n:th closest lattice point. 
To write my program one needs to realize two clever things:

##### Clever thing one

The corners of the zones are circle centers through (at least) 3 points.
One realizes this by the fact that at the corner of a zone,
the origin must be exactly as far away as two other points.
Thus, the corner must lie on the center of the circle intersecting these three points.
Further, if the circle contains k points in it's interior, 
the circle center is an outer corner of the k+1:th brillouin zone.
The reverse, however, does not hold and one needs to realize

##### Clever thing two

The number of lattice points on the edge of the circle is the number of brillouin zones the point is a (inner or outer) corner of.
The (outer) corners of the k:th brillouin zone thus contains the centers of circles with k-1 points inside as well as the 
centers of circles containing k-2 points, 
containing k-3 points and touching at least 4 lattice points, 
containing k-4 points and touching at least 5 lattice points, and so on.

### Pretty figures

Animation with first vector ```[1, 0]``` second vector ```[x, 1]```, where x goes from 0 to 1.
![sqr](https://user-images.githubusercontent.com/40766399/110310207-900d1b00-8002-11eb-8872-3e400ae2840e.gif)

60 brillouin zones of a hexagonal lattice.
![hex60](https://user-images.githubusercontent.com/40766399/110309704-f80f3180-8001-11eb-89e8-8165f33d1b34.png)

60 brillouin zones of an oblique lattice.
![oblique60](https://user-images.githubusercontent.com/40766399/110309759-065d4d80-8002-11eb-8d62-a651acc0c318.png)

100 brillouin zones of a square lattice.
![square100](https://user-images.githubusercontent.com/40766399/110309779-0b220180-8002-11eb-963f-6bf20f00cc85.png)

60 brillouin zones of a rectangular lattice.
![theEye_0-4_60](https://user-images.githubusercontent.com/40766399/110309836-1ffe9500-8002-11eb-9390-371047a4e8f0.png)
