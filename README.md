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

First vector ```[1, 0]``` second vector ```[x, 1]```, where x goes from 0 to 1.
![Lovely gif](sqr.gif)

60 brillouin zones of a hexagonal lattice.
![Hexagonal Lattice](hex60.pdf)

60 brillouin zones of an oblique lattice.
![Oblique Lattice](oblique60.pdf)

100 brillouin zones of a square lattice.
![Square Lattice](square100.pdf)
