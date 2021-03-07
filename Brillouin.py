import numpy as np
import math

def getBrillouinZonePoints(v=[0, 1], n=3):
    """
    getBrillouinZonePoints returns a dictionary with integers from
    1 to n (inclusive) as keys and lists of points, sorted by angle
    to positive x-axis, as values. The list of points with key i are 
    the outer boundry of the i:th Brillouin zone.
    The inner boundry of the i:th Brillouin zone would, of course, be
    the outer boundry of the i-1:th Brillouin zone.

    The argument v specifies the second primitive vector of the lattice.
    The first is always taken to be [1, 0].
    """

    # Set the tolerance
    tol = 1e-10

    # Set the primitive basis vectors. 
    # The vector v is transformed to have positive y and x between -0.5 and 0.5.
    u = np.array([1, 0])
    v = np.array(v)
    if v[1] < 0:
        v = -v
    v[0] = v[0] - math.floor(v[0] + 0.5)


    # bzones will be the return value
    bzones = dict()

    # If we are only interested in circles containing up to n points,
    # there is no need to check lattice points beyond some k_u u coordinate
    # and k_v v coordinate.
    # A circle touching u coordinate k_u has diameter at least k_u
    # and thus a s = k_u / sqrt(2) sided square fits inside it.
    # This square contains at least floor(s / v_y) * floor(s / u_x)
    # points.
    # Thus, floor(s / v_y) * floor(s / u_x) < n
    # => (s / v_y - 1) * (s / u_x - 1) < n
    # => s^2 - s*v_y - s*u_x + v_y*u_x < n
    # => s < (v_y + u_x) / 2 + sqrt(((v_y + u_x) / 2)^2 + n - v_y * u_x
    # => k_u = (v_y + u_x) / sqrt(2) + sqrt((v_y + u_x)^2 / 2 + 2 n - 2 v_y*u_x)
    k_u = int((v[1] + 1)/math.sqrt(2) + math.sqrt((v[1] + 1)**2 / 2 +
                                                        2 * n -
                                                        2 * v[1]))
    # Likewise for k_v but with s = k_v * |v| / sqrt(2) gives
    # k_v = 1/|v| * ((v_y + u_x) / sqrt(2) + sqrt((v_y + u_x)^2 / 2 + 2 n - 2 v_y*u_x))
    k_v = int(((v[1] + 1)/math.sqrt(2) +
                    math.sqrt((v[1] + 1)**2 / 2 +
                              2 * n - 2 * v[1])) / np.linalg.norm(v))

    print(f"k_u: {k_u}, k_v: {k_v}")

    # To avoid considering the same circle multiple times we'll
    # keep a dictionary with all points seen so far as key and
    # sets of the points thet they have been on the edge of a
    # circle with as the keys.
    seenWith = dict()

    # We can choose the first point to have positive u-coordinate 
    # and smaller absolute value of u-coordinate, or equal absolute 
    # value of u coordinate and smaller or equal absolute value of v-coordinate
    # The first because the circle centers stemming from the other case
    # are just minus the ones from this case. The second to avoid counting
    # the pair of points p, q AND q, p, which of course yields the same circle.
    for a1 in range(k_u + 1):
        for b1 in range(-k_v, k_v+1):

            p = a1 * u + b1 * v
            itera2 = list(range(-k_u, -a1 + 1)) + list(range(a1, k_u+1))

            for a2 in itera2:

                if abs(a2)==a1:
                    iterb2 = list(range(-k_v, -abs(b1) + 1)) + list(range(abs(b1), k_v+1)) 
                else:
                    iterb2 = range(-k_v, k_v+1)

                for b2 in iterb2:

                    #if a1==1 and b1==-1 and a2==-1 and b2==-2:
                        #breakpoint()

                    # If this q has already been in a circle with p
                    # The center has already been added and there is no need
                    # to do it again.
                    if (a1, b1) in seenWith and (a2, b2) in seenWith[(a1, b1)]:
                        continue
                    q = a2 * u + b2 * v

                    # If p and q are parallel, that would be no circle
                    if np.abs(np.cross(p, q)) < tol:
                        continue

                    # Calculate the circle center
                    center = findCircleCenter(p, q)

                    # If the circle obviously contains more than n points, skip it
                    # The largest rectangle that fits in a circle of radius r
                    # is a square with sidelength a = 2 * r / sqrt(2). 
                    # The smalles number of points inside such a square is
                    # floor(a / v_y) * floor(a / u_x) 
                    a = np.linalg.norm(center) * np.sqrt(2)
                    if np.floor(a / v[1]) * np.floor(a) > n:
                        continue

                    # Calculate the number of points inside the circle and 
                    # a list of the points on the edge of the circle.
                    noPointsInside, pointsOnEdge = latticePointsInCircle(u, v, center, tol)

                    # If the number of points inside is greater than n, skip it
                    if noPointsInside + 1 > n:
                        continue


                    # For each point on the edge, at all the other points
                    # to it's set of points it's been seen with.
                    # Also add the negative versions.
                    for point1 in pointsOnEdge:

                        mpoint1 = (-point1[0], -point1[1])

                        if point1 not in seenWith:
                            seenWith[point1] = set()
                            seenWith[mpoint1] = set()

                        for point2 in pointsOnEdge:
                            if point1 != point2:
                                mpoint2 = (-point2[0], -point2[1])
                                seenWith[point1].add(point2)
                                seenWith[mpoint1].add(mpoint2)


                    # Now we get to the good stuff :)
                    # The center of a circle containing k points with j points on the edge
                    # is an outer corner of the k:th, k+1:th, ..., k+j-1:th zone
                    # and an inner corner of the k+j:th zone.
                    for zone in range(noPointsInside + 1, noPointsInside + len(pointsOnEdge)):
                        if zone in bzones:
                            bzones[zone].append(center)
                            bzones[zone].append(-center)
                        else:
                            bzones[zone] = [center, -center]

    # bzones are popluated, now we sort them.
    for zone in bzones:
        bzones[zone].sort(key=angle)

    # And return the dict
    return bzones

def latticePointsInCircle(u, v, center, tol=1e-8):
    """
    Return the number of lattice points inside the circle
    and a list of the lattice coordinates of the points on the edge
    """
    # First, find the smallest parallellogram enclosing the circle
    radius = np.linalg.norm(center)

    # Calculate a bounds
    # The distance from the line a * u + lambda * v to the center point
    # is | v x (a * u - center) | / |v|
    # where x denotes crossproduct.
    # We thus want to solve
    # | v x (a * u - center) | / |v| < radius
    # for the smallest (in abs value) integer a 
    # | a * v x u - v x center | < radius * |v|
    # If rightside positive
    # a * v x u - v x center < radius * |v|
    # a * v x u < radius * |v| + v x center
    # a > (radius * |v| + v x center) / v x u     
    #   inequality switch as v x u is negative
    # If rightside negative
    # v x center - a * v x u < radius * |v|
    # a * v x u > v x center - radiuss * |v|
    # a < (v x center - radius * |v|) / v x u
    # Once again, the inequality switches because v x u is negative
    aLo = (np.cross(v, center) + (radius + tol) * np.linalg.norm(v)) / np.cross(v, u)
    aHi = (np.cross(v, center) - (radius + tol) * np.linalg.norm(v)) / np.cross(v, u)

    # Calculate b bounds
    # The distance from the line b * v + lambda * u to the center point
    # is | u x (b * v - center) | / |u|
    # Like before
    # | u x (b * v - center) | < radius * |u|
    # If positive
    # b u x v < radius * |u| + u x center
    # b < (radius * |u| + v x center) / u x v
    # If negative
    # b u x v > (u x center - radius * |u|)
    # b > (u x center - radius * |u|) / u x v
    bLo = (np.cross(u, center) - (radius + tol) * np.linalg.norm(u)) / np.cross(u, v)
    bHi = (np.cross(u, center) + (radius + tol) * np.linalg.norm(u)) / np.cross(u, v)

    # These are the most extreme bounds. Since the lattice coordinates are integers,
    # We can ceil / floor them accordingly. 
    aLo = math.ceil(aLo)
    bLo = math.ceil(bLo)
    aHi = math.floor(aHi)
    bHi = math.floor(bHi)

    # For each latticePoint in this parallellogram, check if it is inside the circle
    # If so, increment insideCount
    # If instead it is on the boundry, increment edgeCount
    insideCount = 0
    edgePoints = []
    for a in range(aLo, aHi + 1):
        for b in range(bLo, bHi + 1):
            d = np.linalg.norm(a * u + b * v - center)
            if d < radius - tol:
                insideCount += 1
            elif abs(d - radius) < tol:
                edgePoints.append((a, b))

    return insideCount, edgePoints


def findCircleCenter(p, q):
    """
    Returns the coordinates of the center of the circle going through the origin
    and points p and q as a np array
    """
    A = 2 * np.row_stack((p,q))
    b = np.array([np.dot(p, p), np.dot(q, q)])
    return np.linalg.solve(A, b)

def angle(point):
    """
    Returns the angle the point makes with the positive x-axis
    """
    return np.arctan2(point[1], point[0])
