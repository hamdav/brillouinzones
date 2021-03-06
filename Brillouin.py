import numpy as np
import math

class BZ:
    def __init__(self, v2=[0, 1], n=3):
        """
        v2 is the second basis vector. The first is always [1, 0]

        We assume that v2 is such that there was no smaller v2 that would result in the same lattice.
        """
        self.tol = 1e-8
        self.v1 = np.array([1, 0])
        self.v2 = np.array(v2)
        self.toLatticeCoords = np.linalg.inv(np.column_stack((self.v1, self.v2)))
        self.bzones = dict()
        self.n = n

        self.latticePoints = []
        for i in range(-self.n, self.n+1):
            for j in range(-self.n, self.n+1):
                if not i == j == 0:
                    self.latticePoints.append(i * self.v1 + j * self.v2)

        print(f"Calculated latticePoints: length {len(self.latticePoints)}")
        # 
        circleCenters = []
        seenCenters = set()
        for i, p in enumerate(self.latticePoints[:-1]):
            if i % 100 == 0:
                print(i)
            for q in self.latticePoints[i+1:]:
                # If q is a multiple of p
                if np.abs(np.cross(p, q)) < self.tol:
                    continue

                center = self.findCircleCenter(p, q)

                # If the circle obviously contains more than n points
                # The largest rectangle that fits in a circle of radius r
                # is a square with sidelength a = 2 * r / sqrt(2). 
                # The smalles number of points inside such a rectangle is
                # floor(a / v2_y) * floor(a / v1_x) 
                a = np.linalg.norm(center) * np.sqrt(2)
                if np.floor(a / self.v2[1]) * np.floor(a) > self.n:
                    continue

                centerTuple = tuple(np.around(center, decimals=8))
                if centerTuple in seenCenters:
                    continue
                #for seenCenter in circleCenters:
                    #if np.allclose(seenCenter, center):
                        #break

                else:
                    seenCenters.add(centerTuple)
                    seenCenters.add((-centerTuple[0], -centerTuple[1]))
                    circleCenters.append(center)
                    circleCenters.append(-center)

        print(len(circleCenters))
        #for center in filter(lambda c: np.linalg.norm(c) < self.n/2, circleCenters):
        for center in circleCenters:
            pointsInside, pointsOnEdge = self.latticePointsInCircle(center)
            for zone in range(pointsInside + 1, pointsInside + pointsOnEdge):
                if zone in self.bzones:
                    self.bzones[zone].append(center)
                else:
                    self.bzones[zone] = [center]

    def latticePointsInCircle(self, center):
        """
        Return the number of lattice points inside the circle
        and the number of lattice points on the edge
        """
        # First, find the smallest parallellogram enclosing the circle
        radius = np.linalg.norm(center)
        if radius > self.n:
            return -1

        # Calculate a1 bounds
        # The distance from the line a1 * v1 + lambda * v2 to the center point
        # is | v2 x (a1 * v1 - center) | / |v2|
        # where x denotes crossproduct.
        # We thus want to solve
        # | v2 x (a1 * v1 - center) | / |v2| < radius
        # for the smallest (in abs value) integer a1 
        # | a1 * v2 x v1 - v2 x center | < radius * |v2|
        # If rightside positive
        # a1 * v2 x v1 - v2 x center < radius * |v2|
        # a1 * v2 x v1 < radius * |v2| + v2 x center
        # a1 > (radius * |v2| + v2 x center) / v2 x v1     
        #   inequality switch as v2 x v1 is negative
        # If rightside negative
        # v2 x center - a1 * v2 x v1 < radius * |v2|
        # a1 * v2 x v1 > v2 x center - radiuss * |v2|
        # a1 < (v2 x center - radius * |v2|) / v2 x v1
        # Once again, the inequality switches because v2 x v1 is negative
        a1Lo = (np.cross(self.v2, center) + (radius + self.tol) * np.linalg.norm(self.v2)) / np.cross(self.v2, self.v1)
        a1Hi = (np.cross(self.v2, center) - (radius + self.tol) * np.linalg.norm(self.v2)) / np.cross(self.v2, self.v1)

        # Calculate a1 bounds
        # The distance from the line a2 * v2 + lambda * v1 to the center point
        # is | v1 x (a2 * v2 - center) | / |v1|
        # Like before
        # | v1 x (a2 * v2 - center) | < radius * |v1|
        # If positive
        # a2 v1 x v2 < radius * |v1| + v1 x center
        # a2 < (radius * |v1| + v2 x center) / v1 x v2
        # If negative
        # a2 v1 x v2 > (v1 x center - radius * |v1|)
        # a2 > (v1 x center - radius * |v1|) / v1 x v2
        a2Lo = (np.cross(self.v1, center) - (radius + self.tol) * np.linalg.norm(self.v1)) / np.cross(self.v1, self.v2)
        a2Hi = (np.cross(self.v1, center) + (radius + self.tol) * np.linalg.norm(self.v1)) / np.cross(self.v1, self.v2)

        # These are the most extreme bounds. Since the lattice coordinates are integers,
        # We can ceil / floor them accordingly. 
        a1Lo = math.ceil(a1Lo)
        a2Lo = math.ceil(a2Lo)
        a1Hi = math.floor(a1Hi)
        a2Hi = math.floor(a2Hi)

        # For each latticePoint in this parallellogram, check if it is inside the circle
        # If so, increment insideCount
        # If instead it is on the boundry, increment edgeCount
        insideCount = 0
        edgeCount = 0
        for a1 in range(a1Lo, a1Hi + 1):
            for a2 in range(a2Lo, a2Hi + 1):
                d = np.linalg.norm(a1 * self.v1 + a2 * self.v2 - center)
                if d < radius - self.tol:
                    insideCount += 1
                elif abs(d - radius) < self.tol:
                    edgeCount += 1

        return insideCount, edgeCount

        
    def findCircleCenter(self, p, q):
        A = 2 * np.row_stack((p,q))
        b = np.array([np.dot(p, p), np.dot(q, q)])
        return np.linalg.solve(A, b)


    def getBrillouinZone(self, n):
        return sorted(self.bzones[n], key=angle)

def angle(point):
    return np.arctan2(point[1], point[0])
