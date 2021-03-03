import numpy as np

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
                # floor(a / v2_y) * floor(a / v1_x) maybe - 4 because the corners
                # shouldn't count
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

        # Calculate a1 lower bound
        d = 0
        a1Lo = 0
        while abs(d) <= radius + self.tol:
            d = np.cross(self.v2, a1Lo * self.v1 - center) / np.linalg.norm(self.v2)
            a1Lo -= 1
        # Calculate a1 high bound
        d = 0
        a1Hi = 0
        while abs(d) <= radius + self.tol:
            d = np.cross(self.v2, a1Hi * self.v1 - center) / np.linalg.norm(self.v2)
            a1Hi += 1
        # Calculate a2 lower bound
        d = 0
        a2Lo = 0
        while abs(d) <= radius + self.tol:
            d = np.cross(self.v1, a2Lo * self.v2 - center) / np.linalg.norm(self.v1)
            a2Lo -= 1
        # Calculate a2 high bound
        d = 0
        a2Hi = 0
        while abs(d) <= radius + self.tol:
            d = np.cross(self.v1, a2Hi * self.v2 - center) / np.linalg.norm(self.v1)
            a2Hi += 1

        # For each latticePoint in this parallellogram, check if it is inside the circle
        # If so, increment count
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
