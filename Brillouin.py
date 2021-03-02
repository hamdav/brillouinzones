import numpy as np

class BZ:
    def __init__(self, v2=[0, 1]):
        """
        v2 is the second basis vector. The first is always [1, 0]

        We assume that v2 is such that there was no smaller v2 that would result in the same lattice.
        """
        self.v1 = np.array([1, 0])
        self.v2 = np.array(v2)
        self.bzones = dict()
        
        self.calculate1BZ()

    def lineFromLatticePoint(self, a1, a2):
        # The line is defined as begin the points
        # equidistant from the origin and the lattice point
        point = (a1 * self.v1 + a2 * self.v2) / 2

        # vec is orthogonal to point
        vec = np.array([-point[1], point[0]])

        return Line(point, vec)

    def calculate1BZ(self):
        # While calculating 1BZ we need only consider lines from the points
        # that have at most one step in any direction
        # I.e. points (0, 1), (1, 0), (1, 1), (-1, 0), ...
        # This is because of the assumption that v2 is the smallest equivalent vector. 
        self.lines = [self.lineFromLatticePoint(0, 1),
                      self.lineFromLatticePoint(1, 1),
                      self.lineFromLatticePoint(1, 0),
                      self.lineFromLatticePoint(1, -1),
                      self.lineFromLatticePoint(0, -1),
                      self.lineFromLatticePoint(-1, -1),
                      self.lineFromLatticePoint(-1, 0),
                      self.lineFromLatticePoint(-1, 1)]


        # The point 1/2 * v2 on the line from the point v2 (e.g. (0, 1))
        # will always be on the boundry of 1BZ because of the assumtion that
        # v2 is the smallest equivalent vector.
        p0 = self.v2 / 2
        l0 = self.lines[0]
        p = p0
        l = l0
        self.bzones[1] = [p0]

        while True:
            p, l = getNextIntersect(l, p, self.lines)

            if p is None:
                raise RuntimeError
            elif l == l0:
                break
            else:
                self.bzones[1].append(p)
        

    def populateBzones(n):
        pass


    def getBrillouinZone(n):
        pass

class Line:
    def __init__(self, point, vector):
        self.point = np.array(point)
        self.vector = np.array(vector)

        
def getNextIntersect(line, point, lines):
    """
    Imagine you are hanging out on @line at point @point.
    Then you start wandering along @line in whichever direction is clockwise
    around the origin. This function returns the first intersection of @line
    with any line in @lines that you encounter along with the line it intersects.
    If there is none, returns None

    @line is a Line,
    @point is a np.array([x, y])
    @lines is a list of Lines
    @returns tuple(np.array([x, y]), Line)
    """
    
    # Calculate all intersections
    # The intersection of line1 and line2 is given by the equation
    intersections = [(getIntersect(line, l), l) for l in lines if getIntersect(line, l) is not None]

    # The intersection point that is the closest to p but still clockwise
    # is the one that has the largest cross product with p that is negative
    crossProducts = [(np.cross(point, ipoint), ipoint, l) for ipoint, l in intersections]
    candidates = filter(lambda x: x[0] < 0, crossProducts)

    if not candidates:
        return (None, None)

    closestIntersect = max(candidates, key=lambda x: x[0])

    return (closestIntersect[1], closestIntersect[2])


def getIntersect(line1, line2):
    """
    Returns intersection of line1 and line2 on the form np.array([x, y])
    if line1 == line2, return None
    """
    if line1 == line2:
        return None

    try:
        c1, c2 = np.linalg.solve(np.column_stack((line1.vector, line2.vector)), line1.point - line2.point)
    except np.linalg.LinAlgError:
        return None

    return line1.vector * c1 + line1.point
