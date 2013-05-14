class dyn_data(object):
    simu = [[]]
    times = []
    objects = []
    t = 0
    tstep = float()

# static data:
def force(r):
    - q²/(4*math.pi*Epsilon) * r**-2

obj = pub_data.objects[random.random(0, len(pub_data.objects))]

class point(object):
    def __init__(self, *args):
        self.coords=args
        self.dim = len(args)
    def __add__(self, pt):
        if self.dim != pt.dim:
            raise ValueError
        ret = []
        for i, e in enumerate(self.coords):
            ret.append(e + pt.coords[i])
        return point(*ret)
    def __div__(self, num):
        for e in self.coords:
            ret = []
            ret.append(e/float(num))
        return point(*ret)
    def __mul__(self, k):
        ret = []
        for e in self.coords:
            ret.append(e * k)
        return point(*ret)

import math

class vector(point):
    def __mul__(self, q):
        if type(p) == vector:
            ret=0
            for i, e in enumerate(self.coords):
                ret += e * q.coords[i]
            return ret
        else:
            return super(vector, self).__mul__(q)
    def norm(self):
        return math.sqrt(self * self)


























