
class Graticule:

    def __init__(self, point1, point2, _id):

        self.p0 = point1
        self.p1 = point2
        self.id = _id
        self.p0_intersect_id = None
        self.p1_intersect_id = None
        self.poly_id = None

    # oh rite change directly
    #
    #def set_p0_intersect_id(self, _id):
    #    self.p0_intersect_id = _id

    #def set_p1_intersect_id(self, _id):
    #    self.p1_intersect_id = _id


# id for this one is y
class Lat( Graticule ):

    pass 

# id for this one is x
class Lon(Graticule):

    pass

# id for this one is x - y
class Diag( Graticule ):

    pass

class EdgeSlice():
    def __init__(self):
        self.p0 = (0.0, 0.0) 
        self.p1 = (0.0, 0.0)
        self.clipped = False

    def __init__(self, p0, p1):
        self.p0 = p0 
        self.p1 = p1

        self.clip_start = 0.0
        self.clip_end = 1.0
        self.clipped = False

    def set_start_end(self, s, e):

        assert(s >= 0 and s < 1)
        assert(e > 0 and e <= 1)
        self.clip_start = s
        self.clip_end = e
        self.clipped = True