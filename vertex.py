## class vertex

class Vertex:
    """ record one vertex of brain surface
    """

    def __init__(self, rc):
        """ initializer
        rc is the record of one vertex in comma splited text originally read
        from raw csv file, which has fields in the following order:

            id: index of the vertex (skip)
            lb: label of 56 anotomical region
            sulc:          is this vertex located in a sulcus
            ps:     3-dimensional position of this vertex, a 
            area:           area of nearby surface
            curv:      mean curvature of the nearby surface
            tdep:            travel depth of the nearby surface
            gdep:            geodesic depth of the nearby surface
            cnvx:            freesurfer_convexity
            thck:            freesurfer_thickness
        """
        self.id = int(rc[0])
        self.lb = int(rc[1])
        self.sulc = int(rc[2])
        self.ps = tuple(eval(rc[3]))
        self.area = float(rc[4])
        self.curv = float(rc[5])
        self.tdep = float(rc[6])
        self.gdep = float(rc[7])
        self.cnvx = float(rc[8])
        self.thck = float(rc[9])

    def __str__(self):
        s = "{:6d} {:4d} {:3d}".format(self.id, self.lb, self.sulc)
        s += "[{:7.4f},{:7.4f},{:7.4f}]".format(*self.ps)
        s += "{:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f} {:7.4f}".format(
            self.area, self.curv, self.tdep,
            self.gdep, self.cnvx, self.thck)
        return s
