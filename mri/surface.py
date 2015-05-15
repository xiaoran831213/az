## brain surface
from vertex import Vertex
import sys
import csv

class Surface(list):
    """ brain surface

    An list of vertices sorted by position
    max_xyz: maximum vertex position
    min_xyz: minimum vertex position
    sqn: sample sequence number
    """
    def __init__(self, sqn, lcsv, rcsv):
        """ initializer
        lcsv: file opened on the L.himersphere surface vertices csv file
        rcsv: file opened on the R.himersphere surface vertices csv file
        """
        self.sqn = sqn
        self.min_xyz = None
        self.max_xyz = None

        ## For left himersphere
        ## read header line, then match vertex phenotype key names with
        ## the header.
        ## cols: pairs of (column index, python type) describing vertex
        ## phenotype columns
        ## ipos: index of the vertex position column
        lcsv = csv.reader(lcsv)
        head = lcsv.next()
        cols = [(head.index(k), t) for k, t in Vertex.KEY_PHE]
        ipos = head.index(Vertex.KEY_POS)

        ## iterate through rest of the csv, create one vertex per line;
        for line in lcsv:
            phe = [t(line[i]) for i, t in cols]   # phenotypes
            pos = tuple(eval(line[ipos]))         # the position
            self.append(Vertex(phe, pos))     # new Vertex

        ## Do the same for right himersphere, though the header line is
        ## skipped without processing.
        rcsv = csv.reader(rcsv)
        rcsv.next()             # skip header
        for line in rcsv:
            phe = [t(line[i]) for i, t in cols]   # phenotypes
            pos = tuple(eval(line[ipos]))         # the position
            self.append(Vertex(phe, pos))     # new Vertex

        ## record boundary on 3 dimension each time a vertex is created.
        min_x = min_y = min_z = sys.float_info.max
        max_x = max_y = max_z = -sys.float_info.max
        for pos in (v.pos for v in self):
            min_x = min(min_x, pos[0])
            min_y = min(min_y, pos[1])
            min_z = min(min_z, pos[2])
            max_x = max(max_x, pos[0])
            max_y = max(max_y, pos[1])
            max_z = max(max_z, pos[2])
        self.min_xyz = (min_x, min_y, min_z)
        self.max_xyz = (max_x, max_y, max_z)
        self.sort(key = lambda v: v.pos)

    def get_pos(self, fr = None, to = None):
        """ return a shallow copy of vertex positions """
        return [v.pos for v in self[fr : to]]

    def get_vtx(self, fr = None, to = None):
        """ return a shallow copy of vertices """
        return self[fr:to]

if __name__ == "__main__":
    import vertex
    reload(vertex)

    with open("vertex.test.L.csv") as lcsv:
        with open("vertex.test.R.csv") as rcsv:
            sur = Surface('0000', lcsv, rcsv)
