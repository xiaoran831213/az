## brain surface
from vertex import Vertex
import sys
import csv

class Surface1(list):
    """ brain surface

    An list of vertices sorted by position
    sn: sample 64bit integer serial number
    fp: file object pointing the opened surface vertex csv
    """
    def __init__(self, sn, fp):
        """ initializer
        fp: file opened on the surface vertex csv file
        """
        self.sn = sn

        ## read the header and match with vertex phenotype key names
        ## cols: pairs of (column index, python type) describing vertex
        ## phenotype columns
        ## ipos: index of the vertex position column
        fcsv = csv.reader(fp)
        head = fcsv.next()
        cols = [(head.index(k), t) for k, t in Vertex.KEY_PHE]
        ipos = head.index(Vertex.KEY_POS)

        ## iterate through rest of the csv, create one vertex per line;
        for line in fcsv:
            phe = [t(line[i]) for i, t in cols]   # phenotypes
            pos = tuple(eval(line[ipos]))         # the position
            self.append(Vertex(phe, pos))         # new Vertex

    def get_pos(self, fr = None, to = None):
        """ return a shallow copy of vertex positions """
        return [v.pos for v in self[fr : to]]

    def get_vtx(self, fr = None, to = None):
        """ return a shallow copy of vertices """
        return self[fr:to]

    def get_bnd(self, fr = None, to = None):
        """ return bounds on 3-dimensional
        vertex positions """
        min_x = min_y = min_z = sys.float_info.max
        max_x = max_y = max_z = -sys.float_info.max
        for pos in [v.pos for v in self[fr : to]]:
            min_x = min(min_x, pos[0])
            min_y = min(min_y, pos[1])
            min_z = min(min_z, pos[2])
            max_x = max(max_x, pos[0])
            max_y = max(max_y, pos[1])
            max_z = max(max_z, pos[2])
        return ((min_x, max_x), (min_y, max_y), (min_z, max_z))

if __name__ == "__main__":
    with open("vertex.test.csv") as f:
        sur = Surface1('0000', f)
