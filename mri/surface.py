## brain surface
from vertex import Vertex
import sys

class Surface:
    """ brain surface

    An list of vertices sorted by position
    max_pos: maximum vertex position
    min_pos: minimum vertex position
    """
    def __init__(self):
        self.vertices = []
        self.min_pos = (
            sys.float_info.max, sys.float_info.max, sys.float_info.max)
        self.max_pos = (
            sys.float_info.min, sys.float_info.min, sys.float_info.min)

    def add_vtx(self, vtx):
        """ add one vertex to the surface """
        self.vertices.append(vtx)

    def sort_by_pos(self):
        self.vertices.sort(key = lambda vex: vex.pos)

    def span_by_pos(self):
        lambda p, q: [(min(a, b), max(a,b)) for p, q in zip(p, q)]

if __name__ == "__main__":
    import csv
    import tarfile
    import os
    import vertex
    reload(vertex)

    sur = Surface()
    with open("vertex.test.csv") as fcsv:
        reader = csv.reader(fcsv)
        head = reader.next()              # head line

        ## get pairs of column index - column type for vertex phenotypes
        cols = [(head.index(k), t) for k, t in Vertex.KEY_PHE]

        ## index of vertex position column
        ipos = head.index(Vertex.KEY_POS)

        ## iterate through lines of the csv, create one vertex per line
        for line in reader:
            phe = [t(line[i]) for i, t in cols]
            pos = tuple(eval(line[ipos]))
            vtx = Vertex(phe, pos)
            sur.add_vtx(vtx)
