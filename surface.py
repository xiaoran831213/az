## brain surface
import vertex
import sys

class Surface:
    """ brain surface

    An list of vertices sorted by coordinates
    """
    def __init__(self):
        self.vertices = []
        self.max_cord = (
            sys.float_info.max, sys.float_info.max, sys.float_info.max)
        self.min_cord = (
            sys.float_info.min, sys.float_info.min, sys.float_info.min)

    def add(self, vtx):
        """ add one vertex to the surface """
        self.vertices.append(vtx)


if __name__ == "__main__":
    import csv
    import tarfile
    import os
    import vertex
    reload(vertex)
    
    wd = os.getcwd()
    t1 = tarfile.open('../dat/t1.tar.gz')
    f1 = t1.extractfile('d1.csv')
    c1 = csv.reader(f1)

    s1 = Surface();
    c1.next()                   # skip header
    for row in c1:
        v = vertex.Vertex(row)
        s1.add(v)

    v1 = s1.vertices[1]
    f1.close()
    t1.close()

    os.chdir(wd)
