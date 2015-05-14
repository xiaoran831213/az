## brain surface
import vertex
import sys

class Sample:
    """ MRI sample
    """
    def __init__(self, phe, vxl, vxr)
        """ initializer
        phe:    the clinical phenotype data of the sample in one line of text
        vxl:    the CSV file object of left hemisphere vertices
        vxr:    the CSV file object of right hemisphere vertices
        """
        self.

        class Clinic:
    """ clinical data
    """
    ## Clinical phenotype key - type tuples. Columns matching these keys will be
    ## extracted from CSV table and parsed to corresponding types
    KEY_PHE = []
    def __init__(self, row)
        """ initializer
        row: the clinical phenotype data of a sample in one line of text,
        from which data will be picked and parsed according to key - type
        pairs stored in KEY_PHE
        """
    


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
    
    t1 = tarfile.open('dat/t1.tar.gz')
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
