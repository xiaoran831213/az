## class vertex

class Vertex:
    """ one vertex of brain surface
    """
    ## Vertex phenotype key - type tuples. Columns matching these keys will be
    ## extracted from CSV table and parsed to corresponding types
    KEY_PHE = [
        ("Index", int),
        ("Label", int),
        ("Sulcus", int),
        ("area", float),
        ("mean curvature", float),
        ("travel depth", float),
        ("geodesic depth", float),
        ("FreeSurfer convexity", float),
        ("FreeSurfer thickness", float)]

    ## Vertex position key. Column mataching this will be extracted as
    ## 3D postion of the vertex. The type is [float, float, float]
    KEY_POS = "coordinates"

    def __init__(self, phe, pos):
        """ initializer
        row: one line of text recording one vertex read from csv file,
        from which phenotype data will be picked and parsed according to
        key - type pairs stored in KEY_PHE, and location of the vertex
        be extracted according to KEY_POS
        """
        self.phe = phe
        self.pos = pos

    def __str__(self):
        mapping = zip([k for k, t in Vertex.KEY_PHE], phe);
        mapping.append((Vertex.KEY_POS, pos))
        return repr(self)
        
    def __repr__(self):
        mapping = zip([k for k, t in Vertex.KEY_PHE], self.phe);
        mapping.append((Vertex.KEY_POS, self.pos))
        return repr(mapping)


        
if __name__ == "__main__":
    import csv

    pcsv = "vertex.test.csv";
    with open(pcsv) as fcsv:
        reader = csv.reader(fcsv)
        head = reader.next()              # head line

        ## get column index - column type pair
        cols = [(head.index(k), t) for k, t in Vertex.KEY_PHE]

        ## index of vertex position column
        ipos = head.index(Vertex.KEY_POS)

        ## create a few vertices and print them
        c = 0
        for line in reader:
            if c > 10:
                break
            c += 1
            phe = [t(line[i]) for i, t in cols]
            pos = tuple(eval(line[ipos]))
            vtx = Vertex(phe, pos)
            print vtx

