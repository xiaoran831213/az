## the whole study data set
import csv
import tarfile
import os
import cPickle
import gzip
import glob
import surface

class MRIS(list):
    """ all MRI scaned surfaces in a study

    mir.Scan is a collection of mri.Surface
    """

    def __init__(self, src):
        """ initializer

        Load scanned surfaces from source directory.
        src: the source directoy containing surface csv files
        """
        src = os.path.join(src, "*")
        for fi in glob.glob(src):
            sn = os.path.splitext(os.path.basename(fi))[0]
            fp = open(fi, 'rb')
            self.append(surface.Surface(sn, fp))
            fp.close()
                                        
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        src = sys.argv[1]
    else:
        src = "ext"

    if len(sys.argv) > 2:
        dst = sys.argv[2]
    else:
        dst = None

##    mris = MRIS(src)
