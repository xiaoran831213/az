## the whole study data set
import csv
import tarfile
import os
import cPickle
import gzip

class Scan(list):
    """ all MRI scaned surfaces in a study

    mir.Scan is a collection of mri.Surface
    """
    def __init__(self, src):
        """ initializer

        Load scanned surfaces from source directory.
        src: the source directoy containing surface csv files
        """
        pass

def extract(tar, dst = "ext"):
    """ extract surface files from source

    tar: original tar.gz from data provider
    dst: destination directory to place extracted files
    """
    if not os.path.exists(dst):
        os.mkdir(dst)

    ## this flag denote the encounter of left surface csv,
    ## which means the next file is the right surface csv
    ## of the same MRI experiment subject
    lcsv = False

    tar = tarfile.open(fileobj = tar)
    for fi in tar:
        if not fi.path.endswith('vertices.csv'):
            continue

        ex = tar.extractfile(fi)          # csv to extract
        if not lcsv:
            sn = fi.path.split(os.path.sep)[1] # serial number
            fo = open(os.path.join(dst, sn + ".csv"), "wb")
            fo.writelines(ex)
            lcsv = True
        else:
            ex.readline()                 # skip header
            fo.writelines(ex)
            fo.close()
            lcsv = False
            print fo
        ex.close()
    tar.close()

if __name__ == "__main__":
    with open("vertices.test.tar.gz") as tgz:
        extract(tgz)
        
