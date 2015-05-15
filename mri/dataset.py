## the whole study data set
import csv
import tarfile

class DataSet(list):
    """ all MRI scaned surfaces in a study
    """

    def __init__(self, tar):
        """ initializer

        Load sample surfaces from tarball.
        tar: the file object open on the tar ball containing sample
        surfaces created from MRI scanning data
        """
        with tarfile.open(fileobj = tar) as tar:
            
        



