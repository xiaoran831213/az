## the whole study data set
import tarfile
import os

def tar2csv(src, dst = None):
    """ extract MRI scanned brain surface vertices from source
    tarball.
    The surface is stored in CSV table, one line per vertex.

    src: file opened on the source tarball with surface vertex
    data (optionally gzipped), or the filename of the tarball.
    
    dst: destination directory to place extracted csv files
    """
    if not dst:
        dst = "csv"
    
    if not os.path.exists(dst):
        os.mkdir(dst)

    ## open tarball
    if isinstance(src, file):
        tf = tarfile.open(fileobj = src)
    elif isinstance(src, str):
        tf = tarfile.open(name = src)
    else:
        raise IOError(src, 'not a file object nor path')

    ## this flag denote the encounter of left surface csv,
    ## which means the next csv file in the tarball is the
    ## right surface of the same subject
    lcsv = False

    # fi represent files and directories packed in the tarball
    for fi in tf:               
        if not fi.path.endswith('vertices.csv'):
            continue

        ex = tf.extractfile(fi)          # csv to extract
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
    else:
        tf.close()

def test():
    import sys
    import time

    start = time.time()
    if len(sys.argv) > 1:
        src = sys.argv[1]
    else:
        src = "dat/s10.tar.gz"

    if len(sys.argv) > 2:
        dst = sys.argv[2]
    else:
        dst = None

    tar2csv(tar, dst)
    print time.time() - start
    
if __name__ == "__main__":
    test()
