## the whole study data set
import tarfile
import os

def extract(tar, dst = ''):
    """ extract surface files from source

    tar: original tar.gz from data provider
    dst: destination directory to place extracted files
    """
    if not dst:
        dst = "ext"
    
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
    import sys
    if len(sys.argv) > 1:
        src = sys.argv[1]
    else:
        src = "s10.tar.gz"

    if len(sys.argv) > 2:
        dst = sys.argv[2]
    else:
        dst = None

    with open(src) as tgz:
        extract(tgz, dst)
