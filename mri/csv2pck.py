## create surface object from csv tables
## save surfaces into python pickles
import surface
import os
import glob

def csv2pck(src = None, dst = None):
    """ extract surface files from source

    src: source directory containning table of surface vertices in
    csv format, the default is 'ext', which is the default destination
    of tar ball extraction
    dst: destination directory to save created Surface objects in
    python pickle format, the default is 'pck'
    """

    if not src:
        src = 'ext'
    if not dst:
        dst = 'pck'
    
    if not os.path.exists(dst):
        os.mkdir(dst)

    for p in glob.glob(os.path.join(src, '*')):
        sf = surface.from_csv(p)
        surface.save_pk(sf)

def test():
    import sys
    import time
    reload(surface)
    
    if len(sys.argv) > 1:
        src = sys.argv[1]
    else:
        src = None

    if len(sys.argv) > 2:
        dst = sys.argv[2]
    else:
        dst = None

    start = time.time()
    csv2pck(src, dst)
    end = time.time()
    print end - start

if __name__ == "__main__":
    test()
