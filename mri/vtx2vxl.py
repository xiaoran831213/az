## convert floating point vertex postion to integer voxels

import dataset

## field index of vertex position
IDX_POS = [c for i, c, t in dataset.surface.VFD].index("coordinates")

## align surface 3D coordinate to voxel
def pos2vox(sf, sz = 0.5):
    """ algin float vertex position to integer voxcel coordinate

    sf: the surface to operate on
    sz: size of a voxel in powers of 2, some examples are:
        sz =  1: align to integer, 78.6 -> 79, 65.4 -> 65
        sz =  2: align to double, 78.6 -> 39, 79.6 -> 40
        sz = .5: align to half, 78.6 -> 157, 65.4 -> 131
    by default sz is 0.5
    """
    for i, v in enumerate(sf):
        v = list(v)
        vx, vy, vz = v[IDX_POS]
        v[IDX_POS] = (
            int(round(vx/sz)),
            int(round(vy/sz)),
            int(round(vz/sz)))
        
        sf[i] = tuple(v)

## unit test
def test():
    import sys
    import time
    reload(dataset)

    ## parse command line
    if len(sys.argv) > 1:
        src = sys.argv[1]
    else:
        tag = "dat/ds0"

    if len(sys.argv) > 2:
        dst = sys.argv[2]
    else:
        dst = "dat/vx0"

    ## load source dataset
    print 'load ds0'
    start = time.time()
    ds0 = dataset.make(tag)
    end = time.time()
    print "\n".join(ds0)
    print "loading time = ", end - start, "\n"

    ## fetch the first surface and try voxel alignment
    print 'align voxel for surface #0, fetching:',
    start = time.time()
    sf0 = ds0.fetch(0)
    end = time.time()
    print sf0
    print "fetching time = ", end - start, "\n"

    print 'perform voxel alignment:'
    start = time.time()
    pos2vox(sf0, 2)
    end = time.time()
    print sf0
    print "alignment time = ", end - start, "\n"

    ## try vocel alignment on the dataset
    print 'align voxel for ds0 and save as ds1'
    start = time.time()
    ds1 = dataset.make(dst)
    print ds1
    dataset.apply(pos2vox, ds0, ds1, sz = 0.5)
    end = time.time()
    print "total time = ", end - start, "\n"

    print "compre pre and post alignment coordinate",
    print "of surface #1"
    print ds0.fetch(1)
    print ds1.fetch(1)


## unit test invocation
if __name__ == "__main__":
    test()
