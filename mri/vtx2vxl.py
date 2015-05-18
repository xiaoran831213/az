## convert floating point vertex postion to integer voxels

import dataset
import surface

## field index of vertex position

IX, IY, IZ = [i for i, c, t in surface.VTX_COL if c in 'xyz']
class Voxface(list):
    """ 3D surface represented by voxels

    Voxel are fix sized volumns in 3D space.
    Voxel based surface is created from vertex based surface.
    surface vertices fall in the same voxel is to be aggregated
    later
    """

    def __init__(self, sf, sz = 1):
        """ create voxel based surface by algin vertex in floating
        point coordinates to integer voxel cells

        sf: the raw surface to be aligned
        sz: size of a voxel in the raw surface measure, examples:
            sz =  1: align to rinteger, 78.6 -> 79, 65.4 -> 65
            sz =  2: align to double, 78.6 -> 39, 79.6 -> 40
            sz = .5: align to half, 78.6 -> 157, 65.4 -> 131
        by default sz is 1
        """

    
## align surface 3D coordinate to voxel
def vtx2vxl(sf, sz = 1):
    """ algin float vertex position to integer voxcel coordinate

    sf: the surface to operate on
    sz: size of a voxel in the original vertex masure, some examples
    are:
        sz =  1: align to rinteger, 78.6 -> 79, 65.4 -> 65
        sz =  2: align to double, 78.6 -> 39, 79.6 -> 40
        sz = .5: align to half, 78.6 -> 157, 65.4 -> 131
    by default sz is 1
    """
    for v in sf:
        v[IX] = int(round(v[IX]/sz))
        v[IY] = int(round(v[IY]/sz))
        v[IZ] = int(round(v[IZ]/sz))

## unit test
def test():
    import sys
    import time
    reload(dataset)
    reload(surface)

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
    vtx2vxl(sf0, 2)
    end = time.time()
    print sf0
    print "alignment time = ", end - start, "\n"

    ## try vocel alignment on the dataset
    print 'align voxel for ds0 and save as ds1'
    start = time.time()
    ds1 = dataset.make(dst)
    dataset.apply(vtx2vxl, ds0, ds1, sz = 2)
    end = time.time()
    print "total time = ", end - start, "\n"

    print "compre the coordinates before and after voxel",
    print "alignment"
    print "surface #1 before alignment:"
    print ds0.fetch(1)
    print "after alignment:"
    print ds1.fetch(1)


## unit test invocation
if __name__ == "__main__":
    test()
