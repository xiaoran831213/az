import numpy as np
## csv to surface
CSV = [
    (0, "Index", int),
    (1, "Label", int),
    (2, "Sulcus", int),
    (3, "coordinates", lambda w: tuple(eval(w))),
    (4, "area", float),
    (5, "mean curvature", float),
    (6, "travel depth", float),
    (7, "geodesic depth", float),
    (8, "FreeSurfer convexity", float),
    (9, "FreeSurfer thickness", float)]

__F3D = np.dtype([('x', '<f4'), ('y', '<f4'), ('z', '<f4')])
__I3D = np.dtype([('x', '<u1'), ('y', '<u1'), ('z', '<u1')])

NPY = np.dtype([
    ('idx', '<i4'),
    ('lbl', '<i2'),
    ('slc', '<i1'),
    ('pos', __F3D),
    ('crv', '<f4'),
    ('are', '<f4'),
    ('tdp', '<f4'),
    ('gdp', '<f4'),
    ('cnv', '<f4'),
    ('tck', '<f4')])

VOX = np.dtype([
    ('idx', '<i4'),
    ('lbl', '<i2'),
    ('slc', '<i1'),
    ('pos', __I3D),
    ('crv', '<f4'),
    ('are', '<f4'),
    ('tdp', '<f4'),
    ('gdp', '<f4'),
    ('cnv', '<f4'),
    ('tck', '<f4')])

BND = np.dtype([
    ('ssn', np.str, 32),
    ('min', '<f4', (3,)),
    ('max', '<f4', (3,)),
    ('len', '<f4', (3,))])

VLM = np.dtype([
    ('lbl', '<i2'),
    ('slc', '<u1'),
    ('crv', '<f4'),
    ('are', '<f4'),
    ('tdp', '<f4'),
    ('gdp', '<f4'),
    ('cnv', '<f4'),
    ('tck', '<f4')])

D32 = np.array((32,)*3, dtype = __I3D)
D36 = np.array((36,)*3, dtype = __I3D)
D64 = np.array((64,)*3, dtype = __I3D)
D96 = np.array((96,)*3, dtype = __I3D)

