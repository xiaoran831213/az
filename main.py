import csv
import os
import mri.tar2csv as tar2csv
import mri.work as work

def test():
    work.extract_region('dat/npy', 2010)
    pass
    # work.csv2npy('dat/csv', 'dat/npy', ovr = 0)
    # work.vtx2grd('dat/npy', 'dat/grd', ovr = 0, sz = 1)
    # work.srt_pos('dat/grd', 'dat/srt', ovr = 0)
    # work.cmb_pos('dat/srt', 'dat/cmb', ovr = 0)
    # work.sfr2vlm('dat/cmb', 'dat/vlm', ovr = 0, dim = (64,)*3)

def main():
    pass
    # tar2csv.tar2csv('raw/vtx.tar.gz', 'dat/csv')
    
if __name__ == "__main__":
    test()
