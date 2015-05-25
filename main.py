import csv
import os

def test():
    import mri.work as mw
    import mri.tar2csv as t2csv

    mw.extract_region('dat/npy', 1003, 0)
    mw.extract_region('dat/npy', 1035, 1)
    mw.extract_region('dat/npy', 2003, 2)
    mw.extract_region('dat/npy', 2035, 3)

    mw.vlm2trn('dat/1003/vlm', 'dat/1003/trn')
    mw.vlm2trn('dat/1035/vlm', 'dat/1035/trn')
    mw.vlm2trn('dat/2003/vlm', 'dat/2003/trn')
    mw.vlm2trn('dat/2035/vlm', 'dat/2035/trn')
    pass

def main():
    pass
    # tar2csv.tar2csv('raw/vtx.tar.gz', 'dat/csv')
    
if __name__ == "__main__":
    test()
