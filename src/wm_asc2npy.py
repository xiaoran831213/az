#!/usr/bin/env python
## transform ascii formated vertex into Numpy binary, also combine two hemisphere during
## the process
## the ascii filename is supposedly {sbj_id}.{lh|rh}.asc
import pdb
import numpy as np
import os
import os.path as pt
import argparse
import sys

def main():
    arg = parser.parse_args()

    ## source table
    arg.sbj = pt.normpath(arg.sbj)
    arg.sbj = pt.expandvars(arg.sbj)
    arg.sbj = pt.expanduser(arg.sbj)

    ## destination
    arg.dst = pt.normpath(arg.dst)
    arg.dst = pt.expandvars(arg.dst)
    arg.dst = pt.expanduser(arg.dst)
    if pt.isdir(arg.dst):
        npz = pt.basename(arg.sbj) + '.npz'
        arg.dst = pt.join(arg.dst, npz)

    if pt.exists(arg.dst):
        print arg.dst + ": exists"
        return

    lh = np.genfromtxt(arg.sbj + '.lh.asc', delimiter = arg.dlm, names = arg.hdr)
    rh = np.genfromtxt(arg.sbj + '.rh.asc', delimiter = arg.dlm, names = arg.hdr)
    np.savez_compressed(arg.dst, lh=lh, rh=rh)
    print arg.dst + ": created"
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = """ Transform ascii formated vertex into Numpy binary,
        also combine two hemisphere during the process.
        The ascii data file name is supposedly "{sbj_id}.{lh|rh}.asc".
        """)
    
    parser.add_argument(
        'sbj',
        help = """ The subject root name. """)
    parser.add_argument(
        '--dst', '-d', default = '.',
        help = """ Destination to store the Numpy binary.
        """)
    parser.add_argument(
        '--delimiter', '--dlm', dest = 'dlm', 
        help = """
        the column delimiter. the default is white spaces.
        """)
    parser.add_argument(
        '--header', '--hdr', action = 'store_const', dest = 'hdr',
        const = True,
        help = """
        tells the program that the first row should be treated as column header,
        consequently the data fields should be named accordingly.
        """)

    if len(sys.argv) < 2:
        parser.print_help()
    else:
        main()    
