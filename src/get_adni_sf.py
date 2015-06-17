import sys
import csv
import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import xml.etree.ElementTree as ET

__NPY = np.dtype([
    ("img", "<i4"),
    ("sbj", "a16"),
    ("grp", "a16"),
    ("sex", "a08"),
    ("age", "<i1"),
    ("vst", "<i1"),
    ("mod", "a16"),
    ("dsc", "a32"),
    ("typ", "a16"),
    ("acq", "a16"),
    ("fmt", "a08")])
    
__CSV = [
    ("Image Data ID", int),
    ("Subject", str),
    ("Group", str),
    ("Sex", str),
    ("Age", int),
    ("Visit", int),
    ("Modality", str),
    ("Description", str),
    ("Type", str),
    ("Acq Date", str),
    ("Format", str)]

def get_manifest(manifest = None):
    """ get columns from text table """
    if manifest:
        fi = open(manifest, 'rb')
    else:
        fi = sys.stdin
        
    ## prepare csv reader and writers
    ## deduce field name from header
    ls = []
    for line in csv.DictReader(fi,fieldnames = None):
        ls.append(tuple([t(line[k]) for k, t in __CSV]))
    ar = np.array(ls, dtype = __NPY)
        
    if manifest:
        fi.close()
    return ar

def get_image_dir(adni, record):
    adni = pt.join(adni, 'ADNI', record['sbj'])
    adni = pt.join(adni, record['dsc']).replace(' ', '_')
    fmt = record['fmt'].lower()
    for pwd, dirs, files in os.walk(adni):
        if len(files) < 1:      # not file
            continue
        
        for f in files:
            ## ignore non-image and non-slice (usually some XML)
            if not f.endswith(fmt):
                continue
            ## find wrong image ID in image/slice file name
            if f.find('I' + str(record['img'])) <0:
                break

            ## find right image ID in image/slice file name
            return pwd
    else:
        return None

def get_recon_import_dcm_cmd(adni, manifest):
    ## sort by subject id
    mf = manifest[manifest['sbj'].argsort()]

    ## get unique subject's starting indices
    S = np.unique(mf['sbj'], return_index = True)[1]

    ## container to hold the combined recon-all command
    C = []

    for sbj in np.split(mf, S[1:]):
        cmd = ''
        for img in sbj:
            img_dir = get_image_dir(adni, img)
            if img_dir == None:
                continue
            cmd += ' -i ' + gg(pt.join(img_dir, '*'))[0]
        if len(cmd) < 1:
            continue
        cmd = 'recon-all' + cmd + ' -s ' + img['sbj']
        C.append(cmd)
    return C

def main():
    import os
    import sys
    if len(sys.argv) < 2:
        print "Usage: ", sys.argv[0], " <adni_csv> <dst>" 
        return
    
    src = None
    if len(sys.argv) > 1:
        src = sys.argv[1]
        
    dst = None
    if len(sys.argv) > 2:
        dst = sys.argv[2]

    get(src, dst)
        
if __name__ == "__main__":
    main()
