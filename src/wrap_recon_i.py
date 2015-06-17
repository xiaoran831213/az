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

def __read_manifest__(manifest):

    fi = open(manifest, 'rb')
        
    ## prepare csv reader, deduce field name from header
    ls = []
    for line in csv.DictReader(fi,fieldnames = None):
        ls.append(tuple([t(line[k]) for k, t in __CSV]))
    mf = np.array(ls, dtype = __NPY)
    fi.close()
    return mf

def __get_image_dir__(adni, record):
    whr = pt.join(adni, record['sbj'], record['dsc'])
    whr = whr.replace(' ', '_')
    fmt = record['fmt'].lower()
    for pwd, dirs, files in os.walk(whr):
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

def get_recon_iscp(adni, manifest):

    ## resolve ADNI location
    adni = pt.expanduser(adni)
    adni = pt.expandvars(adni)
    adni = pt.abspath(adni)
    if not adni.endswith("ADNI"):
        adni = pt.join(adni, "ADNI")
    
    ## load manifest CSV and sort by subject id
    mf = __read_manifest__(manifest)
    mf = mf[mf['sbj'].argsort()]

    ## get unique subject's starting indices
    S = np.unique(mf['sbj'], return_index = True)[1]

    ## container to hold the combined recon-all command
    script = []

    for sbj in np.split(mf, S[1:]):
        cmd = ''
        for img in sbj:
            img_dir = __get_image_dir__(adni, img)
            if img_dir == None:
                continue
            cmd += ' -i ' + gg(pt.join(img_dir, '*'))[0]
        if len(cmd) < 1:
            continue
        cmd = 'recon-all' + cmd + ' -s ' + img['sbj']
        cmd = cmd + '|| echo -n'
        script.append(cmd)
    return script

def main():
    import os
    import sys
    if len(sys.argv) < 3:
        print "Usage: ", sys.argv[0], " <ADNI> <manifest> <*dst>" 
        return
    
    src = sys.argv[1]
    lst = sys.argv[2]
        
    dst = sys.stdout
    if len(sys.argv) > 2:
        dst = file.open(sys.argv[2], 'wb')

    scp = get_recon_iscp(src, lst)
    dst.write(scp)
        
if __name__ == "__main__":
    pass
