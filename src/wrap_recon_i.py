import sys
import csv
import pdb
import numpy as np
import os
import os.path as pt
from glob import glob as gg
import sys
import hlp

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
    """ read ADNI downloaded image manifest CSV file"""
    fi = open(manifest, 'rb')
        
    ## prepare csv reader, deduce field name from header
    ls = []
    for line in csv.DictReader(fi,fieldnames = None):
        ls.append(tuple([t(line[k]) for k, t in __CSV]))
    mf = np.array(ls, dtype = __NPY)
    fi.close()
    return mf

def __get_image_dir__(adni, record):
    """
    given subject record from the download manifest, locate the
    directory containing mri slices of the subject """
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

def write_recon_script(adni, manifest, dst = "hpc", ncpu = 8, batch_size = 64):
    hlp.mk_dir(dst)

    ## resolve ADNI and manifest location
    adni = hlp.resolve_path(adni)
    adni = pt.abspath(adni)
    adni = pt.join(adni, "ADNI")
    
    ## load manifest CSV and sort by subject id
    manifest = hlp.resolve_path(manifest)
    mf = __read_manifest__(manifest)
    mf = mf[mf['sbj'].argsort()]

    ## get unique subject's starting indices
    U, S = np.unique(mf['sbj'], return_index = True)

    ## container to hold the combined recon-all command
    i_bat = 0
    i_cmd = 0
    cmd = 'recon-all {0} -s {1} 1> {1}.out 2> {1}.err || echo -n\n'
    bat = '{}/{:03d}.qs'
    mem = ncpu * 1.0
    wtm = batch_size / ncpu * 0.2
    for scans in np.split(mf, S[1:]):
        ## find mri images in each scan group
        imgs = ''
        for img in scans:
            dr = __get_image_dir__(adni, img)
            if dr == None:
                continue
            imgs += ' -i ' + gg(pt.join(dr, '*'))[0]
        if not imgs:
            continue

        ## new batch
        if i_cmd % batch_size == 0:
            f = open(bat.format(dst, i_bat), 'wb')
            hlp.write_hpcc_header(f, mem = mem, walltime = wtm, nodes = ncpu)
            f.write('\n')
            i_cpu = 0
            
        ## new cpu
        if i_cmd % ncpu == 0:
            f.write('## node {:02d}\n'.format(i_cpu))
            f.write('(\n')
            
        ## write new command
        f.write(cmd.format(imgs, img['sbj']))
        i_cmd += 1

        ## end of one cpu
        if i_cmd % ncpu == 0:
            f.write(')&\n\n')
            i_cpu += 1
        
        ## end of one batch
        if i_cmd % batch_size == 0:
            f.write('wait\n')
            f.close()
            i_bat += 1

    ## write submitor
    if f:
        f.close()
    f = open(pt.join(dst, 'tsk.sh'), 'wb')
    for i in xrange(i_bat):
        f_bat = pt.abspath(bat.format(dst, i))
        f.write('qsub {}\n'.format(f_bat))
        

##g et_recon_script('~', '../raw/ADNI_WGS_VS1.csv', '../import.qs')

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
