#!/bin/bash

# where am I, and my nearby resources?
if [ -h $0 ]; then
    app=$(readlink $0);
else
    app=$0
fi
# mapping from vertex SN to annotation SN
v2a=${app%/*}/v2a.txt	     

# freesurfer subject directory and  altalas 'fsaverage'
dir="$SUBJECTS_DIR"
avg="$FREESURFER_HOME/subjects/fsaverage"

# script options
sid=				# subject id
dst=.				# target file
skp=0				# skip existing
log=

## helper
function help()
{
    echo "Usage: $0 -i <SID> -o <OUT> -hpk"
    echo "Options:"
    echo "  -i subject id"

    echo "  -o output path. if the path leads to an existing directory,"
    echo "     the output is '<OUT>/<SID>.txt."
    
    echo "  -h show this help."
    echo "  -p do a dry run, only print out commands."
    echo "  -k skip existing output."
    echo "  -l log file, def=STDOUT."
}

## get parameters
function opts()
{
    while getopts ":hpki:o:l:" opt; do
	case $opt in
            h)help;exit 0;;
            p)dry=1;;
            k)skp=1;;
            i)sid="$OPTARG";;
            o)dst="$OPTARG";;
	    l)log="$OPTARG";;
	    z)gzp=1;;
            \?)
		echo "Invalid option: -$OPTARG" >&2
		help;exit 1;;
	esac
    done
}

## argument check and report
function args()
{
    ## check logging
    exec 6>&1		# save stdout
    exec 7>&2		# save stderr
    if [ ! -z "$log" ]; then
	exec 1>"$log"
	exec 2>"$log"
	echo "LOG=$log"
    fi

    echo "APP=$app"
    
    ## check the vertex to annotation map
    if [ ! -e "$v2a" ]; then
	echo "the vertex to annotation mapping is invalid."
	exit 1
    fi
    echo "V2A=$v2a"
    ## check FreeSurfer subject directory
    if [ ! -d "$dir" ]; then
	echo "the subject directory '$dir' is invalid, please check" >&2
	echo "the FreeSurfer installation and configuration." >&2
	exit 1
    fi
    echo "DIR=$dir"

    ## check FreeSurfer average brain
    if [ ! -e "$avg" ]; then
	echo "the average subject '$avg' is invalid, please check" >&2
	echo "the FreeSurfer installation and configuration." >&2
	exit 1
    fi
    echo "AVG=$avg"

    ## check subject id
    if [ -z "$sid" ]; then
        echo "must specify subject id" >&2
	help
	exit 1
    fi
    echo "SID=$sid"
    if [ ! -z "$log" ]; then	# SID is always screen printed
       echo "SID=$sid" >&6
    fi

    ## check destination
    if [ -z "$dst" ]; then
        dst="."                 # use default destination
    fi
    if [ -d "$dst" ]; then	# dist is a directory
       dst="${dst%/}/$sid.txt"	# make up a full path
    fi
    echo "DST=$dst"

    ## check skipping
    if [ $skp -ne 0 ]; then
       echo "SKP=$skp"
       if [ -e "$dst" ]; then
	   echo "$dst exists, skip." >&2
	   exit 0
       fi
    fi
}

# get surface vertex values aligned to the average brain
function main()
{
    # create output file, write table header
    hdr=""
    hdr+="hemi\t"		# hemisphare (0=LH, 1=RH)
    hdr+="vseq\t"		# vertex SN (1~163842)
    hdr+="aseq\t"		# annotation SN (0~36)
    hdr+="x\ty\tz\t"		# vertex coordinates
    hdr+="area\t"		# area around vertex
    hdr+="curv\t"		# curvature
    hdr+="sulc\t"		# sulcity
    hdr+="thck"			# WM thickness

    hms=(lh rh)
    vls=(txyz area curv sulc thickness)

    # extract two hemispheres
    for hm in ${hms[@]}; do

	# extract area, curvature, convexity and thickness
	for t in ${vls[@]:1}; do
	    # extract surface values and paint them to the average brain
	    # make sure the symbolic link to "fsaverage" is present in $SUBJECTS_HOME
	    pt="/tmp/$hm.$sid.$t"
	    fo=$pt
	    if [ -e $fo ]; then
		echo "xt: $fo exists."
	    elif mri_surf2surf --hemi $hm --srcsubject $sid --sval $t --src_type curv \
			       --trgsubject fsaverage --tval $fo --trg_type curv; then
 		echo "xt: $fo created."
	    else
		exit 1
	    fi

	    # write the painted values to ASCII file. It doesn't matter which surface
	    # is being painted on (e.g., white, inflated, sphere, sphere.reg), as long
	    # as the surface is of fsaverage, since we only care about the mapping from
	    # vertices of native subject to those in the altlas 
	    fi="$fo"
	    fo="$pt.asc"
	    if [ -e $fo ]; then
		echo "xt: $fo exists."
	    elif mris_convert -c $fi $avg/surf/$hm.white $fo; then
		echo "xt: $fo created." 
	    else
		exit 1
	    fi

	    # cut the value column(the 5th) from ascii files
	    fi="$fo"
	    fo="$pt.val"
	    if [ -e $fo ]; then
		echo "xt: $fo exists."
	    else
		cut "$fi" -d' ' -f5 > "$fo"
 		echo "xt: $fo created."
	    fi
	done

	# extract talirah xyz
	pt="/tmp/$hm.$sid.txyz"
	fo="$pt"
	if [ -e $fo ]; then
	    echo "xt: $fo exists."
	elif mri_surf2surf --hemi $hm --srcsubject $sid --sval-tal-xyz white \
			   --trgsubject fsaverage --tval-xyz --tval $fo; then
	    echo "xt: $fo created."
	else
	    exit 1
	fi

	# write talirah xyz to ascii
	fi="$fo"
	fo="$pt.asc"
	if [ -e $fo ]; then
	    echo "xt: $fo exists."
	elif mris_convert "$fi" "$fo"; then
	    echo "xt: $fo created."
	else
	    exit 1
	fi

	# find out # vertices and triangles
	n=($(sed $fo -ne '2,2p; 2q'))
	n_vtx=${n[0]}
	n_tri=${n[1]}

	# extract vertex coordinate
	fi="$fo"
	fo="$pt.val"
	if [ -e $fo ]; then
	    echo "xt: $fo exists."
	elif sed "$fi" -n \
		 -e "3,$ s/^\([^ ]*\) *\([^ ]*\) *\([^ ]*\) *0$/\1\t\2\t\3/p" \
		 -e "$((2+n_vtx))q" >"$fo"; then
	    echo "xt: $fo created."
	else
	    exit 1
	fi

	# paste vertex data:
	pt="/tmp/$hm.$sid"
	fo=$pt.val
	if [ -e $fo ]; then
	    echo "xt: $fo exists."
	elif paste $pt.{'txyz','area','curv','sulc','thickness'}.val > $fo; then
	    echo "xt: $fo created."
	else
	    exit 1
	fi
    done

    # concatinate vertext data of both hemisphere (LH first, then RH), then
    # pre-append the vertex to annotation mapping
    if cat <(echo -e "$hdr") \
	   <(paste <(tail -n+2 $v2a) \
		   <(cat /tmp/{lh,rh}.$sid.val)) \
	   > "$dst"; then
	echo "xt: $dst created."
    else
	exit 1
    fi
	
    # cleanup
    echo "xt: clean up"
    rm -rf /tmp/{lh,rh}.$sid.*
    echo "xt: success"
}

opts $@
args
main
