#!/bin/bash

dir="$SUBJECTS_DIR"
avg="$FREESURFER_HOME/subjects/fsaverage"
sbj=
dst=
log=

## helper
function help()
{
    echo "Usage: $0 -s <SBJ> -d <DST> -hpf"
    echo "Options:"
    printf "%4s  %s\n" "-s" "subject id"
    printf "%4s  %s\n" "-d" "destination directory"
    printf "%4s  %s\n" "-h" "show this help"
    printf "%4s  %s\n" "-p" "do a dry run, which only print qsub commands."
    printf "%4s  %s\n" "-f" "force conversion even if the output exists."
}

## get parameters
function opts()
{
    while getopts ":hpfs:d:" opt; do
	case $opt in
            h)help;exit 0;;
            p)dry=1;;
            f)frc=1;;
            s)sbj="$OPTARG";;
            d)dst="$OPTARG";;
            \?)
		echo "Invalid option: -$OPTARG" >&2
		help;exit 1;;
	esac
    done
}

## argument check and report
function args()
{
    if [ ! -d "$dir" ]; then
	echo "the subject directory '$dir' is invalid." >&2
	exit 1
    fi

    if [ ! -e "$avg" ]; then
	echo "the average subject '$avg' does not exists" >&2
	exit 1
    fi

    ## check subject list
    if [ -z $sbj ]; then
        echo "must specify subject id" >&2
	help
	exit 1
    fi

    ## check destination directory
    if [ -z "$dst" ]; then
        dst="."                 # use default destination
    else
        dst=${dst%/};           # remove trailing slash
    fi
    if [ -f $dst ]; then
        echo "$dst is a file, not a directory." >&2
	help
        exit 1
    fi
    if [ ! -e $dst ] && ! mkdir "$dst"; then
        echo "failed to create directory $dst." >&2
	help
        exit 1
    fi
    log="$dst/$sbj.log"
}

## get surface vertex values aligned to the average brain
function main()
{
    hms=(lh rh)
    vls=(area curv sulc thickness)
    echo -n "" > "$log"
    for hm in ${hms[@]}; do
	fo="$dst/$sbj.$hm.avg"
	if [ -e $fo ]; then
	    echo "xt: $fo exists." | tee -a "$log"
	    continue
	fi

	## extract area, curvature, convexity and thickness
	for t in ${vls[@]}; do
	    ## extract surface values and paint them to the average brain
	    ## make sure the symbolic link to the fsaverage in freesurfor home is in the
	    ## subject home directory
	    pt="/tmp/$hm.$sbj.$t"
	    fo=$pt
	    if [ -e $fo ]; then
		echo "xt: $fo exists." | tee -a "$log"
	    elif mri_surf2surf --hemi $hm --srcsubject $sbj --sval $t --src_type curv \
		--trgsubject fsaverage  --tval $fo --trg_type curv >> "$log"; then
 		echo "xt: $fo created." | tee -a "$log"
	    else
		exit 1
	    fi

	    ## write the painted values to ASCII file
	    fi="$fo"
	    fo="$pt.asc"
	    if [ -e $fo ]; then
		echo "xt: $fo exists." | tee -a "$log"
	    elif mris_convert -c $fi $avg/surf/$hm.white $fo; then
		echo "xt: $fo created." | tee -a "$log"
	    else
		exit 1
	    fi

	    ## cut the value column(the 5th) from ascii files
	    fi="$fo"
	    fo="$pt.val"
	    if [ -e $fo ]; then
		echo "xt: $fo exists." | tee -a "$log"
	    else
		cut "$fi" -d' ' -f5 > "$fo"
 		echo "xt: $fo created." | tee -a "$log"
	    fi
	done

	## extract talirah xyz
	pt="/tmp/$hm.$sbj.txyz"
	fo="$pt"
	if [ -e $fo ]; then
	    echo "xt: $fo exists." | tee -a "$log"
	elif mri_surf2surf --hemi $hm --srcsubject $sbj --sval-tal-xyz white \
	    --trgsubject fsaverage  --tval-xyz --tval $fo >> "$log"; then
	    echo "xt: $fo created." | tee -a "$log"
	else
	    exit 1
	fi

	## write talirah xyz to ascii
	fi="$fo"
	fo="$pt.asc"
	if [ -e $fo ]; then
	    echo "xt: $fo exists." | tee -a "$log"
	elif mris_convert "$fi" "$fo"; then
	    echo "xt: $fo created." | tee -a "$log"
	else
	    exit 1
	fi

	## find out # vertices and triangles
	n=($(sed $fo -n -e '2,2p'))
	n_vtx=${n[0]}
	n_tri=${n[1]}

	## extract vertex coordinate
	fi="$fo"
	fo="$pt.val"
	if [ -e $fo ]; then
	    echo "xt: $fo exists." | tee -a "$log"
	elif sed "$fi" -ne "3,$((2+n_vtx)) s/^\([^ ]*\) *\([^ ]*\) *\([^ ]*\) *0$/[\1,\2,\3]/p" \
	    >"$fo"; then
	    echo "xt: $fo created." | tee -a "$log"
	else
	    exit 1
	fi

	## concatinate every thing
	pt="/tmp/$hm.$sbj"
	fo="$pt.pst"
	if [ -e $fo ]; then
	    echo "xt: $fo exists." | tee -a "$log"
	elif paste $pt.{'txyz','area','curv','sulc','thickness'}.val > "$fo"; then
	    echo "xt: $fo created." | tee -a "$log"
	else
	    exit 1
	fi

	## final output
	fi="$fo"
	fo="$dst/$sbj.$hm.avg"
	echo -e "xyz\tare\tcrv\tsul\tthk" > "$fo"
	cat "$fi" >> "$fo"
	
	## cleanup
	echo "xt: clean up" | tee -a "$log"
	rm /tmp/$hm.$sbj*
    done
    echo "xt: success"
}

opts $@
args
main
