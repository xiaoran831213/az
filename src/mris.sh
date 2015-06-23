#!/bin/bash

dir="$SUBJECTS_DIR"
avg="$FREESURFER_HOME/subjects/fsaverage"

## extract surface vertices with respect to the the average brain
get_vtx_by_ata()
{
    ## find destination
    sbj=$1
    if [ -z "$2" ]; then
	dst="."
    else
	dst="$2"
    fi
    if [ -d "$dst" ]; then
	log="$dst/$sbj.$hm.log"
    else
	echo "xt: $dst is not a directory" | tee "$log"
	exit 1
    fi

    ## work on both hemesphere
    for hm in {'lh','rh'}; do
	fo="$dst/$sbj.$hm.avg"
	if [ -e $fo ]; then
	    echo "xt: $fo exists." | tee -a "$log"
	    continue
	fi

	## extract area, curvature, convexity and thickness
	for t in {'area','curv','sulc','thickness'}; do
	    ## extract surface values and paint them to the average brain
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
	n=$(sed $fo -n -e '2,2p')
	n=($n)
	n_vtx=${n[0]}
	n_tri=${n[1]}

	## extract vertex coordinate
	fi="$fo"
	fo="$pt.val"
	if [ -e $fo ]; then
	    echo "xt: $fo exists." | tee -a "$log"
	elif sed "$fi" -ne \
	    "3,$((2+n_vtx)) s/^\([^ ]*\) *\([^ ]*\) *\([^ ]*\) *0$/[\1,\2,\3]/p" \
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
	echo "clean up" | tee -a "$log"
	rm /tmp/$hm.$sbj*
    done
}

