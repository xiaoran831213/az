#!/bin/bash

dir="$SUBJECTS_DIR"
avg="$FREESURFER_HOME/subjects/fsaverage"

## extract surface vertices with respect to the the average brain
align_vertex()
{
    sbj=$1
    for hm in {'lh','rh'}; do
	## extract area, curvature, convexity and thickness
	for t in {'area','curv','sulc','thickness'}; do
	    ## extract surface values and paint them to the average brain
	    pt="/tmp/$hm.$sbj.$t"
	    fo=$pt
	    if [ -e $fo ]; then
		echo "xt: $fo exists."
	    else
 		echo "xt: paint $t to $fo"
		mri_surf2surf --hemi $hm --srcsubject $sbj --sval $t --src_type curv \
		    --trgsubject fsaverage  --tval $fo --trg_type curv
	    fi

	    ## write the painted values to ASCII file
	    fi="$fo"
	    fo="$pt.asc"
	    if [ -e $fo ]; then
		echo "xt: $fo exists."
	    elif [ mris_convert -c $fi $avg/surf/$hm.white $fo ]; then
		echo "xt: $fo created."
	    else
		echo "xt: $fo failed."
	    fi

	    ## cut the value column(the 5th) from ascii files
	    fi="$fo"
	    fo="$pt.val"
	    if [ -e $fo ]; then
		echo "xt: $fo exists."
	    else
		cut "$fi" -d' ' -f5 > "$fo"
 		echo "xt: $fo created"
	    fi
	done

	## extract talirah xyz
	pt="/tmp/$hm.$sbj.txyz"
	fo="$pt"
	if [ -e $fo ]; then
	    echo "xt: $fo exists."
	elif [ mri_surf2surf --hemi $hm --srcsubject $sbj --sval-tal-xyz white \
		--trgsubject fsaverage  --tval-xyz --tval $fo ]; then
	    echo "xt: $fo created."
	else
	    echo "xt: $fo failed."
	    exit 1
	fi

	## write talirah xyz to ascii
	fi="$fo"
	fo="$pt.asc"
	if [ -e $fo ]; then
	    echo "xt: $fo exists."
	elif [ mris_convert "$fi" "$fo" ]; then
	    echo "xt: $fo created."
	else
	    echo "xt: $fo failed."
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
	    echo "xt: $fo exists."
	else
	    sed "$fi" -ne "3,$((2+n_vtx)) s/^\([^ ]*\) *\([^ ]*\) *\([^ ]*\) *0$/[\1,\2,\3]/p" >"$fo"
	    echo "xt: $fo created."
	fi

	## concatinate every thing
	pt="/tmp/$hm.$sbj"
	fo="$pt.pst"
	if [ -e $fo ]; then
	    echo "xt: $fo exists."
	else
	    paste $pt.{'txyz','area','curv','sulc','thickness'}.val > "$fo"
	    echo "xt: $fo created."
	fi

	## write to destination
	if [ -z "$2" ]; then
	    dst="."
	else
	    dst="$2"
	fi
	if [ -d "$dst" ]; then
	    dst="$dst/$sbj.$hm.avg"
	else
	    echo "xt: $dst is not a directory"
	    exit 1
	fi

	fi=$fo
	if [ -e "$dst" ]; then
	    echo "xt: $dst exists."
	else
	    echo -e "xyz\tare\tcrv\tsul\tthk" > $dst
	    cat "$fo" >> $dst
	fi
    done
}

# for i in ${!info[@]}
# {
#     echo $i ${info[i]}
# }
