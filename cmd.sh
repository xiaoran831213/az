#!/bin/bash

dir="$SUBJECTS_DIR"
avg="$FREESURFER_HOME/subjects/fsaverage"

## extract surface vertices with respect to the the average brain
get_avg_vtx()
{
    sbj=$1
    for hm in {'lh','rh'}; do
	## extract area, curvature, convexity and thickness
	for t in {'area','curv','sulc','thickness'}; do
	    ## extract surface values and paint them to the average brain
	    fo="/tmp/$hm.$sbj.$t"
	    if [ -e $fo ]; then
		echo "xt: $fo exists."
	    else
 		echo "xt: paint $t to $fo"
		mri_surf2surf --hemi $hm --srcsubject $sbj --sval $t --src_type curv \
		    --trgsubject fsaverage  --tval $fo --trg_type curv
	    fi

	    ## write the painted values to ASCII file
	    if [ -e "$fo.asc" ]; then
		echo "xt: $fo.asc exists."
	    else
 		echo "xt: write $fo.asc"
		mris_convert -c $fo $avg/surf/$hm.white "$fo.asc"
	    fi
	done

	## 
	fo="/tmp/$hm.$sbj.xyz"
	if [ -e $fo ]; then
	    echo "xt: $fo exists."
	else
	    echo "xt: $fo <-- paint xyz to average $hm"
	    mri_surf2surf --hemi $hm --srcsubject $sbj --sval-xyz white \
		--trgsubject fsaverage  --tval-xyz --tval $fo
	fi

	if [ -e $fo".asc" ]; then
	    echo "xt: $fo.asc exists."
	else
	    echo "xt: $fo.asc <-- write $fo to ascii"
	    mris_convert $fo "$fo.asc"
	fi
    done

    n=$(sed $fo.asc -n -e '2,2p')
    n=($n)
    n_vtx=${n[0]}
    n_tri=${n[1]}
    
}

# for i in ${!info[@]}
# {
#     echo $i ${info[i]}
# }
