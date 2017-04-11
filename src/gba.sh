## seperate genes from raw VCF files.
## make sure the symbolic link to the decoded dbGap files is created in 'raw'
wd=$(pwd)			# remember root directory

# genes separated from raw VCF files, version #1
sg=dat/gs1

# genomic analysis
dd=out/gs2/1_ent

for f in $sg/*.rds
do
    g=${f##*/}
    if [ -e $dd/$g ]; then
	continue
    else
	echo "time ./gba.R $f $g"
    fi
done | hpcwp - -d $dd --wtm 4 -m8 -n1 -q30 --ln dat --ln src --ln raw --cp src/gba.R --md GNU/4.9 --md R/3.3.0

## find out failed qsub command in a parallel computation directory
for s in $(comm <(ls -1 log) <(cat pbs/* | sed -n 's/^.*log\/\(.*\)$/\1/p') -3 | cut -f2 | xargs -I{} grep -l {} pbs/* | sort | uniq)
do
    echo $s
done

## find out jobs hitting walltime limit and resubmit them
grep walltime std/* | sed 's|^std/\([0-9][0-9][0-9].sh\).*$|pbs/\1|' | xargs -I{} qsub {}
