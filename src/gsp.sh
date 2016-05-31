## separate genes from VCFs

tpd=/tmp/xtong
mkdir -p $tpd

# download gene list in refFlat format
gls=raw/refGen19.txt

if [ ! -e $gls ]; then
    fi=http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz
    fo=$tpd/refFlat19.txt.gz
    wget $fi -O $fo

    fi=$fo
    fo=$tpd/gl1.txt
    zcat $fi | cut -f1,3,5,6 > $fo

    fi=$fo
    fo=$tpd/gl2.txt
    sed $fi -ne 's/\tchr\([0-9XYM]\{,2\}\)\t/\t\1\t/p' > $fo

    fi=$fo
    fo=$tpd/gl3.txt
    sed $fi -e 's/\tX\t/\t23\t/' | sed 's/\tY\t/\t24\t/' | sort -k2n -k3n -k4n > $fo

    fi=$fo
    fo=$tpd/gl4.txt
    uniq $fi > $fo

    # move gene list to the project direcotry, remove temporary files
    mv $fo $gls
fi

# extract uniqe genes
fi=$gls
fo=$tpd/gs0.txt
smb_last=0
chr_last=0
bp0_last=0
bp1_last=0
s=0
while read smb chr bp0 bp1; do
    # is the previous gene overlapped to much with the current one?
    if [ $chr_last -eq $chr -a $((bp1_last-bp0)) -gt $((bp1-bp1_last)) ]; then
	mrg=1
    elif [ $smb_last = $smb ]; then
    	mrg=1
    else
	mrg=0
    fi

    # merge, or take new gene
    if [ $mrg -eq 1 ]; then
	smb_last=$smb
	if [ $bp1 -gt $bp1_last ]; then
	    bp1_last=$bp1
	fi
    else
	echo -e "$(printf G%04X $s)\t$chr_last\t$bp0_last\t$bp1_last\t$smb_last"
	smb_last=$smb
	chr_last=$chr
	bp0_last=$bp0
	bp1_last=$bp1
	((s+=1))
    fi
done < $fi > $fo
echo -e "$(printf G%04X $s)\t$chr_last\t$bp0_last\t$bp1_last\t$smb_last" >> $fo
tail -n+2 $fo | sed 's/\t23\t/\tX\t/; s/\t24\t/\tY\t/' > raw/gs1.txt


#rm -rf $tpd

# separate vcf files
while read seq chr bp0 bp1 smb; do
    cmd="vcf2dsg.R --rgn $chr:$bp0-$bp1 --wnd 5000 --smb \"$smb\" wgs/$chr.vcf.gz $seq.rds"
    echo "$cmd"
done < raw/gs1.txt | hpcc_wrapper - -d dat/gs1 --wtm 4 --mem 0.5 -n1 -q2896 --cp raw/wgs --md R/3.2.0
