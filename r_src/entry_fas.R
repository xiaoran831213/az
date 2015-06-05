source('src/gno.R');
source('src/utl.R');
source('src/gfa.R');

## wnd --- segment window size
## cnd --- number of candidate segments to pick
## vcf --- VCF file to read. (Varient Call Format)
## idv --- IDV file to read. (Individual)
## seg --- SEG file to read. (Segment file, e.g. Genes)
go <- function(vcf, idv='dat/idv.EUR', seg='dat/gen', wnd=5000L, cnd=1L)
{
    ## the return list
    ret<-list();
    
    ## load segment list & genotype data.
    seg<-GNO$seg(whr = seg, chr = 3L);
    
    ## load individual list, they are shared among all candicate segments
    idv<-read.table(file=idv, header = F, as.is = T);
    
    ## pick candidate segments
    cnd<-sample.int(n = nrow(seg), size = cnd, replace = F)
    cnd<-sort(cnd);
    
    ## fit function for genotype accross candicate segments
    for(i in cnd)
    {
        ## the range
        r<-seg[i,, drop=T];
        
        ## extract genotype
        g<-GNO$vcf(vcf = vcf, chr = r$chr, bp1 = r$bp1 - wnd, bp2 = r$bp2 + wnd, idv = idv);
        if(!is.null(g$err))
        {
            cat(g$err, '\n', file = stderr());
            next;
        }
        g<-GNO$imp(g); ## impute missing genotype values
        
        ndv<-nrow(g$idv);
        ngv<-nrow(g$map);
        gmx<-g$gmx;

        key<-sprintf('G%s.%s', r$seq, r$gen);
        ret[[key]]<-g;
    }
    ret;
}
