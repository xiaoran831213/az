source('src/img.R')

cml.img <- function()
{
    argv <- commandArgs(trailingOnly = TRUE)
    if(length(argv) < 1L)
        return(NULL)
    
    library(argparser)
    p <- arg_parser('AZ image processing.')
    p <- add_argument(
        p, 'src', help = 'source directory of surface archives.')
    p <- add_argument(
        p, 'ssn', help = 'name of the surface sample to be processed')
    p <- add_argument(
        p, '--dst', help = 'target directory to store processed images.',
        default = '.')
    p <- add_argument(
        p, '--act', help = 'action to be done with the sample.',
        default='ini.img')
    
    opt <- parse_args(p, argv)
    attach(opt)
    act <- getFunction(act)
    img <- readRDS(file.path(src, paste(ssn, 'rds', sep='.')))
    img <- act(img)
    saveRDS(img, file.path(dst, paste(ssn, 'rds', sep='.')))
    detach(opt)
}

## run the command line parser
#opt <- cml.img()
