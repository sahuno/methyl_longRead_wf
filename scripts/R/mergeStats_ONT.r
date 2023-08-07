#samuel ahuno
#date - july 11th 2023
#purpose - merge all aggregate stats datasets

#load libraries
library(optparse)
library(data.table)

#capture arguments
option_list <- list(make_option(c("-i","--inputFiles"), default=c("/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/per_read_aggregate/Spectrum-OV-009_N/chr8/data/Spectrum-OV-009_N.chr8.data_aggregate_stats.txt",
                                                            "/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/spectrum_ont_methyl/results/per_read_aggregate/Spectrum-OV-009_N/chr18/data/Spectrum-OV-009_N.chr18.data_aggregate_stats.txt"),
                                                    action="store_true", dest="verbose"),
                    make_option(c("--outputFiles"), default=NULL,action="store_true", dest="verbose"))


parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options
args <- arguments$args
##print(opt)
#print(opt)

#print(str(opt))
#print(str(args))
#print(names(args))


myfilelist <- strsplit(opt$inputFiles, ",")

print(myfilelist)
# print("\nargs")
# print (args)


# opt <- parse_args(OptionParser(option_list = option_list))
# print(opt)
# message(opt$inputFiles)

#read files and merge
# l <- lapply(opt$inputFiles, fread, sep="\t")
# head(l)
# dt <- rbindlist(l)
# dt <- setnames(dt, "chrm", "seqnames")
# setkey(dt , seqnames, pos) #setkey for fast sorting
# dt <- dt[order(seqnames, pos)]


# #save on disk
# fwrite(dt, file = opt$outputFiles, sep="\t")

