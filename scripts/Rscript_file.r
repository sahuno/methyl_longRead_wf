###!bin

#template Rscript
#Samuel Ahuno
#date; May 25nd 2023

#script takes, 2 arguments
#path to sample file,  output file

#load libraries
library(optparse)
library(tidyverse)
library(magrittr)
library(BiocParallel)#used for parallel analysis

#print out date
d1<-date()
message(paste0("date time started - "))
message(d1)


########################################################
######accept arguments from environment
########################################################
option_list_val <- list(make_option(c("-i","--input_file"), type="character", default=NULL, help = "methylation coverage file"),
                        make_option(c("-o","--output_file"), type="character", default=NULL, help = "file name for processed file. Extention is .rds"))

parseObject <- OptionParser(option_list = option_list_val)
opt <- parse_args(parseObject)
print(opt)


########################################################
#sanity check
if (is.null(opt$sample)){
  print_help(opt_parser)
  stop("Please supply sample to work on", call.=FALSE)
}

####################################################################
message("reading files")
#read files from directory
df_input <- fread(bs_unsmoothed, file = opt$output_file)

#process file
proccessed_file <- summary(df_input)

#save file to disk
saveRDS(proccessed_file, file = opt$output_file)