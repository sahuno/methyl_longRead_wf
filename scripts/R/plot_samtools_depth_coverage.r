#### samuel ahuno
####plot depth of coverage at each cpg site

library(data.table)
library(optparse)
library(tidyverse)

option_list <- list(make_option("--depth_file", type="character", default=NULL, help="depth file from samtools depth"),
make_option("--paths_bams_txt", type="character", default=NULL, help="list of bam file (.txt) used to generate depth file"))
opt <- parse_args(OptionParser(option_list = option_list))

#load data
dt <- fread(opt$depth_file, header = TRUE, sep = "\t")
dt_p <- fread(opt$paths_bams_txt, header = F)
dt_p[, samples := gsub(".*/megalodon/|/mappings.sorted.bam", "", V1)]
names(dt) <- c("chr","pos" ,dt_p$samples)


total_cpg_by_chr <- dt[,.N,by=chr] #get expected number of CpGs by chr
mapped_reads_cpg_sites_dt <- dt[, lapply(.SD, function(x) sum(x>=1)), .SDcols = !c("chr","pos"), by = chr][total_cpg_by_chr, on = "chr"]

mapped_reads_cpg_sites_df <- mapped_reads_cpg_sites_dt %>% dplyr::rowwise() %>% 
        mutate(across(!starts_with("N"), ~ (.x/N))) %>% 
                dplyr::select(-N) %>% 
                    pivot_longer(!chr, names_to = "sample", values_to = "percent_mapped_reads_not_0x") 


mapped_reads_cpg_sites_ls <- split(mapped_reads_cpg_sites_df, mapped_reads_cpg_sites_df$chr)
function_plot <- function(df,chr){
        plt_per_mapped  <- ggplot(df , aes(x=sample, y=percent_mapped_reads_not_0x)) + geom_col() + 
        theme(axis.ticks.x=element_blank(), 
              axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1)) + 
        scale_y_continuous(labels = scales::percent) + 
        labs(title = paste0("percentage of number of CpGs with >=1 mapped reads  - ",chr).
        y = "100 * (n_CpGs >=1 mapped reads)/total number of CpGs", x = "sample")
}


mapped_reads_cpg_sites_plots_ls <- imap(mapped_reads_cpg_sites_ls, ~function_plot(df=.x, chr=.y))
pdf(file="percent_mapped_reads_not_0x_epi_triplicate_mm10.pdf", width=10, height=10)
print(mapped_reads_cpg_sites_plots_ls)
dev.off()



# #plot histogram of mapped reads
# plt_hist <- ggplot(dt, aes(`D-Q-3`))+ geom_histogram() + 
# geom_histogram(alpha = 0.5, position="identity") +
# scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                   labels = scales::trans_format("log10", scales::math_format(expr = 10^.x)))
# ggsave(plt_hist, file="histogram_mapped_reads_epi_triplicate_mm10.pdf", width=10, height=10)
