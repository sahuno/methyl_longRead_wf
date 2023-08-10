#make sample list files
library(tidyverse)
library(data.table)
library(yaml)


path_megalodon_dir <- "/juno/work/greenbaum/projects/TRI_EPIGENETIC/megalodon"
files_mod_bases <- list.files(path_megalodon_dir, full.names= TRUE, pattern="er_read_modified_base_calls.txt", recursive = TRUE) %>% as_tibble()
files_mod_bases

files_mod_bases2yaml <- files_mod_bases %>% mutate(samples = gsub("/per_read_modified_base_calls.txt","",gsub(".*megalodon/","",value))) %>% as.data.table()
#names(files_mod_bases2yaml) <- NULL

# cat(yaml::as.yaml(list(samples = split(replace(files_mod_bases2yaml, "samples", NULL), f = files_mod_bases2yaml$samples))))
# cat(yaml::as.yaml(list(samples = split(files_mod_bases2yaml, f = files_mod_bases2yaml$samples))))
# split(files_mod_bases2yaml, f = files_mod_bases2yaml$samples)


files_mod_bases2yaml <- files_mod_bases2yaml %>% mutate(samples= paste0(samples, ": ", value )) %>% dplyr::select(-value) %>% as.data.table()
write_yaml(files_mod_bases2yaml, file="files_mod_bases2yaml.yaml")





# files_mod_bases2yaml <- split(files_mod_bases2yaml, f = files_mod_bases2yaml$samples)
# names(files_mod_bases2yaml) <- NULL
# files_mod_bases2yaml <- files_mod_bases2yaml %>% split(by="samples") %>%
#   lapply(function(x) x[,sample := NULL] %>% .[])

# cat(yaml::as.yaml(files_mod_bases2yaml))



