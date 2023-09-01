#test results

library(data.table)
dt <- fread("/juno/work/greenbaum/users/ahunos/methyl_SPECTRUM/scripts/workflows/methyl_PARP_BrCan/results/addHeader/BRCA_13135_P_2/data/BRCA_13135_P_2.chr1.per_read_modified_base_calls.header.txt")

#how many unique reads are there?
dt[mod_base == "m",.N,by=.(read_id)]
#dt[,.N,by=.(read_id)]
dt_m_sort <- dt[mod_base == "m",][order(c(chrm,pos,strand))][!is.na(read_id),][,`:=`(exp_mod_log_prob = exp(mod_log_prob))]
#dput(head(dt_m_sort[,.(read_id,chrm,pos,strand,exp_mod_log_prob)], 100))
# dt[mod_base == "m" & read_id=="9c15ec94-4529-4675-93b3-1e8acb034c2f",]

dt2 <- head(dt_m_sort[,.(read_id,chrm,pos,strand,exp_mod_log_prob)], 100)
# Filter rows where exp_mod_log_prob > 0.9
dt2 <- dt2[exp_mod_log_prob > 0.9]

# Create a column representing consecutive groups of rows with the same read_id
dt2[, grp := rleid(read_id)]

# Calculate the mean of exp_mod_log_prob for each group of 4 consecutive rows
result <- dt2[, .(mean_exp_mod_log_prob = mean(exp_mod_log_prob)), by = .(chrm, grp)]#[, .(mean_mean_exp_mod_log_prob = mean(mean_exp_mod_log_prob))]

# Print the result
print(result)
# sum(lag(1:5))