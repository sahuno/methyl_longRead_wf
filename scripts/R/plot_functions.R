##plots

dot_plots <- function(data_in, x_var, y_var){
p <- ggplot(data_in, aes(x={{x_var}}, y={{y_var}})) + 
  geom_point(size=1, alpha=0.3) 
return(p)
}

gglzy <- function(ggobj,file_p){
  is.null(file_p){
    message("saving at current dir")
      ggsave("lzw_plot.tiff", plot = ggobj, width=300, height=225, units="mm", dpi=300, compression = "lzw")
  }else{
ggsave(filename=file_p, plot = ggobj, width=300, height=225, units="mm", dpi=300, compression = "lzw")
  }
}

dens_plt <- function(data_in, x_var, y_var, groups_var){
plt_density <-  ggplot(data_in, 
                                    aes(x={{x_var}}, color = {{groups_var}})) + 
                                  geom_density() + 
                                  theme(legend.position="bottom") #+ 
                                  # facet_wrap(~methyl_state, scales = "free_x")
                                  #labs(title = paste0("density of proportions motif reads at genomic site", ls_names_to_append[1]))
# ggsave(plt_density_stats_meth, file="plot_density_mod_prop.pdf")
return(plt_density)
}

#histogram
hist_plt <- function(data_in, x_var, y_var, groups_var){
plt_histo <- ggplot(data_in, aes(x={{x_var}}, color = {{groups_var}})) + 
                  geom_histogram(alpha = 0.5, position="identity") +
                  theme(legend.position="bottom")
return(plt_histo)
}


#plot ecdf
ecdf_plt <- function(data_in, x_var, y_var, groups_var){
plt_ecdf <- ggplot(data_in, aes(x=x={{x_var}}, color = {{groups_var}})) + 
                                stat_ecdf() + 
                                        theme(legend.position="bottom")
return(plt_ecdf)                                    
}
