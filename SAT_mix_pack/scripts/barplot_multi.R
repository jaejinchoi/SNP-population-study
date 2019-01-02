#!/usr/bin/Rscript --no-save

barplot_draw <- function(load_path, k_size)
{
  plot_table <- read.table(load_path, header=F)
  save_path <- paste(load_path, ".png", sep='')
  
  png(filename = save_path, width=700, height=300, pointsize=14) #png use pixels as size, open device

  barplot(t(as.matrix(plot_table)), col=rainbow(k_size), xlab="Individual #", ylab="Ancestry")
  
  dev.off()
  
}



args <- commandArgs(TRUE)

#args <- '../lab_storage/r_project/test_cat'
#args <- '/home/jjc/Desktop/service_aquirement_asignment/papgi/admixture_linux-1.23/135017_snps/135017_snps.2.Q'

for (load_path in args)
{ 
  path_element <- as.matrix(unlist(strsplit(load_path, "[.]")))
  k_size <- strtoi(path_element[nrow(path_element)-1])
  
  barplot_draw(load_path, k_size)
  
}

