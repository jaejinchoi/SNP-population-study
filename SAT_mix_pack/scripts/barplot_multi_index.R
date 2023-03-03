#!/usr/bin/Rscript --no-save
require(optparse)
require(dplyr)

opt_list = list(
	make_option(c("-n", "--name"), type="character", action="store", default="", help="Individual name containing file")
)


barplot_draw <- function(name_load_path, value_load_path, k_size)
{
  name_table <- read.table(name_load_path, header=F) 
  plot_table <- read.table(value_load_path, header=F)

  row.names(plot_table) <- name_table$V2
  item_size <- nrow(plot_table)

  #sort by the first value (the first column)
  #plot_table <- arrange(plot_table, V1)

  #sort by all columns available, and apply to row order
  plot_table <- plot_table[do.call(order, plot_table),]  
 
  save_path <- paste(value_load_path, ".png", sep='')
  
  #print(plot_table)
  #png(filename = save_path, width=700, height=300, pointsize=14) #png use pixels as size, open device
  png(filename = save_path, width = 12*item_size, height = 800, pointsize=12) #png use pixels as size, open device, resize figure width based on a number of individuals or items

  barplot(t(as.matrix(plot_table)), col=rainbow(k_size), las=2, cex.names=0.6, xlab=paste0("Individuals of ", item_size), ylab=paste0("Ancestry admixture using k = ", k_size))
  
  dev.off()

  #print(plot_table)
  #print(t(as.matrix(plot_table)))
}


opt <- parse_args(OptionParser(option_list=opt_list), positional_arguments = TRUE) #for trailing others
#args <- commandArgs(TRUE)

args <- opt$args
{
  name_load_path <- opt$options$name

  for (value_load_path in args)
  { 
    path_element <- as.matrix(unlist(strsplit(value_load_path, "[.]")))
    k_size <- strtoi(path_element[nrow(path_element)-1])
  
    barplot_draw(name_load_path, value_load_path, k_size)
  
  }
}
