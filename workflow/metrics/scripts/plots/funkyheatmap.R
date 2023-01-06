input_file <- snakemake@input$tsv
output_file <- snakemake@output$tsv

install.packages('funkyheatmap', repos = 'http://cran.us.r-project.org')
library(funkyheatmap)
library(data.table)

dt <- fread(input_file)
print(dt)

png(output_file)
funkyheatmap(dt)
dev.off()
