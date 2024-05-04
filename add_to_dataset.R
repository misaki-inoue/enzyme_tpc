data_string <- "16.77,63.12
25.09,95.94
28.94,100.00
33.29,88.28
41.74,67.97
50.19,59.84"

test <- read.table(text = data_string, sep = ",", header = FALSE, col.names = c("Column1", "Column2"))
