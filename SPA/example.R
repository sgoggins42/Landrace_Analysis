myfile.geno <- sample(x = c(0:2), replace = TRUE, size = 20)
myfile.geno <- matrix(data=myfile.geno, nrow = 10, ncol = 2)
write.table(myfile.geno, sep = " ", file = "myfile.geno", quote = FALSE, row.names = FALSE, col.names = FALSE)
