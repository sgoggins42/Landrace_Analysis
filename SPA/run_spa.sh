#bin/bash

##written by Shawn Goggins
## 10/09/15
# <spa_geno.df> is file containing 0, 1, 2 for number of minor alleles. 
# <spa_geno.df> SNPs must be biallelic 

# program available at: http://genetics.cs.ucla.edu/spa

./spa --gfile spa_geno.df --location-output spa_geno.df --model-output spa_geno.model

