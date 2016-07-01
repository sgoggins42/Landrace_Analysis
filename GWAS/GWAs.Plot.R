# No factors Setting
options(stringsAsFactors = FALSE)
# Import data:
#   Genotype data
#   geno.dat <- read.csv("~/Landrace_Analysis/Worked_Datasets/Fumi.Cleaned.Genotypes.csv")

# Outlier data:
Bio.this <- as.data.frame(Bioclim[,c("Taxa", "IC1")])
SNP.this <- as.data.frame(geno.dat[,c("Accession.ID", "SCRI_RS_219749")])
This.SNP.dat <- as.data.frame(merge(x = Bio.this, 
                                    y = SNP.this,
                                    by.x = "Taxa", by.y = "Accession.ID"))
This.SNP.dat$Taxa <- as.character(This.SNP.dat$Taxa)
box.dat <- boxplot(This.SNP.dat$IC1~This.SNP.dat$SCRI_RS_219749,
                   xlim = c(-1, 2), ylim = c(-4,4), 
                   main = "GWAS outlier Allele Freqencies")

par(new = TRUE)

plot(This.SNP.dat$IC1~jitter(This.SNP.dat$SCRI_RS_219749,2), 
     col = ifelse(test = This.SNP.dat$SCRI_RS_219749 == 0, yes = "red", no = "blue"),
     xlim = c(-1, 2), ylim = c(-4,4),
     main = "GWAS outlier",
     xlab  = "Minor allele count",
     ylab = "Temperature phenotype")
table(This.SNP.dat$SCRI_RS_219749)
