range(na.omit(location.altitude_60fixed$altitude))
alt.df <- location.altitude_60fixed
alt.df$side <- 0
hist(alt.df$altitude)
alt.df$side <- ifelse(test = alt.df$altitude > 1000, yes = 1, no = 0)
altside.df <- alt.df[, c(1,6)]
alt.geno <- merge(x = altside.df, y = geno.df, by.x = "Accession.ID", by.y= "X")
alt_basicstats <- basic.stats(alt.geno[, -1], diploid = FALSE)
write.csv(alt_basicstats$perloc, file = "alt.FST.csv")
plot(alt.FST$Fst)
high_alt <- subset(alt.df, alt.df$altitude > 3000)
