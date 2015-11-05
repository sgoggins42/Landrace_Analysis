Landrace_Analysis

Altitude
	Fumi's folder of altitude data for landraces
	
Fst
	hiefstat.R
		documentation for how I partitioned landraces and ran the
		program hierfstat in order to get the basic.stats over them.
		'perloc' is a vector of data in the output from biostats
		that includes the Fst data.
	East_West_biostats.R
		documentation for how I used the hierfstat output to look
		at the differentiation of SNP allele frequencies between
		the East and West landraces.

KML
	RunKML.R
		Load it into your R session
		For the last step, edit the main(file) command to have your file
		in the parenthases. The file should be in the same format as the 
		803_landraces_KML.csv file:
			columns: Accession.ID/Latitude/Longitude/Altitude/Country

LD
	LD_chro2.R
		doumentation for how I calculated the complete LD between the 
		outlier chromosome 2 SNPs that can be identifed by finding the 
		95% of Fst outliers in EW_basicstats.csv that are on chro2.

SPA
	macx: subdirectory with (Yang et al. 2012) documentation
	Practice_SPA: subdirectory with made up dataframes to test spa
	Run_SPA: subdirectory with data_* and spa_*
		data_* files attempt to put in location data - not working
			run with script run_data.sh
		spa_* files have spa infer location - is working
			run with scrfipt run_spa.sh
		write_spa_geno.R
			documentation for how I wrote the .geno file format
		write_spa_loc.R
			documentation for how I wrote the .loc file format **
	spa
		copy of program spa run by beginning the line with ./spa
	Play_SPA.R
		running documentation of exploratory data analysis with spa
	Zhou_outlierSPA.txt
		list of outlier SNPs from Zhou's wild type barley samples

Worked_Datasets
	803_landraces_KML.csv
		Convereted dataset from Ana Poet's Barley_Landraces github in
		KML format
	EW_basicstats.csv
		Output of hierfstat basic.stats$perloc with F-statistics over 
		landraces partitioned by Longitude 48 E, or the Zagros Mountains
	GeneticMap_iSelect_9k.txt
		documentation from Poet's about SNP info
	Land_6152_SNPs_AB.txt
		Genotype data in "A/B" format
	Private_allele_ALLpops.txt
		Dataset of private alleles
	tgeno.csv
		Inverted dataframe of Land_6152_SNPs_AB.txt for convenience of
		merging and messing with dataframes
		
			
