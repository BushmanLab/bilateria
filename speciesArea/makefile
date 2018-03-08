all: out/bayesSpeciesArea.pdf out/dadaSpeciesArea.pdf out/speciesAbundanceFit.pdf out/dadaSpeciesAbundanceFit.pdf

work/speciesAbund.Rdat: readData.R
	Rscript readData.R

work/dadaAbund.Rdat: readData.R
	Rscript readData.R

work/dadaRareN.csv: readData.R
	Rscript readData.R

work/rareN.csv: readData.R
	Rscript readData.R

work/rareN2.csv: readData.R
	Rscript readData.R

out/bayesSpeciesArea.pdf: speciesArea.R work/rareN.csv work/rareN2.csv
	Rscript speciesArea.R

out/dadaSpeciesArea.pdf: dadaSpeciesArea.R work/dadaRareN.csv
	Rscript dadaSpeciesArea.R

out/speciesAbundanceFit.pdf: sad.R work/speciesAbund.Rdat 
	#Rscript interferes with withTimeout
	R CMD BATCH sad.R

out/dadaSpeciesAbundanceFit.pdf: sad.R work/dadaAbund.Rdat
	Rscript sad.R
