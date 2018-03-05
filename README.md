# bilateria
Code for generating figures and statistics for bilaterian animal microbiome project
The bulk of the code for this project can be found in bilateria.R.

# R libraries
You will need to have the following R libraries installed for various code chunks:

pander
ape
phytools
phylobase
vegan
ggplot2
dplyr
reshape2
pheatmap
dnar
Rtsne
lattice
grid
gridExtra
ggimage
phyloseq
eclectic
ade4

# Set work_dir
Replace "/home/kevin/projects/bilateria" with your working directory on line 24 of bilateria.R

# Define functions
Running the code under sections "Libraries", "Setup", and "Functions" is necessary before running any analysis code.

# Running analysis code
The code from each analysis code chunk (e.g. "Figure 1") can be run independently from each other.  Data is reset (re-loaded from OTU table and sample sheet) within each code chunk.
