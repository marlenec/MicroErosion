# MicroErosion

MicroErosion is a R function calculating the loss of phylogenetic or taxonomic diversity of
a metacommunity (gamma-diversity) due to the disappearance of local patches, based on a user-supplied scenario.

## Requirements
* R version 3.6.0 or higher. Former versions of R might cause problems with entropart and ape packages and their dependencies;
* Entropart R-package  version 1.6-1 or higher, and its dependencies [CRAN page](https://CRAN.R-project.org/package=entropart);
* ape R-package version 5.3 or higher, and its dependencies  [CRAN page](https://CRAN.R-project.org/package=ape).

## Input

* **metacom (required)**: Matrix of abundances or relative abundances of species (columns) in a set of N patches (rows). All values in matrix must be numeric.

* **tree (optional)**: Object of class hclust, phylo, phylog or PPtree corresponding to the phylogenetic tree to be used for computing phylogenetic diversity. The tips of the tree must match the column names in metacom. If not provided, only taxonomic diversity is computed.

* **q (default 0)**: weight given to species abundances relative to taxonomy or phylogeny information
For taxonomic diversity, q=0 will compute species richness-1; q=1: Shannon entropy; q=2: Simpson entropy.
If a tree is provided, phylogenetic diversity will be computed, with for q=0: Faith's PD; q=1: Allen's phylogenetic entropy; q=2: Rao's phylogenetic entropy.

* **nb_eq (default FALSE)**: whether or not diversity should be expressed as equivalent number of species. Will be ignored if q=0 and no tree is provided.
If nb_eq=T and a tree is provided, phylogenetic diversity will be independent from the height of the tree.

* **scenario (optional)**: user-defined probability of extinction of each patches. Must be a vector with names() corresponding to rownames(metacom). Probabilities do no need to sum to 1 over all patches. If not provided, neutral scenario (with a random order of patches extinction) will be performed.

* **r (default 1)**: number of erosion replicates required.

* **example (default FALSE)**: whether or not example data should be simulated. This example includes a random metacommunity table of 10 species and 10 patches, a simulated scenario of erosion and a random pure-birth tree containing the species. If example is TRUE, the arguments metacom, tree and scenario will be ignored and the erosion and diversity indices will be computed on the example data.


## Output

A list, containing 3 objects:

* **$Initial_Gamma_diversity**: single value with gamma diversity value before erosion

* **$Gamma_diversity**: Matrix with gamma-diversity calculated for a given number of patches that went extinct. The rows correspond to the successive number of patches lost (from 1 to N), and the columns correspond to the replicates (i.e. if r=10, there will be 10 columns.)

* **$Patches_Removed**: Matrix containing the name of the patches removed at each step of erosion scenario. Rows and columns are as in $Gamma_diversity.

> **Notes**:  
> Diversity and entropy indices are computed using `Tsallis()` and `ChaoPD()` functions provided in entropart R-package (https://cran.r-project.org/web/packages/entropart/index.html).
> Please note that computing phylogenetic diversity indices could be time consuming (e.g. for a tree with 5000 species, 2.5s per local assemblage (i.e. per step of erosion and per replicate) using a 2.2 GHz Intel Core i7).











