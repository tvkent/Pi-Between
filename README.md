# Pi-Between (between population diversity vs recombination rate)

Using ANGSD and Python/R to calculate and plot pi-between in windows along the genome against recombination rate

### Background

This analysis borrows a method from [Brandvain and Sweigart et al.](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004410) to test for selection against gene flow using correlation between diversity between populations and recombination rate across the genome. These scripts should provide the necessary information to modify the process for different projects. These scripts are used in [this](https://github.com/SidBhadra-Lobo/Rice_project) project.
For both R scripts, a file containing recombination information is needed. The scripts use .RData information from [Corbett-Detig et al.](https://github.com/tsackton/linked-selection). Original R code and explanations from Jeff Ross-Ibarra can be found [here](http://rpubs.com/rossibarra/62904). Jeff's version uses minor allele frequency data from SNPs using VCFtools. This version uses whole genome NGS data including invariant sites to obtain minor allele frequency data using [ANGSD](http://popgen.dk/wiki/index.php/ANGSD).

### ANGSD

ANGSD needs to be run to obtain a .mafs file containing information on major and minor alleles. To obtain this, a version of this

	./angsd -bam -ref -anc -doMaf -doMajorMinor -rf -GL -out

needs to be run.  Reference and/or ancestral fastas may or may not be utilized. There are a few options for the `-doMaf` argument:
* 0 (Calculate persite frequencies '.mafs.gz')
* 1: Frequency (fixed major and minor)
* 2: Frequency (fixed major unknown minor)
* 4: Frequency from genotype probabilities
* 8: AlleleCounts based method (known major minor)  

The correct option for your data should be chosen. The Python and R code explained below use different `-doMaf` arguments, and need to be adjusted depending on which argument you choose. Future analyses for the rice gene flow project will be using `-doMaf 4` and the code presented here will be adjusted as needed.  
Regions may be chosen (i.e. calculate mafs for only a subset of chromosomes as in `-r 1:` or basepairs `-r chr1:1-10000`) or the entire genome may be used. For best results and no hang ups, it seems that suppyling a regions file with the desired chromosomes (even all chromosomes) works the best e.g. `-rf filename`. An example regions file can be found in __Scripts/__. A set of bams should be provided as a bamlist file (example in __Scripts/__).

### Python Version

In order to run large files, I translated the original R code into Python. This version uses the Python packages [numpy](http://www.numpy.org/) and [pandas](http://pandas.pydata.org/index.html).

This script takes .mafs.gz files from two populations, filters for sites with data from a minimum number of individuals and adds major allele frequency for each population. It then merges the two population data tables based on chromosome and position intersections, and calulates PiB. PiB is calulated as

	PiB = (Pop1 maf)*(Pop2 1-maf) + (Pop1 1-maf)*(Pop2 maf)

if major alleles at a given site are the same in both populations. If major alleles at a given site differ, PiB is calculated as

	PiB = (Pop1 maf)*(Pop2 maf) + (Pop1 1-maf)*(Pop2 1-maf)

to correct for major allele differences (i.e. we are comparing the same allele in both populations, and if this allele is segregating differently in each population, we must correct to compare the correct allele, not whatever is segregating as major in either).

From here, values are then binned into 1 Mb windows per chromosome and the average divergence per 1Mb is calculated.

Run this using __pib.sh__ from within the __Scripts/__ directory to take advantage of the file clean up (concatenating per chromosome files to one file containing the whole genome) to input into R. This code should work for per chromosome analysis (the provided loop over chromosome number) just make sure to add 1 to total chromosome number in the loop, as the range is non-inclusive.

To plot, use the __pibpy.R__ script in R. This script requires the packages data.table, dplyr, magrittr, and ggplot2.  

_Note:_ this version uses the `-doMaf 1` ANGSD argument in its current form. Header information in the script needs to be changed to work with different `-doMaf` arguments.

### R Version

This is the original code slightly modified from Jeff Ross-Ibarra's code and requires the packages data.table, dplyr, magrittr, and ggplot2. It is functionally identical to the python version explained above.  
Run this using __pib.R__ script in R.  

_Note:_ this version uses the `-doMaf 2` ANGSD argument in its current form. Header information in the script needs to be changed to work with different `-doMaf` arguments.
