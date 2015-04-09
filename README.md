# Pi-Between
Using ANGSD and Python/R to calculate and plot pi-between in windows along the genome against recombination rate

### Background

For both R scripts, a file containing recombination information is needed. The scripts use .Rdata information from [Sackton et al.](link).  
These scripts should provide the necessary information to modify the process for different projects. These scripts are used in [this](https://github.com/SidBhadra-Lobo/Rice_project) project.

### ANGSD

[ANGSD](http://popgen.dk/wiki/index.php/ANGSD) needs to be run to obtain a .mafs file containing information on major and minor alleles.  
The code used in the R version of these scripts is

### Python Version

In order to run large files, I translated the original R code into Python. This version uses the Python packages [numpy](http://www.numpy.org/) and [pandas](http://pandas.pydata.org/index.html).  
Run this using __pib.sh__ from within the __Scripts/__ directory to take advantage of the file clean up to input into R.  
To plot, use the __pibpy.R__ script in R. This script requires the packages data.table, dplyr, magrittr, and ggplot2.  
_Note:_ this version uses the `-doMaf 1` ANGSD argument in its current form. Header information in the script needs to be changed to work with different `-doMaf` arguments.

### R Version

This is the original code slightly modified from [Jeff Ross-Ibarra's code](https://github.com/rossibarra), and requires the packages data.table, dplyr, magrittr, and ggplot2.  
Run this using __pib.R__ script in R.  
_Note:_ this version uses the `-doMaf 2` ANGSD argument in its current form. Header information in the script needs to be changed to work with different `-doMaf` arguments.

__For more information please see the [wiki](https://github.com/tvkent/Pi-Between/wiki).__
