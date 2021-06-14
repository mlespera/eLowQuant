# eLowQuant
eDNA Low copy number Quantification

This repository contains the Rmarkdown (.Rmd) , R (.R) and csv files used to compute the results in 'A Statistical Model for Calibration and Computation of Detection and Quantification Limits for Low Copy Number Environmental DNA samples' by Lesperance, Allison, Bergman, Hocking, Helbing, *eDNA*, 2021.  It also contains a shorter Rmarkdown file that can be used to analyze other data sets.  See the instructions in the Rmarkdown files.  These scripts are run using Rmarkdown within RStudio.  Install the latest versions of R and Rstudio.

Files:  
	- PoissonCalib-7Apr2021.Rmd (Rmarkdown file for publication computations)  
	- PoissonCalib-7Apr2021.pdf (rendered pdf from Rmd file)  
	- PoissonCalib-7Apr2021.html (rendered html from Rmd file)  
	- PoissonCalib-Functions-7April2021.R (R functions used by Rmd file)   
	- GEDWG_LOD_DATA3.csv (Klymus data)   
	- AssaySummary-ML.csv (Klymus results)  
	- eLowQuant.Rmd (Rmarkdown file for general use)  

Note:  To knit an Rmarkdown file to create a pdf file, your computer must have installed a version of latex or use the package tinytex.  If you do not have this, then you can knit to html.
