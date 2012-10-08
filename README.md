# ExpoRkit 

An R-interface to the Fortran package Expokit for matrix exponentiation
and a little more. 

## Installation from CRAN

Version 0.9 is available from CRAN. It can be installed from R.

    :::R
	install.packages(expoRkit)
	
Then load the package and check out a help file.

	:::R
	library(expoRkit)
	?expv
	
## Installation of the development version from source

The current development version can be installed from the source in this
repository using the `devtools` package.

	:::R
	install.packages("devtools")  ## if 'devtools' not installed
	library(devtools)
	install_bitbucket(repo = "exporkit", username = "nielsrhansen")

If you experience permission problems during installation, try
installing in a local directory, e.g.

	:::R
	install_bitbucket(repo = "exporkit", username = "nielsrhansen",
	                  args = "--library=~/local/R/library")
					  
provided that the directory exists. Then load the package
from the local library

	:::R
	library(expoRkit, lib.loc = "~/local/R/library")



					  
