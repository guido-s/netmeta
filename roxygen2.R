##
## (1) Make R packages available
##
library(devtools)
library(roxygen2)


##
## (2) Create documentation file(s) in subdirectory testroxygen/man
##
# roxygenise("testroxygen") # Only consider files in subdirectory testroxygen/R
document("netmeta") # Also considers datasets in subdirectory netmeta/data


##
## (3) Build R package and PDF file with help pages
##
build("netmeta")
build_manual("netmeta")


##
## (4) Install R package
##
install("netmeta")


##
## (5) Check R package
##
check("netmeta")
