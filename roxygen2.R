##
## (1) Make R packages available
##
library(devtools)
library(roxygen2)


##
## (2) Create documentation file(s)
##
document("../netmeta")


##
## (3) Build R package and PDF file with help pages
##
build("../netmeta")
build_manual("../netmeta")


##
## (4) Install R package
##
install("../netmeta")


##
## (5) Check R package
##
check("../netmeta")


##
## (6) Check examples
##
setwd("..")
run_examples("netmeta", fresh = TRUE,
             run_dontrun = TRUE, run_donttest = TRUE)
warnings()


##
## (6) Check R package (with donttest examples)
##
##setwd("netmeta")
##check("../netmeta", run_dont_test = TRUE)
