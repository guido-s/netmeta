# netmeta: Network Meta-Analysis using Frequentist Methods
Official GitHub repository of R package **netmeta**

[![Build Status](https://travis-ci.org/guido-s/netmeta.svg?branch=master)](https://travis-ci.org/guido-s/netmeta)
[![CRAN Version](http://www.r-pkg.org/badges/version/netmeta)](https://cran.r-project.org/package=netmeta)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/netmeta)](http://cranlogs.r-pkg.org/badges/netmeta)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/netmeta)](http://cranlogs.r-pkg.org/badges/grand-total/netmeta)
[![Research software impact](http://depsy.org/api/package/cran/netmeta/badge.svg)](http://depsy.org/package/r/netmeta)


## Description

A comprehensive set of functions providing frequentist methods for
network meta-analysis:

  - frequentist network meta-analysis (Rücker, 2012);

  - net heat plot and design-based decomposition of Cochran's Q (Krahn
    et al., 2013);

  - measures characterizing the flow of evidence between two
    treatments (König et al., 2013);

  - ranking of treatments based on frequentist analogue of SUCRA
    (Rücker & Schwarzer, 2015);

  - partial order of treatment rankings ('poset') and Hasse diagram
    for 'poset' (Carlsen & Bruggemann, 2014);

  - split direct and indirect evidence to check consistency (Dias et
    al., 2010);

  - league table with network meta-analysis results;

  - automated drawing of network graphs (Rücker & Schwarzer, 2016).

Furthermore, R package **netmeta** provides functions and datasets to
support Schwarzer et al. (2015), Chapter 8 "Network Meta-Analysis",
http://meta-analysis-with-r.org/ .

### References

[Carlsen L, Bruggemann R (2014), Partial order methodology: a valuable tool in chemometrics. *Journal of Chemometrics*, **28**, 226--34.](https://scholar.google.de/scholar?q=Partial+order+methodology%3A+a+valuable+tool+in+chemometrics)
  
[Dias S, Welton NJ, Caldwell DM, Ades AE (2010). Checking consistency in mixed treatment comparison meta-analysis. *Statistics in Medicine*, **29**, 932--44.](https://scholar.google.de/scholar?hl=de&q=Checking+consistency+in+mixed+treatment+comparison+meta-analysis)

[König J, Krahn U, Binder H (2013). Visualizing the flow of evidence in network meta-analysis and characterizing mixed treatment comparisons. *Statistics in Medicine*, **32**, 5414-29.](https://scholar.google.de/scholar?q=König+Krahn+Binder+2013+Visualizing+the+flow+of+evidence+in+network+meta-analysis+and+characterizing+mixed+treatment+comparisons)

[Krahn U, Binder H, König J (2013), A graphical tool for locating inconsistency in network meta-analyses. *BMC Medical Research Methodology*, **13**, 35.](https://scholar.google.de/scholar?q=Krahn+Binder+König+2013+A+graphical+tool+for+locating+inconsistency+in+network+meta-analyses)

[Rücker G (2012), Network meta-analysis, electrical networks and graph theory. *Research Synthesis Methods*, **3**, 312-324.](https://scholar.google.de/scholar?q=Rücker+2012+Network+meta-analysis+electrical+networks+and+graph+theory)

[Rücker G, Schwarzer G (2015), Ranking treatments in frequentist network meta-analysis works without resampling methods. *BMC Medical Research Methodology*, *15*: 58](https://scholar.google.de/scholar?q=Rücker+Schwarzer+2015+Ranking+treatments+in+frequentist+network+meta-analysis+works+without+resampling+methods)

[Rücker G, Schwarzer G (2016), Automated drawing of network plots in network meta-analysis. *Research Synthesis Methods*, **7**: 94-107.](https://scholar.google.de/scholar?q=Rücker+Schwarzer+2016+Automated+drawing+of+network+plots+in+network+meta-analysis)

[Schwarzer G, Carpenter JR and Rücker G (2015), *Meta-Analysis with R (Use-R!)*. Springer International Publishing, Switzerland](http://www.springer.com/gp/book/9783319214153)



## Installation

### Current official [![CRAN Version](http://www.r-pkg.org/badges/version/netmeta)](https://cran.r-project.org/package=netmeta) release:
```r
install.packages("netmeta")
```

### Current development version on GitHub:
```r
install.packages("devtools")
devtools::install_github("guido-s/netmeta", ref = "develop")
```

Note, the development version of **netmeta** always has version
number 6.6-6 and date 6666-06-06.


### Bug Reports:

```r
bug.report(package = "netmeta")
```

The bug.report function is not supported in RStudio. Please send an
email to Guido Schwarzer <sc@imbi.uni-freiburg.de> if you use RStudio.
