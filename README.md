# netmeta: Network Meta-Analysis using Frequentist Methods
Official Git repository of R package **netmeta**

[![License: GPL (>=2)](https://img.shields.io/badge/license-GPL-blue)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
[![CRAN Version](https://www.r-pkg.org/badges/version/netmeta)](https://cran.r-project.org/package=netmeta)
[![GitHub develop](https://img.shields.io/badge/develop-3.1--1-purple)](https://img.shields.io/badge/develop-3.1--1-purple)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/netmeta)](https://cranlogs.r-pkg.org/badges/netmeta)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/netmeta)](https://cranlogs.r-pkg.org/badges/grand-total/netmeta)


## Authors

[Gerta Rücker](https://orcid.org/0000-0002-2192-2560),
Ulrike Krahn,
[Jochem König](https://orcid.org/0000-0003-4683-0360),
[Orestis Efthimiou](https://orcid.org/0000-0002-0955-7572),
[Annabel Davies](https://orcid.org/0000-0003-2320-7701),
[Theodoros Papakonstantinou](https://orcid.org/0000-0002-6630-6817),
[Guido Schwarzer](https://orcid.org/0000-0001-6214-9087)


## Contributors

[Theodoros Evrenoglou](https://orcid.org/0000-0003-3336-8058),
[Krzysztof Ciomek](https://orcid.org/0000-0002-2293-2146)


## Description

R package **netmeta** ([Balduzzi et al., 2023](https://www.doi.org/10.18637/jss.v106.i02)) provides frequentist methods for network meta-analysis and supports [Schwarzer et al. (2015)](https://link.springer.com/book/10.1007/978-3-319-21416-0), Chapter 8 "Network Meta-Analysis".

### Available network meta-analysis models

  - frequentist network meta-analysis ([Rücker, 2012](https://scholar.google.com/scholar?q=Rücker+2012+Network+meta-analysis+electrical+networks+and+graph+theory); [Rücker & Schwarzer, 2014](https://scholar.google.com/scholar?q=Rücker+Schwarzer+2014+Reduce+dimension+or+reduce+weights));

  - additive network meta-analysis for combinations of treatments
    ([Rücker, Petropoulou et al.,
    2020](https://doi.org/10.1002/bimj.201800167));

  - network meta-analysis of binary data using the Mantel-Haenszel
    method or the non-central hypergeometric distribution ([Efthimiou
    et al.,
    2019](https://scholar.google.com/scholar?q=Efthimiou+Rücker+Schwarzer+Higgins+Egger+Salanti+2019+Mantel-Haenszel+model)),
    or penalised logistic regression ([Evrenoglou et al.,
    2022](https://doi.org/10.1002/sim.9562)).


### Methods to present results of a network meta-analysis

  - network graphs ([Rücker & Schwarzer,
    2016](https://scholar.google.com/scholar?q=Rücker+Schwarzer+2016+Automated+drawing+of+network+plots+in+network+meta-analysis));

  - forest plots;

  - league tables with network meta-analysis results;
  
  - tables with network, direct and indirect estimates looking similar to the statistical part of a GRADE table for a network meta-analysis ([Puhan et al., 2014](https://scholar.google.com/scholar?q=puhan+schünemann+murad+2014+grade+network+meta-analysis)).


### Methods to rank treatments

  - rankograms and ranking by the Surface Under the Cumulative RAnking curve (SUCRA) ([Salanti et al., 2011](https://scholar.google.com/scholar?q=salanti+ades+ioannidis+2011+graphical+methods+multiple-treatment+meta-analysis));

  - ranking of treatments by P-scores (frequentist analogue of SUCRAs without resampling)
    ([Rücker & Schwarzer,
    2015](https://doi.org/10.1186/s12874-015-0060-8));

  - partial order of treatment rankings ('poset') and Hasse diagram
    for 'poset' ([Carlsen & Bruggemann,
    2014](https://scholar.google.com/scholar?q=Partial+order+methodology%3A+a+valuable+tool+in+chemometrics);
    [Rücker & Schwarzer,
    2017](https://scholar.google.com/scholar?q=Rücker+Schwarzer+2017+resolve+conflicting+rankings+of+outcomes+in+network+meta-analysis)).


### Methods to evaluate network inconsistency

  - split direct and indirect evidence to check consistency ([Dias et
    al.,
    2010](https://scholar.google.com/scholar?q=Checking+consistency+in+mixed+treatment+comparison+meta-analysis));

  - net heat plot and design-based decomposition of Cochran's Q
    ([Krahn et al., 2013](https://doi.org/10.1186/1471-2288-13-35)).


### Additional methods

 - contribution of direct comparisons to network estimates ([Papakonstantinou et al.,
   2018](https://doi.org/10.12688/f1000research.14770.3); [Davies et al., 2022](https://doi.org/10.1002/sim.9346))
 
  - importance of individual studies measured by reduction of precision if removed from network ([Rücker, Nikolakopoulou et al., 2020](https://doi.org/10.1186/s12874-020-01075-y))

  - 'comparison-adjusted' funnel plot ([Chaimani & Salanti,
    2012](https://scholar.google.com/scholar?q=Chaimani+Salanti+Using+network+meta-analysis+to+evaluate+the+existence+of+small-study+effects+in+a+network+of+interventions));
  
  - measures characterizing the flow of evidence between two
    treatments ([König et al.,
    2013](https://scholar.google.com/scholar?q=König+Krahn+Binder+2013+Visualizing+the+flow+of+evidence+in+network+meta-analysis+and+characterizing+mixed+treatment+comparisons)).


## Installation

### Current official [![CRAN Version](https://www.r-pkg.org/badges/version/netmeta)](https://cran.r-project.org/package=netmeta) release:
```r
install.packages("netmeta")
```

### Current [![GitHub develop](https://img.shields.io/badge/develop-3.1--1-purple)](https://img.shields.io/badge/develop-3.1--1-purple) release on GitHub:

Installation using R package
[**remotes**](https://cran.r-project.org/package=remotes):
```r
install.packages("remotes")
remotes::install_github("guido-s/netmeta",
  ref = "develop", build_vignettes = TRUE)
```


## How to cite netmeta?

[Balduzzi S, Rücker G, Nikolakopoulou A, Papakonstantinou T, Salanti G, Efthimiou O, Schwarzer G (2023): netmeta: An R package for network meta-analysis using frequentist methods. *Journal of Statistical Software*, **106**, 1-40](https://doi.org/10.18637/jss.v106.i02)

A BibTeX entry for LaTeX users is provided by

```
citation(package = "netmeta")
```


### Bug Reports:

You can report bugs on GitHub under
[Issues](https://github.com/guido-s/netmeta/issues).

or using the R command

```r
bug.report(package = "netmeta")
```

(which is not supported in RStudio).


## References

[Balduzzi S, Rücker G, Nikolakopoulou A, Papakonstantinou T, Salanti G, Efthimiou O, Schwarzer G (2023): netmeta: An R package for network meta-analysis using frequentist methods. *Journal of Statistical Software*, **106**, 1-40](https://doi.org/10.18637/jss.v106.i02)

[Carlsen L, Bruggemann R (2014): Partial order methodology: a valuable tool in chemometrics. *Journal of Chemometrics*, **28**, 226-34](https://scholar.google.com/scholar?q=Partial+order+methodology%3A+a+valuable+tool+in+chemometrics)

[Chaimani A, Salanti G (2012): Using network meta-analysis to evaluate the existence of small-study effects in a network of interventions. *Research Synthesis Methods*, **3**, 161-76](https://scholar.google.com/scholar?q=Chaimani+Salanti+Using+network+meta-analysis+to+evaluate+the+existence+of+small-study+effects+in+a+network+of+interventions)

[Davies AL, Papakonstantinou T, Nikolakopoulou A, Rücker G, Galla T (2022): Network meta-analysis and random walks. *Statistics in Medicine*, **41**, 2091-2114](https://doi.org/10.1002/sim.9346)

[Dias S, Welton NJ, Caldwell DM, Ades AE (2010): Checking consistency in mixed treatment comparison meta-analysis. *Statistics in Medicine*, **29**, 932-44](https://scholar.google.com/scholar?q=Checking+consistency+in+mixed+treatment+comparison+meta-analysis)

[Efthimiou O, Rücker G, Schwarzer G, Higgins J, Egger M, Salanti G
(2019): A Mantel-Haenszel model for network meta-analysis of rare
events. *Statistics in Medicine*, 1-21](https://scholar.google.com/scholar?q=Efthimiou+Rücker+Schwarzer+Higgins+Egger+Salanti+2019+Mantel-Haenszel+model)

[König J, Krahn U, Binder H (2013): Visualizing the flow of evidence in network meta-analysis and characterizing mixed treatment comparisons. *Statistics in Medicine*, **32**, 5414-29](https://scholar.google.com/scholar?q=König+Krahn+Binder+2013+Visualizing+the+flow+of+evidence+in+network+meta-analysis+and+characterizing+mixed+treatment+comparisons)

[Krahn U, Binder H, König J (2013): A graphical tool for locating inconsistency in network meta-analyses. *BMC Medical Research Methodology*, **13**, 35](https://doi.org/10.1186/1471-2288-13-35)

[Papakonstantinou T, Nikolakopoulou A, Rücker G, Chaimani A, Schwarzer G, Egger M, Salanti G (2018): Estimating the contribution of studies in network meta-analysis: paths, flows and streams. *F1000Research*, **7**, 610](https://doi.org/10.12688/f1000research.14770.3)

[Puhan MA, Schünemann HJ, Murad MH, Li T, Brignardello-Petersen R, Singh JA, Kessels AG, Guyatt GH, for the GRADE Working Group (2014): A GRADE Working Group approach for rating the quality of treatment effect estimates from network meta-analysis. *BMJ*, **349**, g5630](https://scholar.google.com/scholar?q=puhan+schünemann+murad+2014+grade+network+meta-analysis))

[Rücker G (2012): Network meta-analysis, electrical networks and graph theory. *Research Synthesis Methods*, **3**, 312-24](https://scholar.google.com/scholar?q=Rücker+2012+Network+meta-analysis+electrical+networks+and+graph+theory)

[Rücker G, Schwarzer G (2014): Reduce dimension or reduce weights? Comparing two approaches to multi-arm studies in network meta-analysis. *Statistics in Medicine*, **33**, 4353-69](https://scholar.google.com/scholar?q=Rücker+Schwarzer+2014+Reduce+dimension+or+reduce+weights)

[Rücker G, Schwarzer G (2015): Ranking treatments in frequentist network meta-analysis works without resampling methods. *BMC Medical Research Methodology*, **15**, 58](https://doi.org/10.1186/s12874-015-0060-8)

[Rücker G, Schwarzer G (2016): Automated drawing of network plots in network meta-analysis. *Research Synthesis Methods*, **7**, 94-107](https://scholar.google.com/scholar?q=Rücker+Schwarzer+2016+Automated+drawing+of+network+plots+in+network+meta-analysis)

[Rücker G, Schwarzer G (2017): Resolve conflicting rankings of outcomes in network meta-analysis: Partial ordering of treatments. *Research Synthesis Methods*, **8**, 526-36](https://scholar.google.com/scholar?q=Rücker+Schwarzer+2017+resolve+conflicting+rankings+of+outcomes+in+network+meta-analysis)

[Rücker G, Petropoulou M, Schwarzer G (2020): Network meta-analysis of multicomponent interventions. *Biometrical Journal*, **62**, 808-21](https://doi.org/10.1002/bimj.201800167)

[Rücker G, Nikolakopoulou A, Papakonstantinou T, Salanti G, Riley RD, Schwarzer G (2020): The statistical importance of a study for a network meta-analysis estimate. *BMC Medical Research Methodology*, **20**, 190](https://doi.org/10.1186/s12874-020-01075-y)

[Salanti G, Ades AE, Ioannidis JPA (2011): Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial. *Journal of Clinical Epidemiology*, **64**, 163-71](https://scholar.google.com/scholar?q=salanti+ades+ioannidis+2011+graphical+methods+multiple-treatment+meta-analysis)

[Schwarzer G, Carpenter JR and Rücker G (2015): *Meta-Analysis with R (Use R!)*. Springer International Publishing, Switzerland](https://link.springer.com/book/10.1007/978-3-319-21416-0)
