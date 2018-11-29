\name{netmeta-package}

\alias{netmeta-package}

\docType{package}

\title{netmeta: Brief overview of methods and general hints}

\description{
R package \bold{netmeta} provides frequentist methods for network
meta-analysis and supports Schwarzer et al. (2015), Chapter 8 on
network meta-analysis \url{http://meta-analysis-with-r.org/}.
}

\details{
  R package \bold{netmeta} is an add-on package for \bold{meta}
  providing the following meta-analysis methods:
  \itemize{
  \item frequentist network meta-analysis (function
    \code{\link{netmeta}}) based on Rücker (2012) and Rücker & Schwarzer
    (2014);
  \item net heat plot (\code{\link{netheat}}) and design-based
    decomposition of Cochran's Q (\code{\link{decomp.design}})
    described in Krahn et al. (2013);
  \item measures characterizing the flow of evidence between two
    treatments (\code{\link{netmeasures}}) described in König et
    al. (2013);
  \item ranking of treatments (\code{\link{netrank}}) based on
    frequentist analogue of SUCRA (Rücker & Schwarzer, 2015);
  \item partial order of treatment rankings (\code{\link{netposet}},
    \code{\link{plot.netposet}}) and Hasse diagram (\code{\link{hasse}})
    according to Carlsen & Bruggemann (2014) and Rücker & Schwarzer
    (2017);
  \item split direct and indirect evidence (\code{\link{netsplit}}) to
    check for consistency (Dias et al., 2010);
  \item league table with network meta-analysis results
    (\code{\link{netleague}});
  \item additive network meta-analysis for combinations of treatments
    (\code{\link{netcomb}}, \code{\link{discomb}} for disconnected
    networks) (Rücker et al., 2018);
  \item network meta-analysis of binary data using the Mantel-Haenszel
    or non-central hypergeometric distribution method
    (\code{\link{netmetabin}});
  \item{\sQuote{comparison-adjusted} funnel plot
    (\code{\link{funnel.netmeta}}) to assess funnel plot asymmetry in
    network meta-analysis (Chaimani & Salanti, 2012)}
  \item automated drawing of network graphs (\code{\link{netgraph}})
    described in Rücker & Schwarzer (2016).
  }
  
  Furthermore, functions and datasets from \bold{netmeta} are utilised
  in Schwarzer et al. (2015), Chapter 8 "Network Meta-Analysis",
  \url{http://meta-analysis-with-r.org/}.

  Type \code{help(package = "netmeta")} for a listing of R functions
  available in \bold{netmeta}.
  
  Type \code{citation("netmeta")} on how to cite \bold{netmeta} in
  publications.
  
  To report problems and bugs
  \itemize{
  \item type \code{bug.report(package = "netmeta")} if you do not use
    RStudio,
  \item send an email to Guido Schwarzer
    \email{sc@imbi.uni-freiburg.de} if you use RStudio.
  }

  The development version of \bold{netmeta} is available on GitHub
  \url{https://github.com/guido-s/netmeta}.
}

\references{
  Carlsen L, Bruggemann R (2014),
  Partial order methodology: a valuable tool in chemometrics.
  \emph{Journal of Chemometrics},
  \bold{28}, 226--34, DOI:10.1002/cem.2569
  
  Chaimani A & Salanti G (2012),
  Using network meta-analysis to evaluate the existence of small-study
  effects in a network of interventions.
  Research Synthesis Methods,
  \bold{3}(2), 161--176
  
  Dias S, Welton NJ, Caldwell DM, Ades AE (2010),
  Checking consistency in mixed treatment comparison meta-analysis.
  \emph{Statistics in Medicine},
  \bold{29}, 932--44

  König J, Krahn U, Binder H (2013),
  Visualizing the flow of evidence in network meta-analysis and
  characterizing mixed treatment comparisons.
  \emph{Statistics in Medicine},
  \bold{32}(30), 5414--29
  
  Krahn U, Binder H, König J (2013),
  A graphical tool for locating inconsistency in network
  meta-analyses.
  \emph{BMC Medical Research Methodology},
  \bold{13}, 35
  
  Rücker G (2012),
  Network meta-analysis, electrical networks and graph theory.
  \emph{Research Synthesis Methods},
  \bold{3}, 312--24
  
  Rücker G, Schwarzer G (2014),
  Reduce dimension or reduce weights? Comparing two approaches to
  multi-arm studies in network meta-analysis.
  \emph{Statistics in Medicine},
  \bold{33}, 4353--4369

  Rücker G, Schwarzer G (2015),
  Ranking treatments in frequentist network meta-analysis works
  without resampling methods.
  \emph{BMC Medical Research Methodology},
  \bold{15}, 58, DOI:10.1186/s12874-015-0060-8
  
  Rücker G, Schwarzer G (2016),
  Automated drawing of network plots in network meta-analysis.
  \emph{Research Synthesis Methods},
  \bold{7}, 94--107
  
  Rücker G, Schwarzer G (2017),
  Resolve conflicting rankings of outcomes in network meta-analysis:
  Partial ordering of treatments.
  \emph{Research Synthesis Methods},
  \bold{8}, 526--36, DOI:10.1002/jrsm.1270

  Rücker G, Petropoulou M, Schwarzer G (2018),
  Network meta-analysis of multicomponent interventions,
  Manuscript, submitted.

  Schwarzer G, Carpenter JR and Rücker G (2015),
  \emph{Meta-Analysis with R (Use-R!)}.
  Springer International Publishing, Switzerland.
  \url{http://www.springer.com/gp/book/9783319214153}
}

\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}, Gerta Rücker \email{ruecker@imbi.uni-freiburg.de}}
