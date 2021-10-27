## netmeta, version 2.0-1 (2021-10-27)

### User-visible changes

* print.netsplit(), forest.netsplit():
  - sorting by direct evidence proportion possible (argument
    'sortvar')
	
* print.netmeta(), print.netcomb():
  - print number of pairwise comparisons after number of studies

### Internal changes

* pairwise():
  - keep original order of treatments


## netmeta, version 2.0-0 (2021-10-12)

### Major changes

* Behaviour of print and print.summary functions switched (to be in
  line with other print and print.summary functions in R)

* Annabel Davies <annabel.davies@manchester.ac.uk> is a new co-author
  of R package **netmeta**

* Random walk algorithm implemented to estimate network contributions
  [(Davies et al., 2021)](https://arxiv.org/abs/2107.02886)

* Calculation of standardised mean differences and corresponding
  standard errors in pairwise() is based on [Crippa & Orsini
  (2016)](https://doi.org/10.1186/s12874-016-0189-0), equations (4)
  and (5), providing consistent treatment estimates and standard
  errors for multi-arm studies
  
* By default, reference group is defined by first treatment in network
  meta-analysis

* Renamed arguments:
  - 'fixed' (instead of 'comb.fixed')
  - 'random' (instead of 'comb.random')
  - 'level.ma' (instead of 'level.comb')

* New R function hatmatrix() to derive hat matrices

* New R function netpairwise() to conduct pairwise meta-analyses for
  comparisons with direct evidence

* New R function netcomplex() to calculate effect of arbitrary complex
  interventions in component network meta-analysis

* New R function netcomparison() to calculate comparison effects of
  two arbitrary complex interventions in component network
  meta-analysis

* I2 from pairwise comparisons shown in forest plot with direct
  evidence

* In component network meta-analysis, individual component names
  instead of full treatment names can be abbreviated for complex
  interventions with more than one component

* Do not stop calculations if standard errors of multi-arm studies are
  inconsistent (instead only check for positive variance estimates of
  single treatment arms)

* Generic functions are not exported

### Bug fixes

* netmeta():
  - fix error in calculation of between-study variance using REML or
    ML method

* forest.netbind(), forest.netcomb(), forest.netmeta(),
  forest.netsplit():
  - do not print empty row above descriptions if arguments
    'label.left' or 'label.right' is used

* netcomb(), discomb():
  - value provided for argument 'inactive' must be a single treatment
    component not a combination of treatment components (check was
    missing)

* funnel.netmeta():
  - works now with numeric treatment labels

### User-visible changes

* New functions hatmatrix() and print.hatmatrix() to derive hat
  matrices

* New functions netpairwise(), forest.netpairwise(),
  print.netpairwise(), summary.netpairwise() and
  print.summary.netpairwise() to conduct pairwise meta-analyses for
  all comparisons with direct evidence

* netcontrib():
  - new argument 'method' to select method to estimate network
    contributions
  - new logical argument 'hatmatrix.F1000' to specify whether hat
    matrix given in [Papakonstantinou et
    al. (2018)](https://doi.org/10.12688/f1000research.14770.3) should
    be used
  - argument 'nchar.trts' from netmeta object considered in printouts

* decomp.design(), print.decomp.design() print.netsplit(),
  rankogram(), print.rankogram(), plot.rankogram():
  -  new argument 'nchar.trts' to print abbreviated treatment names

* rankogram():
  - new list elements 'cumrank.matrix.fixed' and
    'cumrank.matrix.random' with cumulative ranking probabilites
  - new argument 'cumulative.rankprob' to show cumulative ranking
    probabilites
  - new argument 'trts' to specify subset of treatments

* plot.rankogram():
  - new argument 'pooled' to select results from fixed effect or
    random effects model
  - new argument 'trts' to specify subset of treatments
  - new argument 'sort' to specify order of treatments
  - new argument 'cumulative.rankprob' to show cumulative ranking
    probabilities

* print.rankogram():
  - new argument 'cumulative.rankprob' to show cumulative ranking
    probabilites
  
* netrank():
  - new arguments 'fixed' and 'random'
  - new list elements 'fixed' and 'random'

* netsplit():
  - new argument 'order' to specify order in comparisons; see help
    page of funnel.netmeta()

* forest.netsplit():
  - print I2 from pairwise comparisons if argument direct = TRUE

* forest.netsplit(), print.netsplit():
  - new argument 'only.reference' to only print comparisons with
    reference group
  - new argument 'sortvar' to sort comparisons

* pairwise():
  - keep all variables from original dataset if variable is missing
    for one group but not the other (so far, only the first variable
    was kept)

* summary.netcomb(), print.summary.netcomb(), print.netcomb():
  - new argument 'show.combs'

* netcomb():
  - new argument 'nchar.comps' to abbreviate component names

* discomb(), print.netcomb(), print.summary.netcomb():
  - new argument 'nchar.comps' replaces argument 'nchar.trts'

* New auxiliary function comps() to create unique comparison labels
  with abbreviated treatment names

### Internal changes

* netmeta():
  - new list elements 'seTE.adj.fixed' and 'seTE.adj.random with
    adjusted standard errors under fixed effect and random effects
    model (content of 'seTE.adj.fixed' is identical to the previously
    existing element 'seTE.adj')
  - new list elements 'H.matrix.fixed' (replacing list element
    'H.matrix') and 'H.matrix.random'
  - new list elements 'Q.direct', 'tau2.direct', 'tau.direct' and
    'I2.direct' with information on the between-study heterogeneity of
    direct comparisons

* netcomb(), discomb():
  - new list elements 'L.matrix.fixed', 'Lplus.matrix.fixed',
    'H.matrix.fixed', 'L.matrix.random', 'Lplus.matrix.random' and
    'H.matrix.random'
	
* compsplit():
  - remove leading and trailing blanks
	
* setref():
  - argument 'reference.group' can be a vector instead of a single value

* nma.additive():
  - new list elements 'L.matrix', 'Lplus.matrix' and 'H.matrix'

* New internal function hatmatrix.aggr() to calculate the aggregated
  hat matrix

* Internal function contribution.matrix() moved to
  contribution.matrix.tpapak()

* New internal functions contribution.matrix.davies()

* New internal function compos() to abbreviate component names

* New internal function createC.full() to create full C matrix for all
  k-fold combinations of n components

* New internal function charfac()

* New wrapper function contribution.matrix()

* New argument 'aggr' in createB() to calculate aggregated B matrix


## netmeta, version 1.5-0 (2021-06-28)

### Major changes

* Theodoros Papakonstantinou <dev@tpapak.com> is a new co-author of R
  package **netmeta**

* Rankograms added

* Surface under the cumulative ranking (SUCRA) can be calculated using
  resampling methods

* Method by [Papakonstantinou et
  al. (2018)](https://doi.org/10.12688/f1000research.14770.3) to
  estimate the contribution of studies to network meta-analysis
  implemented

### Bug fixes

* netgraph.netmeta():
  - use of arguments seq = "optimal" and srt.labels = "orthogonal"
    resulted in wrong rotation of treatment labels

### User-visible changes

* New functions rankogram(), print.rankogram() and plot.rankogram()
  for rankograms

* New functions netcontrib() and print.netcontrib() to calculate
  network contributions

* netrank():
  - new argument 'method' to choose between P-scores and SUCRAs
  - can be used with an R object created with rankogram()

* forest.netmeta():
  - argument 'digits.Pscore' renamed to 'digits.prop' (as this
    argument is also used for SUCRAs)
  - new argument 'nsim' to specify number of simulations for SUCRAs

* netgraph.netmeta():
  - by default, do not mark multi-arm studies (argument multiarm = FALSE)
  - by default, use inverse standard error of random effects estimates
    for argument 'thickness' if random effects network meta-analysis
    was conducted

* netsplit(), forest.netsplit():
  - new option "reference.only" for argument 'show'

* discomb(), netcomb():
  - new check for unidentifiable components implemented
  - new argument 'details.chkident' to print more details on
    unidentifiable components


## netmeta, version 1.4-0 (2021-05-11)
    
### Major changes

* Restricted maximum likelihood and maximum likelihood estimator for
  between-study heterogeneity implemented by calling rma.mv() from R
  package **metafor** internally

* R function pairwise() can be used to generate reduced data set with
  basic comparisons [Rücker & Schwarzer
  (2014)](https://doi.org/10.1002/sim.6236)

* R function discomb() can be used with an R object created with pairwise()

### Bug fixes

* netmeta():
  - use of ginv() instead of solve() to calculate pseudo inverse of
    Laplace matrix resulted in wrong results for some extreme network
    structures (bug was introduced in **netmeta**, version 1.3-0)

* netgraph.netmeta():
  - highlighting more than one comparison (argument 'highlight')
    resulted in an error if argument 'col.highlight' was not of same
    length

### User-visible changes

* netmeta():
  - new argument 'method.tau' to select estimation method for
    between-study variance
  - new arguments 'sd1', 'sd2', 'time1', and 'time2' used to calculate
    variance-covariance matrix for REML and ML estimator of
    between-study heterogeneity
  - do not calculate leverages for multi-arm studies

* netconnection():
  - first argument can be of class 'pairwise'

* pairwise():
  - new arguments 'reference.group' and 'keep.all.comparisons'
  - first argument can be of class 'pairwise'

* netgraph.netconnection():
  - print subnetworks in different colors

* forest.netmeta():
  - new argument 'equal.size' to determine whether square size should
    be equal or proportional to precision of estimates

* forest.netbind(), forest.netsplit():
  - new default for argument 'equal.size', i.e., square sizes are equal

* forest.netbind():
  - new argument 'subset.treatments' to select treatments shown in
    forest plot

* netmeasures():
  - return direct evidence proportion for single design with two
    treatments

* netheat():
  - all designs are shown in net heat plot by default

### Internal changes

* New internal function invmat() to calculate inverse of matrix

* New internal function calcV() to calculate variance-covariance
  matrix used as input to rma.mv()

* nma.ruecker(), multiarm():
  - use solve() instead of ginv() to calculate inverse of matrix L

* netmeta():
  - new list element 'comparisons' with information on direct comparisons

* netconnection():
  - new list elements 'comparisons' and 'subnet.comparisons' with
    information on direct comparisons

* netimpact():
  - consider values for arguments 'tol.multiarm' and 'tol.multiarm.se'
    from network meta-analysis object

* nma.krahn():
  - new argument 'reference.group' to specify the reference group


## netmeta, version 1.3-0 (2021-01-15)

### Major changes

* Treatment labels can be rotated in network graphs

* Additional checks whether net heat plot is feasible

* Calculate direct evidence proportion and indirect treatment
  estimates for network meta-analyses created with netmetabin()

* Network graph for objects created with netconnection(), netcomb() or
  discomb()

* In summaries, print z- and p-values for network estimates if a
  reference group is defined

### Bug fixes

* netmetabin():
  - not overall heterogeneity / inconsistency statistic calculated for
    Mantel-Haenszel method, but, the overall inconsistency statistic
    (treatment effects are aggregated within designs)

* decomp.design(), netheat():
  - do not conduct decomposition of designs for objects created with
    netmetabin(); reported decomposition and net heat plot referred to
    reanalysis using netmeta()

* netmeasures():
  - works for network meta-analyses created with netmetabin();
    reported results referred to reanalysis using netmeta()

* netsplit():
  - SIDDE method works for network meta-analysis created without
    argument 'data' or empty networks after dropping a design
    
### User-visible changes

* print.summary.netmeta():
  - output states "test of heterogeneity" for a single design and
    "test of inconsistency" for Mantel-Haenszel method
  
* netrank():
  - can be used with network meta-analysis object created with
    netcomb()

* netmeta():
  - new argument 'small.values' passed on to netrank(),
    forest.netmeta() and netposet()

* netgraph.netmeta():
  - new argument 'srt.labels' to rotate treatment labels
  
* plot.netrank():
  - get rid of warnings 'Undefined global functions or variables'

* netbind():
  - argument '...' can be a single list of network meta-analysis objects
  
* print.netmeta():
  - new arguments 'truncate' and 'text.truncate' to only show
    selection of individual study results (useful for very long
    printouts or to only show individual results of multi-arm
    studies))
  - print z- and p-values for tests of overall effect against a
    reference group (argument 'reference.group')
  
* print.netconnection():
  - new argument 'distance' in order to print the distance matrix
  - by default, do not print the distance matrix
  
* netsplit():
  - new arguments 'comb.fixed', 'comb.random' and 'backtransf'

### Internal changes

* Rename list elements starting with 'zval.' to 'statistic.'

* netcomb():
  - export covariance matrices 'Cov.fixed' and 'Cov.random'

* netmeta(), netcomb():
  - export study design for each pairwise comparison (list element
    'design')

* discomb():
  - export study designs (list element 'design'), list of designs
    ('design') and number of designs ('d')

* Argument 'single' renamed to 'length' in calls of chkchar(),
  chkcolor(), chklevel() and chknumeric()

* New internal function is.zero() to determine whether a small number
  is essentially zero (i.e., whether *abs(x) < 10 *
  .Machine$double.eps*)

* New internal function updateversion() to update older netmeta
  objects

* New internal function designs() to determine study designs

* New internal function netheat-internal() with auxiliary functions
  for netheat()

* nma.ruecker(), nma.additive():
  - use ginv() to calculate inverse of matrix L
  - set small numbers which are essentially zero to zero in matrix
    Lstar


## netmeta, version 1.2-1 (2020-04-15)

### Major changes

* Function netimpact() can be used with network meta-analysis objects
  created with netmetabin()

* No warning is printed that treatments within comparisons have been
  re-sorted in increasing order (which has been rather a note than a
  warning)

### Bug fixes

* netmetabin():
  - bug fix for studies with 0 events and 0 participants
    
* discomb():
  - return standard errors of individual studies in correct order
    (list element 'seTE')

* netbind():
  - function can be used with network meta-analysis objects created
    with netmetabin()
    
* as.data.frame.netmeta():
  - correct printout if number of studies is equal to the number of
    treatments to the power of 2

* summary.netmeta() exported

### User-visible changes

* netmeta(), discomb():
  - do not print a warning / note that treatments within comparisons
    have been re-sorted in increasing order
  - check whether all studies provide missing estimates and standard
    errors (and stop execution with an informative error message)

* netmeta(), netcomb(), discomb():
  - export the design matrix (new list element 'X.matrix')

* discomb():
  - export the number of subnetworks

* forest.netsplit():
  - print an informative warning (instead of an obscure error
    message) if no comparison in the network provides the requested
    combination of direct and indirect evidence

* Use Markdown for NEWS

### Internal changes

* prepare():
  - call metagen() with argument 'method.tau.ci = ""' to suppress
    calculation of confidence interval for tau2


## netmeta, version 1.2-0 (2019-12-20)

### Major changes

* In multi-arm studies, negative weights (resulting from slightly
  inconsistent standard errors) contributing less than 0.1% to the
  weight of a study are set to 0

* Separate consistency tolerances can be specified for treatment
  estimates and standard errors (which are consistent by design in
  multi-arm studies, however, can be inconsistent due to rounding of
  treatment estimates and standard errors)

* New default for consistency tolerance: 0.001 instead of 0.0005

* Print tau in addition to tau2 in outputs

* Print confidence interval for I2 in outputs

### User-visible changes

* netmeta(), discomb():
  - check for numeric values in arguments 'TE' and 'seTE'
  - new argument 'tol.multiarm.se' to check standard errors in
    multi-arm studies

* netmetabin():
  - check for numeric values in arguments 'event1', 'n1', 'event2',
    and 'n2'
  - new argument 'tol.multiarm.se' to check standard errors in
    multi-arm studies

* print.summary.netmeta(), print.summary.netcomb():
  - new argument 'digits.tau' to specify number of digits for tau
  - new arguments 'text.tau2', 'text.tau', and 'text.I2' to change
    text printed to identify respective heterogeneity measure
    
* Revision of help page examples to reduce runtime below 5 seconds
  (CRAN requirement)

### Internal changes

* multiarm():
  - negative weights contributing less than 0.1% to the weight of a
    study are set to 0

* chkmultiarm():
  - use same order of arguments as in netmeta()
  - new argument 'tol.multiarm.se' to check standard errors

* nma.ruecker(), nma.additive():
  - return lower and upper confidence limits for I2


## netmeta, version 1.1-0 (2019-08-06)

### Major changes

* New functions netimpact(), netgraph.netimpact(), print.netimpact()
  to measure the impact of individual studies to network estimates

* New function plot.netrank() to produce image plot of P-scores

* Equivalence limits can be added to forest plots

* Use **roxygen2** for development of R package **netmeta**

### Bug fixes

* netmeta():
  - no error if argument 'studlab' is missing
  - tackle numerical problems with zero treatment arm variances -
    actually by changes in internal function multiarm()
  - calculate the correct number of patients and events for each
    treatment arm in networks with multi-arm studies

* netmetabin():
  - only drop single treatment arms without any events in multi-arm
    studies instead of the complete design (which is the correct
    behaviour only for two-arm studies)
  - keep the original design for multi-arm studies if a single
    treatment arm without any events has been dropped, e.g., the
    design X:Y:Z is not changed to X:Z if no events were observed in
    treatment arm Y in any multi-arm study
  - no error if argument 'event1' is not an R object created with
    pairwise() and argument 'data' is not used
  - no error in printing of network meta-analysis object if argument
    'event1' is not an R object created with pairwise()
  - calculate the correct number of patients and events for each
    treatment arm in networks with multi-arm studies

* netgraph():
  - argument 'seq' equal to "optimal" works for a single design
  - multi-arm studies are highlighted if argument 'labels' is used
  - comparisons defined by argument 'highlight' are marked in the
    network graph if argument 'labels' is used

* pairwise():
  - no error if function is used without argument 'data'

### User-visible changes

* netgraph() is a generic function and original function netgraph()
  renamed to netgraph.netmeta()

* netgraph.netmeta():
  - new argument 'bg.points'
  - new argument 'scale.highlight' which replaces argument
    'lwd.highlight' (which actually was ignored)
  - arguments 'col.highlight' and 'scale.highlight' can be vectors
    with same length as argument 'highlight'
  - print highlighted comparisons in the background of the graph

* forest.netmeta():
  - new argument 'labels' to provide treatment names - similar to
    netgraph()
  - labels set automatically for columns "pscore", "k" and
    "prop.direct" in argument 'leftcols' and 'rightcols'

### Internal changes

* multiarm():
  - use ginv() to calculate inverse of matrix Lt
  - set very small negative values in matrix W equal to zero
  - object k instead of n denotes number of treatment arms
    (same notation as Rücker & Schwarzer, 2014)
  - more liberal cutpoint to set negative variances to zero

* chkmultiarm():
  - more liberal check for negative variances

* new auxiliary internal functions for forest plots


## netmeta, version 1.0-1 (2019-01-02)

### User-visible changes

* Revision of help page examples
  - to reflect changes in netmeta, versions 0.9-8 and 1.0-0
  - to reduce runtime below 5 seconds (CRAN requirement)


## netmeta, version 1.0-0 (2018-12-20)

### Major changes

* Orestis Efthimiou <oremiou@gmail.com> is a new co-author of R
  package **netmeta**

* New function netmetabin() for network meta-analysis of binary data
  using the Mantel-Haenszel method or the non-central hypergeometric
  distribution

* New function funnel.netmeta() for 'comparison-adjusted' funnel
  plots

* New function forest.netcomb() for forest plots of additive network
  meta-analysis

* New functions netbind(), forest.netbind(), and print.netbind() to
  combine network meta-analysis objects and to generate a forest
  plot with results of several network meta-analyses

* Separate Indirect from Direct Design Evidence (SIDDE) method
  implemented in netsplit()

* All network comparisons can be included in a forest plot

* Circular network graphs with minimal number of crossings available

* Default setting for levels of confidence intervals can be
  specified using R function settings.meta(); the default is still
  0.95, i.e., 95% confidence intervals are computed

* New datasets: Gurusamy2011 and Dong2013

### User-visible changes

* netmeta():
  - new list elements:
    - 'k.trts' with number of studies evaluating a treatment
    - 'n.trts' with number of observations receiving a treatment
    - 'events.trts' with number of events observed for a treatment
    - 'n.matrix' with number of observations in direct comparisons
    - 'events.matrix' with number of events in direct comparisons
    - 'designs' with treatment designs
  - correct entry for list element 'designs' for single design
  - only print information on studies with missing treatment effect
    or standard error if argument 'warn' is TRUE

* pairwise():
  - keep original order of studies
  - arguments 'event' and 'n' can be used for generic meta-analysis
    method based on arguments 'TE' and 'seTE'
  - argument 'n' can be used for meta-analysis with count outcomes
    (based on arguments 'event' and 'time')
  - infinite treatment estimates and standard errors set to NA
  - by default, do not print warnings if comparisons will not be
    included in network meta-analysis
  
* netsplit():
  - new argument 'method' to choose approach to split direct and
    indirect evidence

* forest.netmeta():
  - argument 'reference.group' can be a vector in order to include
    several / all network comparisons in a forest plot
  - (invisibly) returns a data frame with information used to
    produce the forest plot

* forest.netsplit() and print.netsplit():
  - argument 'showall' replaced with 'show' (see help pages)

* forest.netsplit():
  - show prediction intervals as colored bars
  - forest plot with layout 'subgroups by comparisons':
    omit treatment estimates in rows with prediction intervals if
    network estimates are also shown

* netcomb() and discomb():
  - argument 'seq.components' replaced with 'seq.comps'

* discomb():
  - arguments 'reference.group' and 'baseline.reference'

* netposet():
  - allow for ties in rankings
    
### Bug fixes

* netsplit():
  - do not reorder a treatment comparison if reference treatment
    (argument 'reference.group') is part of a combined treatment,
    e.g., for reference.group = "A" and treatment comparison "A + B"
    vs "C", the comparison will not be reordered as "C" vs "A + B"

* forest.netsplit():
  - show correct prediction intervals

* netgraph():
  - arguments 'iterate' and 'allfigures' ignored if argument 'seq'
    is equal to "optimal"
    
* netsplit():
  - consider argument 'digits' for treatment estimates, i.e., do not
    always round to two digits

### Internal changes

* netmeta():
  - keep original order of studies in list element 'data'

* netcomb() and discomb():
  - check if argument 'sep.components' is a single character
  - similar list elements as for R objects created with netmeta();
    especially matrices with all network estimates added to output,
    see, for example, new list elements 'TE.fixed' and 'TE.random'

* Internal function nma.additive():
  - calculate matrices with all network estimates

* netleague() is a generic function

* Help pages:
  - examples added for forest.netsplit()


## netmeta, version 0.9-8 (2018-03-23)

### Major changes

* Main function netmeta():
  - by default, results for fixed effects and random effects network
    meta-analysis are reported (only results for fixed effects model
    were reported in older versions)
  - keep dataset used to conduct network meta-analysis
  - number of events and number of observations can be provided (and
    are considered from R objects created with pairwise() function)

* Function pairwise():
  - all variables from the original dataset are kept in the output
    dataset

* League tables
  - show the direct treatment estimates from pairwise comparisons in
    the upper triangle if the table is created for a single network
    meta-analysis
  - report pairwise comparisons of the treatment in the row versus
    the treatment in the column in the lower triangle and column
    versus row in the upper triangle (common presentation for
    network meta-analyses)

* Network graphs are much more flexible, e.g., color of lines /
  edges can be specified for each direct pairwise comparison

* New function discomb() for disconnected networks sharing at least
  one common treatment component to apply the additive network model
  for combinations of treatments

* Treatment separator (argument 'sep.trts') can be special character
  from regular expressions

* Between-study variance tau-squared is reported as NA instead of 0
  in networks without heterogeneity / inconsistency

### User-visible changes

* netmeta():
  - settings for printing of results are defined by settings.meta(),
    i.e., new default is to print results for both fixed effects and
    random effects model
  - new arguments 'event1', 'event2', 'n1', and 'n2' to provide
    number of events and observations for the two treatment groups
  - new argument 'keepdata' to choose whether original dataset
    should be part of network meta-analysis object    
  - new list element 'data' to keep dataset used to conduct network
    meta-analysis (if argument 'keepdata' is TRUE)

* netgraph():
  - argument 'col' can be a matrix, e.g., created with netmatrix(),
    to specify the color of lines / edges for each direct pairwise
    comparison
  - arguments 'col.points', 'cex.points', and 'pch.points' can be
    used to specify the color, size, and plotting symbol for each
    treatment separately
  - new argument 'adj' to specify the adjustment of treatment labels
  - new argument 'pos.number.of.studies' to specify the position of
    the number of treatments on the edges
  - (invisibly) returns a list with information on nodes and edges
    (position, color, etc.) used to produce the network graph

* pairwise():
  - print correct labels in error message for studies with duplicate
    treatments
 - all variables from the original dataset are kept in the output
   dataset

* print.netmeta():
- print uniquely abbreviated treatment names

* netconnection() and print.netconnection():
  - new argument 'sep.trts' to print abbreviated treatment names

* netleague():
  - new argument 'big.mark' to specify character printed as
    thousands separator, e.g., big.mark = "," will result in
    printing of 1,000 for the number 1000
  - new argument 'text.NA' to label missing values, i.e., for
    pairwise comparisons without direct evidence
  - new argument 'direct' to print direct treatment estimates (if
    argument 'y' is not missing)

* print.netsplit():
  - new argument 'legend' to suppress printing of the legend

* print.summary.netcomb():
  - new arguments 'digits.tau2' and 'digits.I2'

* New auxiliary function netmatrix() to create a matrix with
  additional information for pairwise comparisons (e.g., risk of
  bias assessment)

### Bug fixes

* netmeta():
  - list elements 'P.fixed' and 'P.random' contained wrong values
    for network meta-analyses with a single design which resulted in
    an error using forest.netmeta()

* netcomb():
   - use correct between-study variance in random effects additive
     network meta-analysis model
   - use correct order of studies to calculate treatment estimates
   - can be used with single pairwise comparison

* print.summary.netcomb():
  - use correct (abbreviated) names for treatment components

### Internal changes

* new internal functions compmatch() and compsplit() for argument
  'sep.trts' taking special character from regular expression into
  account, e.g., sep.trts = "."

* new internal functions bySummary() which is used in netmeta(),
  netmatrix(), and pairwise()

* createC() can be used with netconnection() objects (for
  disconnected networks)

* netcomb.netmeta() has been renamed to netcomb()

* netmeta():
  - report I2 as value between 0 and 1 (instead of 0 and 100)
  - new list element 'trts' (character vector with treatment names)

* Internal function nma.ruecker():
  - use internal function isquared() from package **meta** to
    calculate I2

* Internal function nma.additive():
  - calcuate between-study variance tau2 and heterogeneity statistic
    I2

* Internal function prcombs():
  - new argument 'seq' to order treatments


## netmeta, version 0.9-7 (2017-12-06)

### Major changes

* Version of R package **meta** must be larger or equal 4.9-0

* New function to produce forest plots with network, direct, and
  indirect evidence

* New function to estimate additive network meta-analysis for
  combinations of treatments

* New default in function netsplit(): treatment comparisons are
  selected from upper treatment estimates matrix, i.e., comparisons
  are "A vs B", "A vs C", and "B vs C" for treatments "A", "B", and
  "C" instead of "B vs A", etc.

* Zero treatment arm variance in multi-arm studies results in a
  warning instead of an error message

* League tables can be exported as CSV or Excel file

* New argument 'backtransf' indicating whether estimates should be
  back-transformed in printouts and plots, e.g., to show results as
  odds ratios instead of log odds ratios

* P-values can be printed in scientific notation

* P-values equal to 0 are actually printed as "0" instead of
  "< 0.0001"

* Thousands separator can be used in printouts and forest plots for
  large numbers

### User-visible changes

* new functions:
  - forest.netsplit() to produce forest plots with direct and
    indirect evidence
  - print.netleague() to print league table
  - treats() to create uniquely abbreviated treatment names
  - netcomb() and netcomb.netmeta() to estimate additive network
    meta-analysis models for combinations of treatments   
  - summary.netcomb(), print.netcomb(), and print.summary.netcomb()
    to print (summaries of) netcomb objects

* print.decomp.design(), print.netmeta(), print.netsplit(), and
  print.summary.netmeta():
  - new argument 'big.mark' to specify character printed as
    thousands separator, e.g., big.mark = "," will result in
    printing of 1,000 for the number 1000
  - new argument 'scientific.pval' to print p-values in scientific
    notation, e.g., 1.2345e-01 instead of 0.12345

* netmeta():
  - new argument 'backtransf' (see above)
  - new argument 'nchar.trts' to abbreviate treatment names in
    printouts

* netleague():
  - function does not print league table, but only generates it
    (necessary for export of league table)
  - new arguments 'bracket' and 'separator' to define layout of
    confidence intervals (see R function cilayout() from R package
    **meta**)

* netsplit():
  - new argument 'upper' to specify whether lower or upper triangle
    of treatment estimate matrix should be used to build comparisons
  - column with comparisons added to data frames with network, direct,
    and indirect estimates
  - additional new arguments (see help file):
    'reference.group', 'baseline.reference', 'sep.trts', 'quote'

* netmeasures():
  - do not round results to four digits

* print.netmeta(), print.summary.netmeta():
  - new arguments 'backtransf' and 'nchar.trts' (see above)
  - argument 'logscale' removed (replaced by argument 'backtransf')

* Dataset Senn2013:
  - new columns 'treat1.long' and 'treat2.long' with full treatment
    names added

* Help pages:
  - new help pages for forest.netsplit(), print.netleague(), and
    treats()
  - updated help pages for netposet() and netsplit()

### Bug fixes

* netsplit():
  - order of treatments in printouts corresponds to treatment
    comparison, e.g., "A:B" means that treatment "A" was compared
    with treatment "B" (and not the other way around).
    Side note, not sure whether this is a bug or a feature as "A:B"
    noted the design "comparison A and B" so far.

* netposet():
  - function works with a ranking matrix that contains missing
    elements, i.e., rankings that do not include all treatments

* netgraph():
  - areas for multi-arm studies are printed at the correct locations
    if argument 'start.layout' is not equal to "circle" and argument
    'seq' defines a specific treatment order (this bug was
    introduced in netmeta, version 0.7-0)

### Internal changes

* chkmultiarm():
  - warning for zero treatment arm variance instead of an error

* new internal function uppertri() to extract elements from the
  upper triangle of a matrix

* new internal function treats() to abbreviate treatment names

* new internal function nma.additive() for estimation of additive
  network meta-analysis models
  
* new internal function createC() to create C matrix used as input
  to nma.additive()
  
* new internal functions prcombs() and prcomps() for printing of
  netcomb objects

* Internal function p.ci() replaced with formatCI() from R package
  **meta**

* Internal function format.TE() replaced with formatN() from R
  package **meta**


## netmeta, version 0.9-6 (2017-08-09)

### Major changes

* Prediction intervals can be calculated for treatment estimates
  from a network meta-analysis

* In netmeta(), Q statistics for heterogeneity and design
  inconsistency are calculated according to Krahn et al. (2013);
  see help page of decomp.design()

* In printouts and forest plots, the reference treatment can be
  considered as treatment of interest or comparator (default),
  i.e., either comparisons of reference vs other treatments or
  other treatments vs reference are reported

* Tests for heterogeneity and design inconsistency are shown in
  printouts

* A biplot can be generated to show partial ordering of treatment
  rankings for more than two outcomes

* Additional checks implemented for multi-arm studies:
  - negative or zero treatment arm variances
  - duplicate treatment comparisons or incomplete sets of
    treatment comparisons within a study

### User-visible changes

* netmeta():
  - new arguments 'prediction' and 'level.predict' to calculate
    prediction intervals
  - list elements 'Q.heterogeneity' and 'Q.inconsistency' based on
    Krahn et al. (2013)
  - new list elements 'prediction', 'lower.predict',
    'upper.predict', 'df.Q.heterogeneity',
    'pval.Q.heterogeneity','df.Q.inconsistency', and
    'pval.Q.inconsistency'
  - list element 'df' renamed to 'df.Q'
  - stop with an informative error message if (i) any treatment arm
    variance derived from the treatment comparison variances is
    negative or zero, or (ii) in case of duplicate comparisons or an
    incomplete set of treatment comparisons within a study
  - argument 'details.tol.multiarm' renamed to 'details.chkmultiarm'

* netmeta(), forest.netmeta(), print.netmeta(),
  print.summary.netmeta(), and summary.netmeta():
  - new argument 'baseline.reference' to print results for
    comparisons between reference and other treatments, or vice versa

* print.netmeta(), print.summary.netmeta(), and summary.netmeta():
  - new argument 'prediction' to print prediction intervals

* print.summary.netmeta():
  - print information on tests for overall heterogeneity and
    inconsistency

* summary.netmeta():
  - arguments 'level' and 'level.comb' removed from R function
    (i.e., one has to re-run the netmeta() command for confidence
    intervals with other coverage levels)

* plot.netposet():
  - new argument 'plottype' to choose between scatter plot or biplot
  - new arguments to modify layout ('cex.text', 'col.text', pch,
    cex.points, col.points)
  
* decomp.design() and netmeasures():
  - new argument 'warn' to suppress printing of warnings

### Internal changes

* New internal function upgradenetmeta() to add missing list
  elements to older netmeta objects

* R function ci() from R package **meta** added to NAMESPACE

* chkmultiarm():
  - additional checks for (i) negative and zero variances as well as
    (ii) duplicate treatment comparisons or incomplete sets of
    treatment comparisons within a study


## netmeta, version 0.9-5 (2017-05-31)

### Major changes

* New function netleague() to print league table with network
  meta-analysis results

* pairwise():
  - zero events for binary outcomes or incidence rates are handled
    correctly in multi-arm studies by adding an increment to all
    treatment arms (in older versions of netmeta inconsistent
    treatment effects for multi-arm studies were possible as
    increments were considered in individual comparisons instead of
    all comparisons for a multi-arm study)
  - print warning and information on treatment comparisons with
    missing treatment estimate or standard error

* forest.netmeta():
  - reference group can be omitted from forest plot
  - treatments can be sorted by treatment estimate (TE), standard
    error (seTE), number of studies in direct comparison (k), and
    proportion of direct information (prop.direct)

* netmeta():
  - additional checks for correct number of comparisons in multi-arm
    studies and more informative error message for uncorrect number
    of comparisons in multi-arm studies due to missing treatment
    effects or standard errors in single comparisons
  - separator used in comparison names to concatenate treatment
    labels can be specified by user (default: ":")

* In decomp.design(), by default, only print designs contributing to
  design-specific decomposition of within-designs Q statistic

* Input to netdistance() can be either a netmeta object or a matrix

### User-visible changes
    
* forest.netmeta():
  - new argument 'drop.reference.group'
  - argument 'sortvar' can be used in the following ways:
    sortvar = TE, sortvar = -TE, sortvar = seTE, sortvar = -seTE,
    sortvar = k, sortvar = -k, sortvar = prop.direct, sortvar = -prop.direct

* print.decomp.design() and netheat():
  - new argument 'showall' which defaults to FALSE

* print.summary.netmeta():
  - print number of designs
  - print preset between-study variance and corresponding
    information if argument 'tau.preset' is not NULL in netmeta()

* pairwise():
  - in multi-arm studies exclude comparisons with missing sample
    size or standard error from calculation of pooled variance for
    standardized mean difference (sm = "SMD")

* plot.netposet():
  - new default for argument 'arrows', i.e., by default, do not show
    arrows in scatter plot

* print.netsplit():
  - number of studies providing direct evidence printed

* netdistance():
  - argument name changed from 'A' to 'x' in order to reflect that
    input of R function can be either a netmeta object or an
    adjacency matrix
  
* Help pages:
  - examples corrected for dataset dietaryfat
  - do not run all examples in forest.netmeta() as CRAN only allows
    a run time below 10 seconds for examples provided on a help page
  - R code to produce forest plot added to examples in dataset
    Wood2010

### Internal changes

* netmeta():
  - new list element 'd' with number of designs
  - new list element 'B.matrix' with the edge-vertex incidence
    matrix

* summary.netmeta():
  - new list element 'd' with number of designs
  - new list element 'tau.preset'

* netsplit():
  - new list element 'k' with number of studies providing direct
    evidence

* netconnection():
  - argument checks added
  - better code documentation

* Internal function decomp.tau():
  - detach all designs (including protuding edges)

* New internal function createB() to calculate edge-vertex incidence
  matrix

* netmeta(), netconnection(), multiarm(), and chkmultiarm():
  - use internal function createB() instead of dedicated R code

* print.summary.netmeta(), nma.ruecker(), and decomp.tau():
  - use command pchisq(..., lower.tail = FALSE) instead of
    1 - pchisq(...)


## netmeta, version 0.9-4 (2017-04-07)

### Bug fix release

* netsplit() used wrong comparison labels if argument
  'reference.group' was used in netmeta()

* netmeasures() ignores value of argument 'reference.group' in
  netmeta object


## netmeta, version 0.9-3 (2017-03-12)

### Major changes

* Calculate indirect treatment estimates based on direct evidence
  proportion

* Ranking of treatments based on fixed effect model added to
  netrank()

* New function netsplit() to split direct and indirect evidence

* New functions netposet(), print.netposet(), and plot.netposet() to
  calculate, print and plot partial ordering of rankings

* New function hasse() to draw Hasse diagram of partially ordered
  treatment rankings

* netmeta():
  - can be used with R objects created with pairwise()
  - checks for consistency of treatment effects and variances in
    multi-arm studies

* Import ginv() from R package **MASS** (for consistency checks)

* Suggested packages added (for Hasse diagram): **hasseDiagram** and
  **grid**

* Bug fixes:
  - netmeta() calculates correct direct evidence estimates under
    random effects model (list components 'TE.direct.random',
    'seTE.direct.random', ..., 'pval.direct.random'); so far results
    from fixed effect model have been used
  - netmeta() excludes a treatment from list component 'seq' if all
    comparisons containing the respecitve treatment are excluded due to
    missing values in treatment effect or standard error
  - netmeasures() does not result in an error if no or only one
    study with two treatments is available

### User-visible changes
    
* New arguments random and tau.preset in netmeasures()

* New functions netsplit() and print.netsplit()

* Consider ordering of treatments in netrank() which is defined by
  argument seq in netmeta()

* For multi-arm studoes, calculate pooled standard deviation in
  pairwise() if means and standard deviations are provided and
  summary measure is equal to "SMD"

### Internal changes

* netmeta():
  - new list element 'k.direct' with number of studies in
    meta-analyses with direct evidence

* nma.ruecker():
  - bug fix such that estimates from random effects model are used
    for direct treatment estimates if argument 'tau.direct' is
    larger than zero

* nma.krahn():
  - bug fix such that use of function does not result in an error if
    either no or only one study with two treatments is available

* pairwise():
  - data.frame commands use argument stringsAsFactors = FALSE

* chkmultiarm():
  new internal function to check consistency of treatment effects
  and variances in multi-arm studies; calls ginv() from MASS library

* new internal function lowertri() to extract elements from the
  lower triangle of a matrix


## netmeta, version 0.9-2 (2016-11-19)

### Major changes

* R package **rgl** moved from imported to suggested packages as
  - 3-D network plots are not essential for network meta-analysis
  - installation of **netmeta** breaks under macOS if XQuartz is not
    available
    
### User-visible changes
    
* Help page of netgraph() updated (information on **rgl** package)

### Internal changes

* Use chkclass() from **meta** package to check for class membership


## netmeta, version 0.9-1 (2016-10-13)

### Major changes

* Number of studies can be added to network graph

* Distance matrix can be provided directly to generate network graph

* shadowtext() from **TeachingDemos** package by Greg Snow added to
  **netmeta** package

* P-scores can be printed in forest plot
  
### User-visible changes
    
* help page with brief overview of **netmeta** package added

* netgraph():
  - new arguments to add number of studies to network graph
    (number.of.studies, cex.number.of.studies,
    col.number.of.studies, bg.number.of.studies)     
  - plastic look retained for highlighted comparisons
  - new argument D.matrix to provide treatment distances directly

* netmeta():
  - function can be used with a single pairwise comparison without
    resulting in an error

* forest.netmeta():
  - argument sortvar can be equal to Pscore, "Pscore", -Pscore, or
    "-Pscore" to sort treatments according to ranking generated by
    netrank()
  - argument leftcols or rightcols can include "Pscore" to add a
    column with P-Scores to the forest plot
  - new arguments small.values and digits.Pscore for P-Scores

* print.netmeta():
  - use correct layout for network meta-analysis with a single
    pairwise comparison

* decomp.design(), netheat(), netmeasures():
  - print a warning and return NULL for network meta-analysis with a
    single design

* netconnection():
  - print sensible error message if argument treat2 is missing or of
    different length than argument treat 1

* netdistance():
  - print sensible error message if argument A is not a matrix

* Help pages updated:
  decomp.design(), print.decomp.design(),
  netgraph(), netheat(), netmeasures()
    
### Internal changes

* New function:
  - shadowtext() to print number of studies

* nma.ruecker():
  - keep dimension of matrices W and B.matrix for network
    meta-analysis with a single pairwise comparison

* nma.krahn():
  - print a warning and return NULL for network meta-analysis with a
    single design

* decomp.tau(), tau.within():
  - return NULL for network meta-analysis with a single design


## netmeta, version 0.9-0 (2016-04-26)

* New functions:
  - netdistance (calculate distance matrix; replacement for internal
  function nodedist)
  - netconnection (Get connectivity information for network)
  - print.netconnection (corresponding print function)

* Internal function nodedist removed (replaced by netdistance function)

* Import functions from R package **rgl** (for 3-D plots)

* New dataset Woods2010 (use long format in pairwise function)

* Function netmeta:
  - check connectivity of network and stop with informative error
  message if network is not fully connected
  - new list components:
  'Cov.fixed' (variance-covariance matrix for fixed effect model)
  'Cov.random' (variance-covariance matrix for random effects model)

* Function pairwise:
  - extension to long data format (see example on help page)

* Function netmeta:
  - new arguments 'dim', 'eig3', and 'zpos' to generate 3-D network
  plots

* Function stress (used internally):
  - extension to generate 3-D network plots
  - use netdistance function instead of nodedist

* Function nma.ruecker (used internally):
  - use of netmeta function does not result in an error for networks
  without heterogeneity / inconsistency, i.e. networks with zero
  degrees of freedom (e.g. a star-shaped network with only a single
  study for each comparison; simple example: single comparisons A-B,
  A-C, A-D)
  - calculate variance-covariance matrix

* Function print.netrank:
  - print title of meta-analysis (if available)

* Function print.summary.netmeta:
  - print "--" instead of "< 0.0001" in networks without heterogeneity
  / inconsistency
  - print "0" instead of "< 0.0001" if tau-squared is zero
  - print 'p-value' instead of 'p.value'

* Function print.decomp.design:
  - print 'p-value' instead of 'p.value'

* Help page of netmeta function:
  - more details on contrast- and arm-based data format
  - reference to book "Meta-Analysis with R" and Rücker & Schwarzer (2014) added
  - add information that hazard ratio is a possible summary measure
  - change error in description of adjustment in random effects model

* Help page of netgraph function:
  - example for 3-D network plot added

* Help page of netrank function:
  - reference to Rücker & Schwarzer (2015) updated

* Help page of pairwise function:
  - description on use of long data format added
  - more information on additional arguments for meta-analysis
    functions

* New help pages:
  - netconnection, print.netconnection
  - netdistance
  - Wooks2010 dataset


## netmeta, version 0.8-0 (2015-06-26)

* New functions netrank and print.netrank:
  - frequentist method to rank treatments in network

* Function netmeta:
  - print less irritating warning if treatment comparisons are
    resorted (as this is more a note than a warning)

* Function print.netmeta:
  - minor change in printout (old: "Data utilised in network
    meta-analysis ..."; new: "Results ...")

* Help pages:
  - new help page for netrank function
  - reference Rücker & Schwarzer (2015) added in help page of netgraph
    function
  - link to pairwise function added in help page of netmeta function


## netmeta, version 0.7-0 (2015-02-04)

* Version of R package **meta** must be larger or equal 4.0-0

* Title of R package changed

* New function pairwise:
  - transforms data that are given in an arm-based format (e.g. input
    for WinBUGS is of this format) to contrast-based format that can
    be read by function netmeta

* New datasets:
  - dietaryfat (dataset with incidence rates as outcomes)
  - parkinson (continuous outcomes)
  - smokingcessation (binary outcomes)

* Function netmeta:
  - implement a general check for correct number of comparisons for
    multi-arm studies
  - use setseq function to check and set value of argument 'seq'
  - use setref function to check and set value of argument
    'reference.group'
  - use chklevel function from R package **meta** to check levels of
    confidence intervals
  - consider attribute 'sm' from R objects generated with R function
    pairwise
  - function can be used for a pairwise meta-analysis (bug fix in
    nma.ruecker function used internally)

* Function netgraph:
  - check that matrix 'thickness' (if provided) has same row and
    column names as argument 'labels'
  - use setseq function to check and set value of argument 'seq'
  - stop with an error message if argument 'seq' or 'labels' is NULL

* Function netheat:
  - no net heat plot produced if (i) the number of designs is equal or
    smaller than 2 or (ii) no between-design heterogeneity exists
  - unintentional warnings omitted

* Function forest.netmeta:
  - print a warning that the first treatment is used as reference if
    the reference group is unspecified instead of producing an error
  - use setseq function to check and set value of argument 'seq'
  - use setref function to check and set value of argument
    'reference.group'

* Function print.summary.netmeta:
  - print "." instead of "0" or "1" for diagonal elements of treatment
    effect and confidence interval matrices
  - print "." instead of "0" or "1" for reference group (if provided)
  - use setref function to check and set value of argument
    'reference.group'
  - use is.relative.effect function from R package **meta** to check
    if a relative effect measure is used (argument 'sm')

* Function print.netmeta:
  - use setref function to check and set value of argument
    'reference.group'
  - use is.relative.effect function from R package **meta** to check
    if a relative effect measure is used (argument 'sm')

* Function summary.netmeta:
  - use setref function to check and set value of argument
    'reference.group'

* Function decomp.tau and tau.within (used internally):
  - bug fix such that no error is produced in decomp.design and
    netheat function for networks without heterogeneity and
    inconsistency

* Function print.decomp.design:
  - omit printing of information on between-designs Q statistic after
    detaching of single designs if no between-design heterogeneity
    exists
  - use format.tau function from R package **meta** to print "0"
    instead of "< 0.0001" if tau-squared is zero

* New functions (used internally):
  - setseq - check and set argument 'seq' (and argument 'sortvar' in
    forest.meta function)
  - setref - check and set argument 'reference.group'
  - chklist - check for a list

* New help pages for function pairwise and datasets dietaryfat,
  parkinson, and smokingcessation.


## netmeta, version 0.6-0 (2014-07-29)

* Function netgraph:
  - complete rewrite of this function (without changing previous
    default settings substantially)
  - list of major new features:
  - additional layouts beside circular presentation (see argument
    'start.layout')
  - implementation of stress majorization algorithm to optimize layout
    (argument 'iterate')
  - additional methods to determine width of lines connecting
    treatments (argument 'thickness')
  - highlight multi-arm studies (arguments 'multiarm' and
    'col.multiarm')
  - possibility to provide a neighborhood matrix to specify
    neighborhood differently than using the adjacency matrix, for
    example content-based (argument 'N.matrix')
  - possibility to provide x- and y-coordinates for network plot
    (arguments 'xpos' and 'ypos')

* Function netmeta:
  - calculate treatment estimates from all direct pairwise treatment
    comparisons (both fixed effect and random effects model)
  - new list components: 'tau.preset', 'TE.direct.fixed',
   'seTE.direct.fixed', 'lower.direct.fixed', 'upper.direct.fixed',
   'zval.direct.fixed', 'pval.direct.fixed', 'TE.direct.random',
   'seTE.direct.random', 'lower.direct.random', 'upper.direct.random',
   'zval.direct.random', 'pval.direct.random'

* Function nma.ruecker (used internally)
  - changed accordingly to reflect changes in netmeta function

* Function forest.netmeta:
  - new argument sortvar (default: sort treatment effect estimates
    according to list component 'seq' of netmeta object)

* New functions stress and nodedist (used internally)
  - auxiliary functions for netgraph function

* Help pages updated accordingly


## netmeta, version 0.5-0 (2014-06-24)

* Functions nma.krahn, netmeasures, netheat, decomp.design, and
  print.decomp.design:
  - random effects network meta-analysis added

* Function netheat:
  - new argument 'random' 

* Functions nma.krahn, decomp.design, and netheat:
  - new argument 'tau.preset'

* Function decomp.design:
  - correct design-specific decomposition of Q statistic in network
    meta-analysis with multi-arm studies
  - list component 'Q.design' renamed to 'Q.het.design'
  - list component 'Q.detach' renamed to 'Q.inc.detach'
  - list component 'residuals' renamed to 'residuals.inc.detach'
  - new list components: 'Q.inc.random', 'Q.inc.random.preset',
    'Q.inc.design.random.preset',
    'residuals.inc.detach.random.preset', 'tau.preset'

* New functions tau.within and decomp.tau (used internally)

* Help pages updated accordingly


## netmeta, version 0.4-4 (2014-05-27)

* Functions netmeta and nma.ruecker:
  - modified such that the estimated tau-squared in random effects
    model considers multi-arm studies

* Function print.netmeta:
  - information on percentage weight not printed as interpretation is
    difficult

* Dataset Senn2013:
  - use of unpooled standard error for each treatment comparison


## netmeta, version 0.4-3 (2014-04-14)

* Function netmeta:
  - numeric values for arguments 'treat1' and 'treat2' not converted
    to character values (only factors converted to characters)
  - check whether treatments are different (arguments 'treat1' and
    'treat2')

* Function print.summary.netmeta:
  - print random effects estimates according to argument 'seq'

* Function forest.netmeta:
  - sort treatment effect estimates according to argument 'seq'

* Function nma.ruecker (used internally):
  - changed such that all treatment effects are calculated
    irregardless of treatment order (some treatment effects remained
    NA depending on order of treatments)


## netmeta, version 0.4-2 (2014-03-31)

* Function netmeasures:
  - bug fix using correct formula to calculate direct evidence
    proportion (variance instead of standard error)


## netmeta, version 0.4-1 (2014-03-21)

* Function netmeta:
  - argument 'seq' added (see also R function netgraph)

* Function netgraph:
  - new default for argument 'seq'

* Help pages updated accordingly

* Some internal code cleaning to improve readability of R functions


## netmeta, version 0.4-0 (2014-03-07)

* New functions added:
  - netgraph (network graph)
  - netheat (net heat graph)
  - netmeasures (measures for network meta-analysis)
  - decomp.design (design-based decomposition of Cochran's Q)
  - print.decomp.design (corresponding print function)
  - p.ci, format.TE, nma.krahn, nma.ruecker (used internally)

* Function netmeta:
  - Check added whether all pairwise comparisons are provided for
    multi-arm studies

* Help pages added for new functions

* Help page of function netmeta updated


## netmeta, version 0.3-1 (2013-08-01)

* Functions netmeta and summary.netmeta:
  - new list component 'n' (number of treatments)

* Function print.summary.netmeta:
  - modified such that number of treatments is printed
  - modified such that argument 'reference.group' works as expected
    for random effects model


## netmeta, version 0.3-0 (2013-07-24)

 * First version released on CRAN
