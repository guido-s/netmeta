\name{Linde2015}
\docType{data}
\alias{Linde2015}
\title{Network meta-analysis of treatments for depression}

\description{
Network meta-analysis of nine classes of antidepressants including
placebo for the primary care setting; partly shown in Linde et al.
(2015), supplementary Table 2.
}

\usage{data(Linde2015)}

\format{A data frame with the following columns:
  \describe{
  \item{id}{Study ID}
  \item{author}{First author}
  \item{year}{Publication year}
  \item{treatment1}{First treatment}
  \item{treatment2}{Second treatment}
  \item{treatment3}{Third treatment}
  \item{n1}{Number of patients receiving first treatment}
  \item{resp1}{Number of early responder (treatment 1)}
  \item{remi1}{Number of early remissions (treatment 1)}
  \item{loss1}{Number of patients loss to follow-up (treatment 1)}  
  \item{loss.ae1}{Number of patients loss to follow-up due to adverse
    events (treatment 1)}
  \item{ae1}{Number of patients with adverse events (treatment 1)}
  \item{n2}{Number of patients receiving second treatment}
  \item{resp2}{Number of early responder (treatment 2)}
  \item{remi2}{Number of early remissions (treatment 2)}
  \item{loss2}{Number of patients loss to follow-up (treatment 2)}  
  \item{loss.ae2}{Number of patients loss to follow-up due to adverse
    events (treatment 2)}
  \item{ae2}{Number of patients with adverse events (treatment 2)}
  \item{n3}{Number of patients receiving third treatment}
  \item{resp3}{Number of early responder (treatment 3)}
  \item{remi3}{Number of early remissions (treatment 3)}
  \item{loss3}{Number of patients loss to follow-up (treatment 3)}  
  \item{loss.ae3}{Number of patients loss to follow-up due to adverse
    events (treatment 3)}
  \item{ae3}{Number of patients with adverse events (treatment 3)}
  }
}

\source{
  Linde K, Kriston L, Rücker G, et al. (2015).
  Efficacy and acceptability of pharmacological treatments for
  depressive disorders in primary care: Systematic review and network
  meta-analysis.
  \emph{Annals of Family Medicine}
  \bold{13}, 69--79.
}

\seealso{\code{\link{pairwise}}, \code{\link{metabin}}, \code{\link{netmeta}}, \code{\link{netposet}}}

\examples{
data(Linde2015)

# Transform data from arm-based format to contrast-based format
# Outcome: early response
p1 <- pairwise(list(treatment1, treatment2, treatment3),
               event = list(resp1, resp2, resp3),
	       n = list(n1, n2, n3),
               studlab = id, data = Linde2015, sm = "OR")

# Define order of treatments
trts <- c("TCA", "SSRI", "SNRI", "NRI",
          "Low-dose SARI", "NaSSa", "rMAO-A", "Hypericum",
          "Placebo")

# Conduct network meta-analysis
net1 <- netmeta(p1, comb.fixed = FALSE, comb.random = TRUE,
                reference = "Placebo",
		seq = trts)
print(summary(net1), digits = 2)
}

\keyword{datasets}
