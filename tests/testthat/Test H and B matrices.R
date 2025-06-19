
runTomsHandBTest <- function(data1){
library(corpcor)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


#----------------------------------------------
# Binary, fixed effects
#----------------------------------------------
data1 <- read.csv("NMA_data_binary_FE.csv")

pairwise1 <- meta::pairwise(
  treat = T,
  event = R,
  n = N,
  studlab = Study,
  data = data1,
  sm = "OR"
)

netmeta1 <- netmeta::netmeta(
  TE = TE,
  seTE = seTE,
  treat1 = treat1,
  treat2 = treat2,
  studlab = studlab,
  common = TRUE,
  data = pairwise1,
  sm = "OR"
)

contributions1 <- netmeta::netcontrib(
  x = netmeta1,
  method = "shortestpath",
  random = FALSE
)

#The current version of CINeMA uses the old hat matrix
# contributions1 <- netmeta::netcontrib(
#   x = netmeta1,
#   method = "shortestpath",
#   random = FALSE,
#   hatmatrix.F1000 = TRUE
# )


#The B matrix for 3-arm studies described in Rucker, Schwarzer
B_matrix_3_arms <- matrix(
  c(1, -1, 0,
    1, 0, -1,
    0, 1, -1),
  byrow = TRUE,
  nrow = 3
)

#Unadjusted weights
weights <- 1 / pairwise1$seTE^2
names(weights) <- paste0(pairwise1$studlab, "_", pairwise1$treat1, "_", pairwise1$treat2)

#Variance matrix for study ABC1
index_ABC1 <- grep("ABC1", pairwise1$studlab)
comparison_variances_ABC1 <- pairwise1$seTE[index_ABC1]^2
names(comparison_variances_ABC1) <- paste0(pairwise1$treat1[index_ABC1], "_", pairwise1$treat2[index_ABC1])
treatments_ABC1 <- data1$T[data1$Study == "StudyABC1"]
variance_matrix_ABC1 <- matrix(0, nrow = 3, ncol = 3)
rownames(variance_matrix_ABC1) <- treatments_ABC1
colnames(variance_matrix_ABC1) <- treatments_ABC1
variance_matrix_ABC1["A", "B"] <- comparison_variances_ABC1["A_B"]
variance_matrix_ABC1["B", "A"] <- variance_matrix_ABC1["A", "B"]
variance_matrix_ABC1["A", "C"] <- comparison_variances_ABC1["A_C"]
variance_matrix_ABC1["C", "A"] <- variance_matrix_ABC1["A", "C"]
variance_matrix_ABC1["B", "C"] <- comparison_variances_ABC1["B_C"]
variance_matrix_ABC1["C", "B"] <- variance_matrix_ABC1["B", "C"]
n_arms_ABC1 <- length(index_ABC1)
#Inverse Laplacian for StudyABC1
pseudo_inverse_laplacian_ABC1 <- -t(B_matrix_3_arms) %*% B_matrix_3_arms %*% variance_matrix_ABC1 %*% t(B_matrix_3_arms) %*% B_matrix_3_arms / (2 * n_arms_ABC1 ^ 2)
#Laplacian for StudyABC1
laplacian_ABC1 <- corpcor::pseudoinverse(pseudo_inverse_laplacian_ABC1)
rownames(laplacian_ABC1) <- treatments_ABC1
colnames(laplacian_ABC1) <- treatments_ABC1

#Variance matrix for study ABC2
index_ABC2 <- grep("ABC2", pairwise1$studlab)
comparison_variances_ABC2 <- pairwise1$seTE[index_ABC2]^2
names(comparison_variances_ABC2) <- paste0(pairwise1$treat1[index_ABC2], "_", pairwise1$treat2[index_ABC2])
treatments_ABC2 <- data1$T[data1$Study == "StudyABC2"]
variance_matrix_ABC2 <- matrix(0, nrow = 3, ncol = 3)
rownames(variance_matrix_ABC2) <- treatments_ABC2
colnames(variance_matrix_ABC2) <- treatments_ABC2
variance_matrix_ABC2["A", "B"] <- comparison_variances_ABC2["A_B"]
variance_matrix_ABC2["B", "A"] <- variance_matrix_ABC2["A", "B"]
variance_matrix_ABC2["A", "C"] <- comparison_variances_ABC2["A_C"]
variance_matrix_ABC2["C", "A"] <- variance_matrix_ABC2["A", "C"]
variance_matrix_ABC2["B", "C"] <- comparison_variances_ABC2["B_C"]
variance_matrix_ABC2["C", "B"] <- variance_matrix_ABC2["B", "C"]
n_arms_ABC2 <- length(index_ABC2)
#Inverse Laplacian for StudyABC2
pseudo_inverse_laplacian_ABC2 <- -t(B_matrix_3_arms) %*% B_matrix_3_arms %*% variance_matrix_ABC2 %*% t(B_matrix_3_arms) %*% B_matrix_3_arms / (2 * n_arms_ABC2 ^ 2)
#Laplacian for StudyABC2
laplacian_ABC2 <- corpcor::pseudoinverse(pseudo_inverse_laplacian_ABC2)
rownames(laplacian_ABC2) <- treatments_ABC2
colnames(laplacian_ABC2) <- treatments_ABC2

#Updated weights taking multi-arm studies into account
updated_weights <- weights
updated_weights["StudyABC1_A_B"] <- -laplacian_ABC1["A", "B"]
updated_weights["StudyABC1_A_C"] <- -laplacian_ABC1["A", "C"]
updated_weights["StudyABC1_B_C"] <- -laplacian_ABC1["B", "C"]
updated_weights["StudyABC2_A_B"] <- -laplacian_ABC2["A", "B"]
updated_weights["StudyABC2_A_C"] <- -laplacian_ABC2["A", "C"]
updated_weights["StudyABC2_B_C"] <- -laplacian_ABC2["B", "C"]

#Contributions to the AB comparison (top row of the contribution matrix)
study_contributions_to_AB <- vector(length = length(unique(data1$Study)))
names(study_contributions_to_AB) <- unique(data1$Study)

#Contributions from direct AB comparison
index_AB <- grep("_A_B", names(updated_weights), value = TRUE)
sum_AB_weights <- sum(updated_weights[index_AB])
study_contributions_to_AB["StudyAB1"] <- contributions1$common["A:B", "A:B"] * updated_weights["StudyAB1_A_B"] / sum_AB_weights
study_contributions_to_AB["StudyAB2"] <- contributions1$common["A:B", "A:B"] * updated_weights["StudyAB2_A_B"] / sum_AB_weights
study_contributions_to_AB["StudyABC1"] <- contributions1$common["A:B", "A:B"] * updated_weights["StudyABC1_A_B"] / sum_AB_weights
study_contributions_to_AB["StudyABC2"] <- contributions1$common["A:B", "A:B"] * updated_weights["StudyABC2_A_B"] / sum_AB_weights

#Contributions from direct AC comparison
index_AC <- grep("_A_C", names(updated_weights), value = TRUE)
sum_AC_weights <- sum(updated_weights[index_AC])
study_contributions_to_AB["StudyABC1"] <- study_contributions_to_AB["StudyABC1"] + contributions1$common["A:B", "A:C"] * updated_weights["StudyABC1_A_C"] / sum_AC_weights
study_contributions_to_AB["StudyABC2"] <- study_contributions_to_AB["StudyABC2"] + contributions1$common["A:B", "A:C"] * updated_weights["StudyABC2_A_C"] / sum_AC_weights
study_contributions_to_AB["StudyAC1"] <- contributions1$common["A:B", "A:C"] * updated_weights["StudyAC1_A_C"] / sum_AC_weights
study_contributions_to_AB["StudyAC2"] <- contributions1$common["A:B", "A:C"] * updated_weights["StudyAC2_A_C"] / sum_AC_weights

#Contributions from direct AD comparison
index_AD <- grep("_A_D", names(updated_weights), value = TRUE)
sum_AD_weights <- sum(updated_weights[index_AD])
study_contributions_to_AB["StudyAD1"] <- contributions1$common["A:B", "A:D"] * updated_weights["StudyAD1_A_D"] / sum_AD_weights
study_contributions_to_AB["StudyAD2"] <- contributions1$common["A:B", "A:D"] * updated_weights["StudyAD2_A_D"] / sum_AD_weights

#Contributions from direct BC comparison
index_BC <- grep("_B_C", names(updated_weights), value = TRUE)
sum_BC_weights <- sum(updated_weights[index_BC])
study_contributions_to_AB["StudyABC1"] <- study_contributions_to_AB["StudyABC1"] + contributions1$common["A:B", "B:C"] * updated_weights["StudyABC1_B_C"] / sum_BC_weights
study_contributions_to_AB["StudyABC2"] <- study_contributions_to_AB["StudyABC2"] +contributions1$common["A:B", "B:C"] * updated_weights["StudyABC2_B_C"] / sum_BC_weights
study_contributions_to_AB["StudyBC1"] <- contributions1$common["A:B", "B:C"] * updated_weights["StudyBC1_B_C"] / sum_BC_weights
study_contributions_to_AB["StudyBC2"] <- contributions1$common["A:B", "B:C"] * updated_weights["StudyBC2_B_C"] / sum_BC_weights

#Final contributions to AB
study_contributions_to_AB


  res <- list(contrs = study_contributions_to_AB,
              weights = updated_weights
              )
  return(res)
}