#' Yellow-bellied toad microsatellite data
#'
#' Data from Prohl et al 2021
#'
#' @format ## `Toad_microsat_freq`
#' A list of 10 matrices of allele frequencies, each with 47 rows (subpopulations) and a number of columns (alleles):
#' \describe{
#'   \item{row}{Subpopulation sampled}
#'   \item{columm}{microsatellite allele}
#' }
#' @source <Prohl et al 2021>
"Toad_microsat_freq"

#' Microsatellite genotypes of the yellow-bellied toad (Pröhl et al. 2021)
#'
#' Genotypes of individual *Bombina variegata* at 10 microsatellite loci,
#' in the format returned by [microsat_genotypes()].
#'
#' @format ## `toad_genotypes`
#' A data frame with 885 rows (individuals) and 22 columns:
#' \describe{
#'   \item{id}{Individual identifier}
#'   \item{pop}{Subpopulation (factor, 47 levels)}
#'   \item{columns 3--22}{Alleles at the 10 loci, two consecutive columns per locus
#'     (fragment length in base pairs, `NA` = missing). The first column of each
#'     pair is named after the locus.}
#' }
#' @source Pröhl et al. (2021) Conservation Genetics, \doi{10.1007/s10592-021-01350-5}
"toad_genotypes"
