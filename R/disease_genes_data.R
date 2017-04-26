#' Disease genes according to GWAS and OMIM
#'
#' Menche and colleagues collected gene-disease associations from the \href{http://www.ebi.ac.uk/gwas/}{GWAS catalog}
#' and the \href{https://www.omim.org/}{Online Mendelian Inheritance in Man} (OMIM) resource to study the topological
#' characteristics of disease modules in the human interactome.
#' \code{disease.genes} contains all genes reported to be associated with human diseases.
#'
#' @docType data
#'
#' @format A character vector with Entrez IDs of genes associated with Human Diseases.
#'
#' @keywords datasets
#'
#' @references Menche J. et al. (2015) Uncovering disease-disease relationships 
#' through the incomplete interactome. Science 347(6224):1257601. DOI:10.1126/science.1257601
#'
#' @source \href{http://science.sciencemag.org/content/suppl/2015/02/18/347.6224.1257601.DC1?_ga=1.180294052.1212156151.1486368973}{Supplementary Material of reference paper.}
#'
#' @examples
#' # Show the first three disease genes
#' disease.genes[1:3]
#' 
"disease.genes"

