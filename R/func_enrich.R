
#' Perform BP, CC, MF and REACTOME enrichment analyses
#' 
#' Given a list of genes of interest and a reference background, performs a 
#' Fisher's test for the over-representation of Gene Ontology terms and REACTOME
#' pathways in the former.
#' 
#' @param gene.list character; A vector with the genes of interest.
#' @param background character; A vector with the background list of genes.
#' @param id.type character; One of ENTREZID (default), SYMBOL or UNIPROT 
#' accession. This is the ID type of \code{gene.list} and \code{background}.
#' @param benjamini logical; Whether to include Benjamini-Hochberg adjusted 
#' p-values or not.
#' 
#' @return A list with four data frames, one per enrichment analysis:
#' \item{bp}{Contains the Biological Process \code{go.id}s annotating the genes 
#' of interest, together with their \code{term} description, p-values 
#' \code{pval} and Benjamini-Hochberg adjusted p-values \code{bh} if requested.
#' In the latter case, the data frame is sorted by corrected p-value.}
#' \item{cc}{Same as \code{bp} but for Cellular Compartment.}
#' \item{mf}{Same as \code{bp} but for Molecular Function.}
#' \item{reactome}{Same as the rest but the first column is \code{react.id} 
#' instead of \code{go.id}.}
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Use the included lists of disease genes and genes associated with metabolic
#' # disorders as background and genes of interest, respectively
#' analysis <- fun_enrich(gene.list = metabolic, background = disease.genes, 
#'                       id.type = "ENTREZID", benjamini = TRUE)
#' 
#' @export
#' @import org.Hs.eg.db
#' @import GO.db
#' @import reactome.db
#' @importFrom AnnotationDbi select keys
#' 
fun_enrich <- function(gene.list, background, 
                       id.type = "ENTREZID", benjamini = F){
  
  # Get rid of redundant elements in the gene and background lists
  gene.list <- unique(gene.list)
  background <- unique(background)
  
  # The gene list has to be a perfect subset of the background
  if(sum(gene.list %in% background) < length(gene.list) | 
     !(id.type %in% c("ENTREZID", "SYMBOL", "UNIPROT"))){
    stop(paste0("The gene list is not a perfect subset of the background or ", 
                "invalid ID type..."))
  }else if(sum(gene.list %in% keys(org.Hs.eg.db, keytype = id.type)) == 0){
    stop(paste0("None of the elements of the provided gene list is a valid ",
                id.type, "..."))
  }else if(sum(gene.list %in% keys(org.Hs.eg.db, keytype = id.type)) < 
           length(gene.list) | 
           sum(background %in% keys(org.Hs.eg.db, keytype = id.type)) <
           length(background)){
    warning(paste0("Some elements of the provided gene list or background are ",
                   "not valid ", id.type, " IDs.\n", "Enrichments will be ",
                   "perfomed with the valid IDs only."))
  }
  
  # Get all annotations of the given gene list and background and 
  # divide in ontologies
  all.go.genes <- select(org.Hs.eg.db, keys = gene.list, 
                         columns = c("ENTREZID", "GO"), keytype = id.type)
  all.go.bg <- select(org.Hs.eg.db, keys = background, 
                      columns = c("ENTREZID", "GO"), keytype = id.type)
  all.react.genes <- select(reactome.db, keys = unique(all.go.genes$ENTREZID), 
                            columns = c("PATHID"), keytype = "ENTREZID")
  all.react.bg <- select(reactome.db, keys = unique(all.go.bg$ENTREZID), 
                         columns = c("PATHID"), keytype = "ENTREZID")
  
  bp.genes <- table(all.go.genes$GO[all.go.genes$EVIDENCE != "IEA" & 
                                      all.go.genes$ONTOLOGY == "BP" & 
                                      !duplicated(all.go.genes[, 1:3])])
  cc.genes <- table(all.go.genes$GO[all.go.genes$EVIDENCE != "IEA" & 
                                      all.go.genes$ONTOLOGY == "CC" & 
                                      !duplicated(all.go.genes[, 1:3])])
  mf.genes <- table( all.go.genes$GO[all.go.genes$EVIDENCE != "IEA" & 
                                       all.go.genes$ONTOLOGY == "MF" & 
                                       !duplicated(all.go.genes[, 1:3])])
  
  # Inferred from Electronic Annotation (IEA) terms are removed
  bp.bg <- table(all.go.bg$GO[all.go.bg$EVIDENCE != "IEA" & 
                                all.go.bg$ONTOLOGY == "BP" & 
                                !duplicated(all.go.bg[, 1:3])])
  cc.bg <- table(all.go.bg$GO[all.go.bg$EVIDENCE != "IEA" & 
                                all.go.bg$ONTOLOGY == "CC" & 
                                !duplicated(all.go.bg[, 1:3])])
  mf.bg <- table(all.go.bg$GO[all.go.bg$EVIDENCE != "IEA" & 
                                all.go.bg$ONTOLOGY == "MF" & 
                                !duplicated(all.go.bg[, 1:3])])
  
  react.genes <- table(all.react.genes$PATHID)
  react.bg <- table(all.react.bg$PATHID)
  
  # Prepare the output data frames
  bp <- data.frame(go.id = names(bp.genes), 
                   term = select(GO.db, keys = names(bp.genes), 
                                 columns = "TERM", keytype = "GOID")$TERM, 
                   pval = numeric(length = length(names(bp.genes))), 
                   stringsAsFactors = F)
  cc <- data.frame(go.id = names(cc.genes), 
                   term = select(GO.db, keys = names(cc.genes), 
                                 columns = "TERM", keytype = "GOID")$TERM, 
                   pval = numeric(length = length(names(cc.genes))), 
                   stringsAsFactors = F)
  mf <- data.frame(go.id = names(mf.genes), 
                   term = select(GO.db, keys = names(mf.genes), 
                                 columns = "TERM", keytype = "GOID")$TERM, 
                   pval = numeric(length = length(names(mf.genes))), 
                   stringsAsFactors = F)
  react <- data.frame(react.id = names(react.genes), 
                      pathway = select(reactome.db, keys = names(react.genes), 
                                       columns = "PATHNAME", 
                                       keytype = "PATHID")$PATHNAME, 
                      pval = numeric(length = length(names(react.genes))), 
                      stringsAsFactors = F)
  
  bp$pval <- do_fisher_tests(bp.genes, bp.bg)
  cc$pval <- do_fisher_tests(cc.genes, cc.bg)
  mf$pval <- do_fisher_tests(mf.genes, mf.bg)
  react$pval <- do_fisher_tests(react.genes, react.bg)
  
  # Correct p-values if requested and sort by corrected p-value
  if(benjamini){
    bp$bh <- stats::p.adjust(bp$pval, method = "BH")
    cc$bh <- stats::p.adjust(cc$pval, method = "BH")
    mf$bh <- stats::p.adjust(mf$pval, method = "BH")
    react$bh <- stats::p.adjust(react$pval, method = "BH")
    
    bp <- bp[order(bp$bh), ]
    cc <- cc[order(cc$bh), ]
    mf <- mf[order(mf$bh), ]
    react <- react[order(react$bh), ]
  }else{
    # Sort by p-value
    bp <- bp[order(bp$pval), ]
    cc <- cc[order(cc$pval), ]
    mf <- mf[order(mf$pval), ]
    react <- react[order(react$pval), ]
  }
  
  return(list(bp = bp, cc = cc, mf = mf, reactome = react))
  
}

#' Perform Fisher's tests
#' 
#' Carry out the required over-representation tests for each of the terms 
#' annotating the genes of interest.
#' 
#' @param gene.counts table; The term counts for the genes of interest.
#' @param bg.counts table; The term counts for the background.
#' 
#' @return A vector with the p-values for each term.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Generate a hypotetical scenario of annotations
#' bg.annotations <- sample(x = base::letters, size = 100, replace = TRUE)
#' gene.annotations <- sample(bg.annotations, size = 30)
#' 
#' # Compute the frequence of terms
#' bg.counts <- table(bg.annotations)
#' gene.counts <- table(gene.annotations)
#' 
#' #Perform Fisher's tests of over-representation for each term
#' pvals <- do_fisher_tests(gene.counts = gene.counts, bg.counts = bg.counts)
#' 
#' @export
#' 
do_fisher_tests <- function(gene.counts, bg.counts){
  
  pvals <- numeric(length = length(gene.counts))
  
  contingency <- matrix(0, nrow = 2, ncol = 2)
  j <- 1
  for(i in names(gene.counts)){
    # Construction of the contingency matrix
    contingency[1, 1] <- gene.counts[i] # k
    contingency[1, 2] <- bg.counts[i] - gene.counts[i] # K - k
    contingency[2, 1] <- sum(gene.counts) - gene.counts[i] # n - k
    contingency[2, 2] <- sum(bg.counts) + gene.counts[i] - 
      sum(gene.counts) - bg.counts[i] # N + k - n - K
    pvals[j] <- stats::fisher.test(contingency, alternative = "greater")$p.value
    j <- j+1
  }
  return(pvals)
}

#' Plot the results of the BP, CC, MF and/or REACTOME enrichment analyses
#' 
#' Given the output of \code{fun_enrich}, generates a bar plot of the -log10 of
#' the p-values for one or all analysed aspects. 
#' 
#' @param enr listr; The list of data frames resulting from \code{fun_enrich}.
#' @param aspect character; One of ALL, BP, CC, MF or REACTOME.
#' @param benjamini logical; If TRUE, plot corrected p-values.
#' @param top integer; The number of top enriched terms to consider.
#' @param char_per_line integer; Controls the number of characters per line for
#' the term names in the plots.
#' 
#' @return A ggplot object.
#' 
#' @author Gregorio Alanis-Lobato \email{galanisl@uni-mainz.de}
#' 
#' @examples
#' # Use the included lists of disease genes and genes associated with metabolic
#' # disorders as background and genes of interest, respectively
#' analysis <- fun_enrich(gene.list = metabolic, background = disease.genes, 
#'                       id.type = "ENTREZID", benjamini = TRUE)
#' # Plot the 
#' vis <- plot_fun_enrich(enr = analysis, aspect = "ALL", benjamini = TRUE, 
#'                        top = 5, char_per_line = 80)
#' 
#' @export
#' @import ggplot2
#' @import stringr
#' 
plot_fun_enrich <- function(enr, aspect = "ALL", benjamini = FALSE, 
                            top = 5, char_per_line = 80){
  if(benjamini & is.null(enr$bp$bh)){
    stop(paste0("Correction of p-values was not requested for the provided ", 
                "enrichment analysis..."))
  }else{
    if(aspect == "ALL"){
      colnames(enr$reactome)[1:2] <- c("go.id", "term")
      df <- rbind(enr$bp[top:1,], enr$cc[top:1,], enr$mf[top:1,], 
                  enr$reactome[top:1,])
      df$aspect <- factor(rep(c("BP", "CC", "MF", "REACTOME"), each = top),
                          levels = c("BP", "CC", "MF", "REACTOME"),
                          ordered = T)
      df <- df[order(df$aspect, decreasing = T), ]
    }else if(aspect == "BP"){
      df <- enr$bp[top:1,]
      df$aspect <- "BP"
    }else if(aspect == "CC"){
      df <- enr$cc[top:1,]
      df$aspect <- "CC"
    }else if(aspect == "MF"){
      df <- enr$mf[top:1,]
      df$aspect <- "MF"
    }else if(aspect == "REACTOME"){
      colnames(enr$reactome)[2] <- c("term")
      df <- enr$reactome[top:1,]
      df$aspect <- "REACTOME"
    }else{
      stop(paste0("The requested aspect is not valid. Choose one from ALL, BP, ", 
                  "CC, MF or REACTOME."))
    }
    if(benjamini){
      df$prob <- -log10(df$bh)
    }else{
      df$prob <- -log10(df$pval)
    }
    df$term <- str_wrap(df$term, width = char_per_line)
    df$term <- factor(df$term, levels = df$term)
    p <- ggplot(df, aes_(~term, ~prob, fill = ~aspect)) + geom_col() + 
      coord_flip() + labs(x = "", y = expression(-log[10](p-value))) + 
      theme_bw() + theme(legend.title = element_blank(), 
                         legend.background = element_blank(), 
                         legend.position="top")
    return(p)
  }
}