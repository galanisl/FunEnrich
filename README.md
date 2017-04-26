FunEnrich
================
Gregorio Alanis-Lobato

BP, CC, MF (i.e. Gene Ontology) and REACTOME enrichment analyses for a list of genes of interest, given a list of background genes.

The package supports ENTREZIDs (default), GENE SYMBOLs and UNIPROT accessions.

Installation
============

1.  Install the `devtools` package from CRAN if you haven't done so:

``` r
install.packages("devtools")
```

1.  Load the `devtools` package:

``` r
library("devtools")
```

1.  Install `FunEnrich` using the `install_github` function:

``` r
install_github("galanisl/FunEnrich")
```

Usage
=====

To start using `FunEnrich`, load the package:

``` r
library("FunEnrich")
```

Let's now use the gene-disease associations included in the package. We will use `disease.genes` as background and genes associated with `metabolic` disorders as genes of interest (type `?metabolic` and `?disease.genes` in R for more information about these datasets):

``` r
analysis <- fun_enrich(gene.list = metabolic, background = disease.genes, 
                       id.type = "ENTREZID", benjamini = FALSE)
```

We can now explore, for example, the Molecular Functions enriched in genes that are specifically associated with metabolic disorders...

``` r
head(analysis$mf)
```

    ##          go.id                                                 term
    ## 334 GO:0005254                            chloride channel activity
    ## 577 GO:0030170                          pyridoxal phosphate binding
    ## 297 GO:0005159          insulin-like growth factor receptor binding
    ## 111 GO:0004029                aldehyde dehydrogenase (NAD) activity
    ## 260 GO:0005007 fibroblast growth factor-activated receptor activity
    ## 314 GO:0005219 ryanodine-sensitive calcium-release channel activity
    ##             pval
    ## 334 7.647321e-05
    ## 577 6.699617e-04
    ## 297 7.688331e-04
    ## 111 2.142760e-03
    ## 260 2.142760e-03
    ## 314 2.142760e-03

... as well as the enriched REACTOME pathways:

``` r
head(analysis$reactome)
```

    ##     react.id
    ## 91   1430728
    ## 885    71387
    ## 670   446193
    ## 672   446203
    ## 180  1660662
    ## 873    70326
    ##                                                                                                                                 pathway
    ## 91                                                                                                             Homo sapiens: Metabolism
    ## 885                                                                                           Homo sapiens: Metabolism of carbohydrates
    ## 670 Homo sapiens: Biosynthesis of the N-glycan precursor (dolichol lipid-linked oligosaccharide, LLO) and transfer to a nascent protein
    ## 672                                                                                     Homo sapiens: Asparagine N-linked glycosylation
    ## 180                                                                                          Homo sapiens: Glycosphingolipid metabolism
    ## 873                                                                                                    Homo sapiens: Glucose metabolism
    ##             pval
    ## 91  4.267708e-13
    ## 885 2.154152e-09
    ## 670 6.281716e-08
    ## 672 8.262870e-08
    ## 180 1.894112e-06
    ## 873 2.034622e-06
