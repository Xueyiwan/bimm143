class15
================

The data for for hands-on session comes from GEO entry: GSE37704, which
is associated with the following publication:

> Trapnell C, Hendrickson DG, Sauvageau M, Goff L et al. “Differential
> analysis of gene regulation at transcript resolution with RNA-seq”.
> Nat Biotechnol 2013 Jan;31(1):46-53. PMID: 23222703

## Import count and metadata

``` r
metaFile <- "GSE37704_metadata.csv"
countFile <- "GSE37704_featurecounts.csv"

# Import metadata and take a peak
colData = read.csv(metaFile, row.names=1)
head(colData)
```

    ##               condition
    ## SRR493366 control_sirna
    ## SRR493367 control_sirna
    ## SRR493368 control_sirna
    ## SRR493369      hoxa1_kd
    ## SRR493370      hoxa1_kd
    ## SRR493371      hoxa1_kd

``` r
# Import countdata
countData = read.csv(countFile, row.names=1)
head(countData)
```

    ##                 length SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ## ENSG00000186092    918         0         0         0         0         0
    ## ENSG00000279928    718         0         0         0         0         0
    ## ENSG00000279457   1982        23        28        29        29        28
    ## ENSG00000278566    939         0         0         0         0         0
    ## ENSG00000273547    939         0         0         0         0         0
    ## ENSG00000187634   3214       124       123       205       207       212
    ##                 SRR493371
    ## ENSG00000186092         0
    ## ENSG00000279928         0
    ## ENSG00000279457        46
    ## ENSG00000278566         0
    ## ENSG00000273547         0
    ## ENSG00000187634       258

remember that we need the countData and colData files to match up so we
will need to remove that odd first column in countData namely
countData$length.

``` r
# Note we need to remove the odd first $length col
countData <- as.matrix(countData[,-1])
head(countData)
```

    ##                 SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ## ENSG00000186092         0         0         0         0         0
    ## ENSG00000279928         0         0         0         0         0
    ## ENSG00000279457        23        28        29        29        28
    ## ENSG00000278566         0         0         0         0         0
    ## ENSG00000273547         0         0         0         0         0
    ## ENSG00000187634       124       123       205       207       212
    ##                 SRR493371
    ## ENSG00000186092         0
    ## ENSG00000279928         0
    ## ENSG00000279457        46
    ## ENSG00000278566         0
    ## ENSG00000273547         0
    ## ENSG00000187634       258

Double check that the colnmes in countData match the id values in the
metadata file.

``` r
colnames(countData)
```

    ## [1] "SRR493366" "SRR493367" "SRR493368" "SRR493369" "SRR493370" "SRR493371"

``` r
rownames(colData)
```

    ## [1] "SRR493366" "SRR493367" "SRR493368" "SRR493369" "SRR493370" "SRR493371"

``` r
all( colnames(countData) == rownames(colData) )
```

    ## [1] TRUE

``` r
# Test how all() function works
all( c(T, T, T) )
```

    ## [1] TRUE

``` r
all( c(T, F, T) )
```

    ## [1] FALSE

## Remove zero count genes

We want to remove genes that have 0 count values in all experiments
(i.e. rows that have 0 across all cols)

``` r
# Filter out
countData = as.matrix(countData[rowSums(countData) != 0, ])
head(countData)
```

    ##                 SRR493366 SRR493367 SRR493368 SRR493369 SRR493370
    ## ENSG00000279457        23        28        29        29        28
    ## ENSG00000187634       124       123       205       207       212
    ## ENSG00000188976      1637      1831      2383      1226      1326
    ## ENSG00000187961       120       153       180       236       255
    ## ENSG00000187583        24        48        65        44        48
    ## ENSG00000187642         4         9        16        14        16
    ##                 SRR493371
    ## ENSG00000279457        46
    ## ENSG00000187634       258
    ## ENSG00000188976      1504
    ## ENSG00000187961       357
    ## ENSG00000187583        64
    ## ENSG00000187642        16

# DESeq analysis

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
    ##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
    ##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    ##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
    ##     setdiff, sort, table, tapply, union, unique, unsplit, which,
    ##     which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## Loading required package: BiocParallel

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

``` r
# Setup the object with our data in the way DESeq wants it
dds = DESeqDataSetFromMatrix(countData=countData,
                             colData=colData,
                             design=~condition)
# Run the analysis
dds = DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

Get our results

``` r
res = results(dds)

res
```

    ## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    ## Wald test p-value: condition hoxa1 kd vs control sirna 
    ## DataFrame with 15975 rows and 6 columns
    ##                         baseMean     log2FoldChange              lfcSE
    ##                        <numeric>          <numeric>          <numeric>
    ## ENSG00000279457 29.9135794276176   0.17925708367269  0.324821565250145
    ## ENSG00000187634 183.229649921658  0.426457118403306  0.140265820376892
    ## ENSG00000188976 1651.18807619944 -0.692720464846366 0.0548465415913946
    ## ENSG00000187961 209.637938486147  0.729755610585225  0.131859899969345
    ## ENSG00000187583 47.2551232589398 0.0405765278756312  0.271892808601774
    ## ...                          ...                ...                ...
    ## ENSG00000273748 35.3026523877463  0.674387102558605  0.303666187454139
    ## ENSG00000278817 2.42302393023632 -0.388988266500022   1.13039377720313
    ## ENSG00000278384 1.10179649846993  0.332990658240634    1.6602614216556
    ## ENSG00000276345 73.6449563127136 -0.356180759105173   0.20771565839825
    ## ENSG00000271254 181.595902546813 -0.609666545167282  0.141320482803512
    ##                               stat               pvalue
    ##                          <numeric>            <numeric>
    ## ENSG00000279457  0.551863246932648    0.581042050747032
    ## ENSG00000187634   3.04034951107421  0.00236303749730996
    ## ENSG00000188976  -12.6301576133481 1.43989540156582e-36
    ## ENSG00000187961   5.53432552849563  3.1242824807768e-08
    ## ENSG00000187583  0.149237223611387    0.881366448669148
    ## ...                            ...                  ...
    ## ENSG00000273748   2.22081723425481    0.026363342804782
    ## ENSG00000278817 -0.344117487502873    0.730757932009185
    ## ENSG00000278384  0.200565196478864    0.841038574220432
    ## ENSG00000276345  -1.71475160732598   0.0863907773559434
    ## ENSG00000271254  -4.31407063627813 1.60275699407084e-05
    ##                                 padj
    ##                            <numeric>
    ## ENSG00000279457      0.6865547778329
    ## ENSG00000187634  0.00515718149494361
    ## ENSG00000188976  1.7654890539332e-35
    ## ENSG00000187961 1.13412993107607e-07
    ## ENSG00000187583    0.919030615571381
    ## ...                              ...
    ## ENSG00000273748   0.0479091179108357
    ## ENSG00000278817    0.809772069001658
    ## ENSG00000278384    0.892653531513564
    ## ENSG00000276345    0.139761501281218
    ## ENSG00000271254 4.53647639305091e-05

``` r
plot( res$log2FoldChange, -log(res$padj) )
```

![](class15_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# Make a color vector for all genes
mycols <- rep("grey", length(res$padj))

# Color red the genes with absolute fold change above 2
mycols[ abs(res$log2FoldChange) > 2] <- "red"

# Color blue those with adjusted p-value less than 0.01 and absolute fold change more than 2
inds <- (res$padj < 0.01) & (abs(res$log2FoldChange) > 2)
mycols[ inds ] <- "blue"

plot( res$log2FoldChange, -log(res$padj), col=mycols, xlab="Log2(FoldChange)", ylab="-Log(P-value)" )
```

![](class15_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# Adding gene zymbols and entrez ids

``` r
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
```

``` r
library("AnnotationDbi")
library("org.Hs.eg.db")
```

    ## 

``` r
# We can translate betwen al the following databse ID systems
columns(org.Hs.eg.db)
```

    ##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
    ##  [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
    ##  [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
    ## [13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
    ## [17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
    ## [21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
    ## [25] "UNIGENE"      "UNIPROT"

``` r
res$symbol = mapIds(org.Hs.eg.db,
                    keys=row.names(countData), # where are my IDs
                    keytype="ENSEMBL",    # what format are my IDs
                    column="SYMBOL",      # the new format I want
                    multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
res
```

    ## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    ## Wald test p-value: condition hoxa1 kd vs control sirna 
    ## DataFrame with 15975 rows and 7 columns
    ##                         baseMean     log2FoldChange              lfcSE
    ##                        <numeric>          <numeric>          <numeric>
    ## ENSG00000279457 29.9135794276176   0.17925708367269  0.324821565250145
    ## ENSG00000187634 183.229649921658  0.426457118403306  0.140265820376892
    ## ENSG00000188976 1651.18807619944 -0.692720464846366 0.0548465415913946
    ## ENSG00000187961 209.637938486147  0.729755610585225  0.131859899969345
    ## ENSG00000187583 47.2551232589398 0.0405765278756312  0.271892808601774
    ## ...                          ...                ...                ...
    ## ENSG00000273748 35.3026523877463  0.674387102558605  0.303666187454139
    ## ENSG00000278817 2.42302393023632 -0.388988266500022   1.13039377720313
    ## ENSG00000278384 1.10179649846993  0.332990658240634    1.6602614216556
    ## ENSG00000276345 73.6449563127136 -0.356180759105173   0.20771565839825
    ## ENSG00000271254 181.595902546813 -0.609666545167282  0.141320482803512
    ##                               stat               pvalue
    ##                          <numeric>            <numeric>
    ## ENSG00000279457  0.551863246932648    0.581042050747032
    ## ENSG00000187634   3.04034951107421  0.00236303749730996
    ## ENSG00000188976  -12.6301576133481 1.43989540156582e-36
    ## ENSG00000187961   5.53432552849563  3.1242824807768e-08
    ## ENSG00000187583  0.149237223611387    0.881366448669148
    ## ...                            ...                  ...
    ## ENSG00000273748   2.22081723425481    0.026363342804782
    ## ENSG00000278817 -0.344117487502873    0.730757932009185
    ## ENSG00000278384  0.200565196478864    0.841038574220432
    ## ENSG00000276345  -1.71475160732598   0.0863907773559434
    ## ENSG00000271254  -4.31407063627813 1.60275699407084e-05
    ##                                 padj       symbol
    ##                            <numeric>  <character>
    ## ENSG00000279457      0.6865547778329           NA
    ## ENSG00000187634  0.00515718149494361       SAMD11
    ## ENSG00000188976  1.7654890539332e-35        NOC2L
    ## ENSG00000187961 1.13412993107607e-07       KLHL17
    ## ENSG00000187583    0.919030615571381      PLEKHN1
    ## ...                              ...          ...
    ## ENSG00000273748   0.0479091179108357           NA
    ## ENSG00000278817    0.809772069001658 LOC102724770
    ## ENSG00000278384    0.892653531513564           NA
    ## ENSG00000276345    0.139761501281218 LOC107987373
    ## ENSG00000271254 4.53647639305091e-05 LOC102724250

``` r
res$entrez = mapIds(org.Hs.eg.db,
                    keys=row.names(countData),
                    keytype="ENSEMBL",
                    column="ENTREZID",
                    multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
res
```

    ## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    ## Wald test p-value: condition hoxa1 kd vs control sirna 
    ## DataFrame with 15975 rows and 8 columns
    ##                         baseMean     log2FoldChange              lfcSE
    ##                        <numeric>          <numeric>          <numeric>
    ## ENSG00000279457 29.9135794276176   0.17925708367269  0.324821565250145
    ## ENSG00000187634 183.229649921658  0.426457118403306  0.140265820376892
    ## ENSG00000188976 1651.18807619944 -0.692720464846366 0.0548465415913946
    ## ENSG00000187961 209.637938486147  0.729755610585225  0.131859899969345
    ## ENSG00000187583 47.2551232589398 0.0405765278756312  0.271892808601774
    ## ...                          ...                ...                ...
    ## ENSG00000273748 35.3026523877463  0.674387102558605  0.303666187454139
    ## ENSG00000278817 2.42302393023632 -0.388988266500022   1.13039377720313
    ## ENSG00000278384 1.10179649846993  0.332990658240634    1.6602614216556
    ## ENSG00000276345 73.6449563127136 -0.356180759105173   0.20771565839825
    ## ENSG00000271254 181.595902546813 -0.609666545167282  0.141320482803512
    ##                               stat               pvalue
    ##                          <numeric>            <numeric>
    ## ENSG00000279457  0.551863246932648    0.581042050747032
    ## ENSG00000187634   3.04034951107421  0.00236303749730996
    ## ENSG00000188976  -12.6301576133481 1.43989540156582e-36
    ## ENSG00000187961   5.53432552849563  3.1242824807768e-08
    ## ENSG00000187583  0.149237223611387    0.881366448669148
    ## ...                            ...                  ...
    ## ENSG00000273748   2.22081723425481    0.026363342804782
    ## ENSG00000278817 -0.344117487502873    0.730757932009185
    ## ENSG00000278384  0.200565196478864    0.841038574220432
    ## ENSG00000276345  -1.71475160732598   0.0863907773559434
    ## ENSG00000271254  -4.31407063627813 1.60275699407084e-05
    ##                                 padj       symbol      entrez
    ##                            <numeric>  <character> <character>
    ## ENSG00000279457      0.6865547778329           NA          NA
    ## ENSG00000187634  0.00515718149494361       SAMD11      148398
    ## ENSG00000188976  1.7654890539332e-35        NOC2L       26155
    ## ENSG00000187961 1.13412993107607e-07       KLHL17      339451
    ## ENSG00000187583    0.919030615571381      PLEKHN1       84069
    ## ...                              ...          ...         ...
    ## ENSG00000273748   0.0479091179108357           NA          NA
    ## ENSG00000278817    0.809772069001658 LOC102724770   102724770
    ## ENSG00000278384    0.892653531513564           NA          NA
    ## ENSG00000276345    0.139761501281218 LOC107987373   107987373
    ## ENSG00000271254 4.53647639305091e-05 LOC102724250   102724250

``` r
res$name =   mapIds(org.Hs.eg.db,
                    keys=row.names(res),
                    keytype="ENSEMBL",
                    column="GENENAME",
                    multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
res
```

    ## log2 fold change (MLE): condition hoxa1 kd vs control sirna 
    ## Wald test p-value: condition hoxa1 kd vs control sirna 
    ## DataFrame with 15975 rows and 9 columns
    ##                         baseMean     log2FoldChange              lfcSE
    ##                        <numeric>          <numeric>          <numeric>
    ## ENSG00000279457 29.9135794276176   0.17925708367269  0.324821565250145
    ## ENSG00000187634 183.229649921658  0.426457118403306  0.140265820376892
    ## ENSG00000188976 1651.18807619944 -0.692720464846366 0.0548465415913946
    ## ENSG00000187961 209.637938486147  0.729755610585225  0.131859899969345
    ## ENSG00000187583 47.2551232589398 0.0405765278756312  0.271892808601774
    ## ...                          ...                ...                ...
    ## ENSG00000273748 35.3026523877463  0.674387102558605  0.303666187454139
    ## ENSG00000278817 2.42302393023632 -0.388988266500022   1.13039377720313
    ## ENSG00000278384 1.10179649846993  0.332990658240634    1.6602614216556
    ## ENSG00000276345 73.6449563127136 -0.356180759105173   0.20771565839825
    ## ENSG00000271254 181.595902546813 -0.609666545167282  0.141320482803512
    ##                               stat               pvalue
    ##                          <numeric>            <numeric>
    ## ENSG00000279457  0.551863246932648    0.581042050747032
    ## ENSG00000187634   3.04034951107421  0.00236303749730996
    ## ENSG00000188976  -12.6301576133481 1.43989540156582e-36
    ## ENSG00000187961   5.53432552849563  3.1242824807768e-08
    ## ENSG00000187583  0.149237223611387    0.881366448669148
    ## ...                            ...                  ...
    ## ENSG00000273748   2.22081723425481    0.026363342804782
    ## ENSG00000278817 -0.344117487502873    0.730757932009185
    ## ENSG00000278384  0.200565196478864    0.841038574220432
    ## ENSG00000276345  -1.71475160732598   0.0863907773559434
    ## ENSG00000271254  -4.31407063627813 1.60275699407084e-05
    ##                                 padj       symbol      entrez
    ##                            <numeric>  <character> <character>
    ## ENSG00000279457      0.6865547778329           NA          NA
    ## ENSG00000187634  0.00515718149494361       SAMD11      148398
    ## ENSG00000188976  1.7654890539332e-35        NOC2L       26155
    ## ENSG00000187961 1.13412993107607e-07       KLHL17      339451
    ## ENSG00000187583    0.919030615571381      PLEKHN1       84069
    ## ...                              ...          ...         ...
    ## ENSG00000273748   0.0479091179108357           NA          NA
    ## ENSG00000278817    0.809772069001658 LOC102724770   102724770
    ## ENSG00000278384    0.892653531513564           NA          NA
    ## ENSG00000276345    0.139761501281218 LOC107987373   107987373
    ## ENSG00000271254 4.53647639305091e-05 LOC102724250   102724250
    ##                                                                     name
    ##                                                              <character>
    ## ENSG00000279457                                                       NA
    ## ENSG00000187634                 sterile alpha motif domain containing 11
    ## ENSG00000188976 NOC2 like nucleolar associated transcriptional repressor
    ## ENSG00000187961                              kelch like family member 17
    ## ENSG00000187583                 pleckstrin homology domain containing N1
    ## ...                                                                  ...
    ## ENSG00000273748                                                       NA
    ## ENSG00000278817                                            protein DGCR6
    ## ENSG00000278384                                                       NA
    ## ENSG00000276345                 39S ribosomal protein L23, mitochondrial
    ## ENSG00000271254                 neuroblastoma breakpoint family member 1

# Pathway analysis

``` r
BiocManager::install( c("pathview", "gage", "gageData") )
```

``` r
library(pathview)
```

    ## ##############################################################################
    ## Pathview is an open source software package distributed under GNU General
    ## Public License version 3 (GPLv3). Details of GPLv3 is available at
    ## http://www.gnu.org/licenses/gpl-3.0.html. Particullary, users are required to
    ## formally cite the original Pathview paper (not just mention it) in publications
    ## or products. For details, do citation("pathview") within R.
    ## 
    ## The pathview downloads and uses KEGG data. Non-academic uses may require a KEGG
    ## license agreement (details at http://www.kegg.jp/kegg/legal.html).
    ## ##############################################################################

``` r
library(gage)
library(gageData)
```

``` r
data(kegg.sets.hs)
data(sigmet.idx.hs)

# Focus on signaling and metabolic pathways only
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

# Examine the first 3 pathways
head(kegg.sets.hs, 3)
```

    ## $`hsa00232 Caffeine metabolism`
    ## [1] "10"   "1544" "1548" "1549" "1553" "7498" "9"   
    ## 
    ## $`hsa00983 Drug metabolism - other enzymes`
    ##  [1] "10"     "1066"   "10720"  "10941"  "151531" "1548"   "1549"  
    ##  [8] "1551"   "1553"   "1576"   "1577"   "1806"   "1807"   "1890"  
    ## [15] "221223" "2990"   "3251"   "3614"   "3615"   "3704"   "51733" 
    ## [22] "54490"  "54575"  "54576"  "54577"  "54578"  "54579"  "54600" 
    ## [29] "54657"  "54658"  "54659"  "54963"  "574537" "64816"  "7083"  
    ## [36] "7084"   "7172"   "7363"   "7364"   "7365"   "7366"   "7367"  
    ## [43] "7371"   "7372"   "7378"   "7498"   "79799"  "83549"  "8824"  
    ## [50] "8833"   "9"      "978"   
    ## 
    ## $`hsa00230 Purine metabolism`
    ##   [1] "100"    "10201"  "10606"  "10621"  "10622"  "10623"  "107"   
    ##   [8] "10714"  "108"    "10846"  "109"    "111"    "11128"  "11164" 
    ##  [15] "112"    "113"    "114"    "115"    "122481" "122622" "124583"
    ##  [22] "132"    "158"    "159"    "1633"   "171568" "1716"   "196883"
    ##  [29] "203"    "204"    "205"    "221823" "2272"   "22978"  "23649" 
    ##  [36] "246721" "25885"  "2618"   "26289"  "270"    "271"    "27115" 
    ##  [43] "272"    "2766"   "2977"   "2982"   "2983"   "2984"   "2986"  
    ##  [50] "2987"   "29922"  "3000"   "30833"  "30834"  "318"    "3251"  
    ##  [57] "353"    "3614"   "3615"   "3704"   "377841" "471"    "4830"  
    ##  [64] "4831"   "4832"   "4833"   "4860"   "4881"   "4882"   "4907"  
    ##  [71] "50484"  "50940"  "51082"  "51251"  "51292"  "5136"   "5137"  
    ##  [78] "5138"   "5139"   "5140"   "5141"   "5142"   "5143"   "5144"  
    ##  [85] "5145"   "5146"   "5147"   "5148"   "5149"   "5150"   "5151"  
    ##  [92] "5152"   "5153"   "5158"   "5167"   "5169"   "51728"  "5198"  
    ##  [99] "5236"   "5313"   "5315"   "53343"  "54107"  "5422"   "5424"  
    ## [106] "5425"   "5426"   "5427"   "5430"   "5431"   "5432"   "5433"  
    ## [113] "5434"   "5435"   "5436"   "5437"   "5438"   "5439"   "5440"  
    ## [120] "5441"   "5471"   "548644" "55276"  "5557"   "5558"   "55703" 
    ## [127] "55811"  "55821"  "5631"   "5634"   "56655"  "56953"  "56985" 
    ## [134] "57804"  "58497"  "6240"   "6241"   "64425"  "646625" "654364"
    ## [141] "661"    "7498"   "8382"   "84172"  "84265"  "84284"  "84618" 
    ## [148] "8622"   "8654"   "87178"  "8833"   "9060"   "9061"   "93034" 
    ## [155] "953"    "9533"   "954"    "955"    "956"    "957"    "9583"  
    ## [162] "9615"

The main gage() function requires a named vector of fold changes, where
the names of the values are the Entrez gene IDs.

``` r
foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)
```

    ##        <NA>      148398       26155      339451       84069       84808 
    ##  0.17925708  0.42645712 -0.69272046  0.72975561  0.04057653  0.54281049

``` r
# Get the results
keggres = gage(foldchanges, gsets=kegg.sets.hs)
```

now lets look at the object returned from gage()

``` r
attributes(keggres)
```

    ## $names
    ## [1] "greater" "less"    "stats"

``` r
# Look at the first few down (less) pathways
head(keggres$less)
```

    ##                                          p.geomean stat.mean        p.val
    ## hsa04110 Cell cycle                   8.995727e-06 -4.378644 8.995727e-06
    ## hsa03030 DNA replication              9.424076e-05 -3.951803 9.424076e-05
    ## hsa03013 RNA transport                1.246882e-03 -3.059466 1.246882e-03
    ## hsa03440 Homologous recombination     3.066756e-03 -2.852899 3.066756e-03
    ## hsa04114 Oocyte meiosis               3.784520e-03 -2.698128 3.784520e-03
    ## hsa00010 Glycolysis / Gluconeogenesis 8.961413e-03 -2.405398 8.961413e-03
    ##                                             q.val set.size         exp1
    ## hsa04110 Cell cycle                   0.001448312      121 8.995727e-06
    ## hsa03030 DNA replication              0.007586381       36 9.424076e-05
    ## hsa03013 RNA transport                0.066915974      144 1.246882e-03
    ## hsa03440 Homologous recombination     0.121861535       28 3.066756e-03
    ## hsa04114 Oocyte meiosis               0.121861535      102 3.784520e-03
    ## hsa00010 Glycolysis / Gluconeogenesis 0.212222694       53 8.961413e-03

Now, let’s try out the pathview() function from the pathview package to
make a pathway plot with our RNA-Seq expression results shown in color.

``` r
pathview(gene.data=foldchanges, pathway.id="hsa04110")
```

    ## 'select()' returned 1:1 mapping between keys and columns

    ## Info: Working in directory C:/Users/wanxu/Desktop/BIMM 143/LE 10/bimm143/class15

    ## Info: Writing image file hsa04110.pathview.png

This downloads the pathway figure data from KEGG and adds our results to
it.

You can play with the other input arguments to pathview() to change the
display in various ways including generating a PDF graph.

``` r
# A different PDF based output of the same data
pathview(gene.data = foldchanges, pathway.id = "hsa04110", kegg.native = FALSE)
```

    ## 'select()' returned 1:1 mapping between keys and columns

    ## Info: Working in directory C:/Users/wanxu/Desktop/BIMM 143/LE 10/bimm143/class15

    ## Info: Writing image file hsa04110.pathview.pdf
