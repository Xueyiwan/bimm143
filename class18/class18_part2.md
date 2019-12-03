class18\_part2
================

# Protein sequences from healthy and tumor tissue

The following sequences resulted from such an NGS analysis of patient
healthy and tumor tissue.

Your task is to identify tumor specific mutations that could potentially
be used for vaccine development.

> Q1. Identify sequence regions that contain all 9-mer peptides that are
> only found in the tumor.

Hint: You will need to first identify the sites of mutation in the above
sequences and then extract the surrounding subsequence region. This
subsequence should encompass all possible 9-mers in the tumor derived
sequence. In other words extract the subsequence from 8 residues before
and 8 residues after all point mutations in the tumor sequence.

``` r
library(bio3d)
seqs <- read.fasta("lecture18_sequences.fa")

# We can optionally align these sequences, but in this case the provided sequences are already aligned.
# seqs <- seqaln(seqs)

## will need to use MUSCLE to align
## download MUSCLE and call it in R terminal
```

Next we calculate identity per equivalent (i.e. aligned) position and
then use this information to find non identical sites that do not
contain gaps (i.e. indels).

``` r
## Calculate positional identity scores
ide <- conserv(seqs$ali, method="identity")
mutant.sites <- which(ide < 1) 

## Exclude gap positions from analysis
gaps <- gap.inspect(seqs)
mutant.sites <- mutant.sites[mutant.sites %in% gaps$f.inds]

mutant.sites
```

    ## [1]  41  65 213 259

We can use these indices in mutant.sites to extract subsequences as
required for the hands-on session. First however we come up with
suitable names for these subsequences based on the mutation. This will
help us later to make sense and keep track of our results.

``` r
## Make a "names" label for our output sequences (one per mutant)
mutant.names <- paste0(seqs$ali["P53_wt",mutant.sites],
                       mutant.sites,
                       seqs$ali["P53_mutant",mutant.sites])

mutant.names
```

    ## [1] "D41L"  "R65W"  "R213V" "D259V"

Now lets extract all 9-mer mutant encompassing sequences for each mutant
site. This is equivalent to finding the sequence region eight residues
before and eight residues after our mutation sites and outputting this
subsequence to a new FASTA file.

``` r
## Sequence positions surounding each mutant site
start.position <- mutant.sites - 8
end.position <-  mutant.sites + 8

# Blank matrix to store sub-sequences
store.seqs <- matrix("-", nrow=length(mutant.sites), ncol=17)
rownames(store.seqs) <- mutant.names

## Extract each sub-sequence
for(i in 1:length(mutant.sites)) {
  store.seqs[i,] <- seqs$ali["P53_mutant",start.position[i]:end.position[i]]
}

store.seqs
```

    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
    ## D41L  "S"  "P"  "L"  "P"  "S"  "Q"  "A"  "M"  "L"  "D"   "L"   "M"   "L"  
    ## R65W  "D"  "P"  "G"  "P"  "D"  "E"  "A"  "P"  "W"  "M"   "P"   "E"   "A"  
    ## R213V "Y"  "L"  "D"  "D"  "R"  "N"  "T"  "F"  "V"  "H"   "S"   "V"   "V"  
    ## D259V "I"  "L"  "T"  "I"  "I"  "T"  "L"  "E"  "V"  "-"   "-"   "-"   "-"  
    ##       [,14] [,15] [,16] [,17]
    ## D41L  "S"   "P"   "D"   "D"  
    ## R65W  "A"   "P"   "P"   "V"  
    ## R213V "V"   "P"   "Y"   "E"  
    ## D259V "-"   "-"   "-"   "-"

Finally lets output all these sequences to a FASTA file for further
analysis with the IEDB HLA binding prediction website.

``` r
## Output a FASTA file for further analysis
write.fasta(seqs=store.seqs, ids=mutant.names, file="subsequences.fa")
```

# Patient HLA typing results and HLA binding prediction

> Q2: Identify 9-mer peptides in the identified sequence regions unique
> to the tumor that can be potentially presented to T cells.

Hint: Use the IEDB HLA binding prediction server above to identify the
top ranked 9-mer peptides for each patient HLA (see above for HLA typing
results).

> Q3: Identify the top peptide for each patient HLA allele (see above).

Hint: You can download a CSV formated result file for all predictions
and use R or a spreadsheet application to identify the top ranked
peptides for each allele.

``` r
HLA_predictions <- read.csv("result.csv")
head(HLA_predictions)
```

    ##        allele seq_num start end length   peptide
    ## 1 HLA-A*02:01       3     1   9      9 YLDDRNTFV
    ## 2 HLA-B*35:01       3     8  16      9 FVHSVVVPY
    ## 3 HLA-B*07:02       1     1   9      9 SPLPSQAML
    ## 4 HLA-A*02:01       2     9  17      9 WMPEAAPPV
    ## 5 HLA-B*07:02       1     3  11      9 LPSQAMLDL
    ## 6 HLA-A*02:01       4     1   9      9 ILTIITLEV
    ##                                   method Percentile.Rank ann_ic50 ann_rank
    ## 1 Consensus (ann/comblib_sidney2008/smm)             0.2     4.05     0.02
    ## 2 Consensus (ann/comblib_sidney2008/smm)             0.2     5.41     0.02
    ## 3 Consensus (ann/comblib_sidney2008/smm)             0.4    43.33     0.18
    ## 4 Consensus (ann/comblib_sidney2008/smm)             0.4     5.78     0.05
    ## 5 Consensus (ann/comblib_sidney2008/smm)             0.5    35.11     0.14
    ## 6 Consensus (ann/comblib_sidney2008/smm)             0.7    35.81     0.39
    ##   smm_ic50 smm_rank comblib_sidney2008_score comblib_sidney2008_rank
    ## 1     3.90      0.2                 3.63e-06                     0.3
    ## 2     6.16      0.2                 0.000698                     7.3
    ## 3    50.11      0.4                 0.000137                     4.6
    ## 4    17.66      0.4                 9.93e-06                     0.7
    ## 5    66.51      0.5                 6.56e-05                     2.2
    ## 6    41.03      0.7                 1.91e-05                     1.1
    ##   netmhcpan_score netmhcpan_rank
    ## 1               -              -
    ## 2               -              -
    ## 3               -              -
    ## 4               -              -
    ## 5               -              -
    ## 6               -              -

``` r
# get the top 10 ranked peptides
HLA_top_ranks <- HLA_predictions[1:10, ]
View(HLA_top_ranks)
```
