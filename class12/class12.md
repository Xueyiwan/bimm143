class12
================
Xueyi Wan
2019/11/7

## Prepare for docking

We want to produce a protein-only PDB file and a drug only PDB file

``` r
library(bio3d)

# download the PDB file
get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

    ## [1] "./1hsg.pdb"

``` r
pdb <- read.pdb("1hsg.pdb")
protein <- atom.select(pdb, "protein", value = TRUE)
write.pdb(protein, file = "1hsg_protein.pdb")
```

and for the ligand

``` r
ligand <- atom.select(pdb, "ligand", value = TRUE)
write.pdb(ligand, file = "1hsg_ligand.pdb")
```

## Using AutoDockTools to setup protein docking input

## Prepare a docking configuration file

## Processing your docking results

``` r
library(bio3d)
res <- read.pdb("all.pdbqt", multi = TRUE)
write.pdb(res, "results.pdb")
```

``` r
ori <- read.pdb("1hsg_ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.649  4.206 11.110 10.529  4.840 10.932 10.993  3.655 10.996 11.222
    ## [11] 10.567 10.372 11.019 11.338  8.390  9.063  8.254  8.978

## Exploring the conformational dynamics of proteins

``` r
library(bio3d)
pdb <- read.pdb("1hel")
```

    ##   Note: Accessing on-line PDB file

``` r
modes <- nma(pdb)
```

    ##  Building Hessian...     Done in 0.01 seconds.
    ##  Diagonalizing Hessian...    Done in 0.11 seconds.

``` r
plot(modes, sse = pdb)
```

![](class12_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# Visualize NMA results
mktrj(modes, mode = 7, file = "nma_7.pdb")
```
