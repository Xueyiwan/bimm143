Class 6 R functions
================
Xueyi Wan
2019/10/17

# This is H1 (heading 1)

This is my work from class 6 in **BIMM143**.

## This is heading 2

### This is heading 3, etc.

``` r
# this is to demo a code chunk
plot(1:10)
```

![](Class06_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

## Practice reading files (again…)

``` r
# directly read online file
test1 <- read.csv("https://bioboot.github.io/bimm143_S19/class-material/test1.txt")

# or save the online file: save link as > save inside the R project
test1 <- read.csv("test1.txt")

test2 <- read.table("test2.txt", sep = "$", header = TRUE)

test3 <- read.delim("test3.txt", sep = "\t", header = FALSE)
```

## function section in class slides

``` r
add <- function(x, y=1) {
# Sum the input x and y
x + y
}
```

## our wee function is vectorized :-)

``` r
# using the function
add(1)
```

    ## [1] 2

``` r
# overridding the y variable with 5
add(5, 5)
```

    ## [1] 10

``` r
# add(5, 5, 5) will give you an error message
# add(5, "barry") will also give error b/c its non-numerical argument

# add up for every elements in the input/vector
add(c(1, 2, 3))
```

    ## [1] 2 3 4

``` r
add(c(1, 2, 3), 4)
```

    ## [1] 5 6 7

A new function to re-scale data

``` r
## You need a “name”, “arguments” and “body" for a function
rescale <- function(x) {
rng <-range(x)
(x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
# Test on a small example where you know the answer
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
# test more: how would you get your functio to work here
rescale(c(1, 2, NA, 3, 10))
```

    ## [1] NA NA NA NA NA

``` r
# what should your function do here
# rescale(c(1, 10, "string"))
# will get error message: non-numeric argument to binary operator
```

``` r
# ignore the NA value in the vector then test again
x <- c(1, 2, NA, 3, 10)
rng <- range(x, na.rm = TRUE)
rng
```

    ## [1]  1 10

``` r
rescale2 <- function(x) {
rng <-range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}
```

``` r
rescale2(c(1, 2, NA, 3, 10))
```

    ## [1] 0.0000000 0.1111111        NA 0.2222222 1.0000000

Too much…

``` r
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
  
  rng <-range(x, na.rm=na.rm)
  print("Hello")
  
  answer <- (x - rng[1]) / (rng[2] - rng[1])
  
  print("is it me you are looking for?")
  
  if(plot) {
    plot(answer, typ="b", lwd=4)
  }
  print("I can see it in ...")
  return(answer)
}
```

``` r
rescale3(1:10, plot = TRUE)
```

    ## [1] "Hello"
    ## [1] "is it me you are looking for?"

![](Class06_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

    ## [1] "I can see it in ..."

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000
