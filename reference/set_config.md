# Set configuration parameters for HaploTraitR

Set configuration parameters for HaploTraitR

## Usage

``` r
set_config(config = list())
```

## Arguments

- config:

  A named list of configuration parameters. See \[list_config()\] for
  the available parameters.

## Value

Invisibly, the previous values of the changed parameters.

## Examples

``` r
set_config(list(dist_threshold = 2000000, dist_cluster_count = 10))
reset_config()
#> ✔ HaploTraitR configuration reset to default.
```
