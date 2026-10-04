# Load configuration from a file

Load configuration from a file

## Usage

``` r
load_config(filename)
```

## Arguments

- filename:

  A file written by \[save_config()\] (\`.txt\` or \`.rds\`).

## Value

Invisibly, the loaded configuration list.

## Examples

``` r
f <- tempfile(fileext = ".txt")
save_config(f)
#> ✔ Configuration saved to /tmp/RtmpuvFqob/file1a62372ba3f9.txt
load_config(f)
#> ✔ Configuration loaded from /tmp/RtmpuvFqob/file1a62372ba3f9.txt
```
