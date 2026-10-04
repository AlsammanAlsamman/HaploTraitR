# Save configuration to a file

Writes a human-readable \`name=value\` text file and an \`.rds\` copy
that preserves the value types.

## Usage

``` r
save_config(filename)
```

## Arguments

- filename:

  The name of the text file to save the configuration to

## Value

Invisibly, the file name.

## Examples

``` r
f <- tempfile(fileext = ".txt")
save_config(f)
#> ✔ Configuration saved to /tmp/RtmpuvFqob/file1a6227245468.txt
```
