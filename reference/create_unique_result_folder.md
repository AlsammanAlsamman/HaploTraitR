# Create a unique result folder for HaploTraitR

Create a unique result folder for HaploTraitR

## Usage

``` r
create_unique_result_folder(
  base_folder_name = "haplotraitR_run",
  location = NULL
)
```

## Arguments

- base_folder_name:

  The base name for the result folder

- location:

  The location where the folder should be created (default: the current
  working directory)

## Value

The path to the created folder

## Examples

``` r
result_folder <- create_unique_result_folder(location = tempdir())
```
