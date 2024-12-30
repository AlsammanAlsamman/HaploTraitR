#' Create a unique result folder for HaplotraitR
#' @param base_folder_name The base name for the result folder
#' @param location The location where the folder should be created
#' @return The path to the created folder
#' @export
#' @examples
#' result_folder <- create_unique_result_folder()
#' cat("Results will be stored in:", result_folder, "\n")
create_unique_result_folder <- function(base_folder_name = "haplotraitR_run", location = NULL) {
  # If location is not provided, use the home directory and notify the user
  if (is.null(location)) {
    location <- path.expand("~")
    message("No location provided. Using the home directory: ", location)
  }

  # Initialize the folder name and counter
  folder_name <- base_folder_name
  counter <- 1

  # Generate the full path to the folder
  result_folder <- file.path(location, folder_name)

  # Check if the folder already exists, and if so, increment the counter
  while (dir.exists(result_folder)) {
    counter <- counter + 1
    folder_name <- paste(base_folder_name, counter, sep = "_")
    result_folder <- file.path(location, folder_name)
  }

  # Create the directory
  dir.create(result_folder, showWarnings = FALSE)

  # Return the path to the created directory
  return(result_folder)
}

