#' @title Read individual growth plate
#'
#' @description This function receives a file name that contains plate reader OD data. It returns a data.frame with the parsed data.
#'
#' @param plate_file A character variable that is the file name for the plate OD data. Should contain a Time column and a column for each well.
#' @param path_target The path to plate_file.
#' @param time_vector A numeric vector that indicates the Time at which each measurement was taken.
#' @param sep The delimiter character that separates columns in plate_file.

#' @return A data.frame with the data from plate_file. The Time column contains the values in time_vector, and the wells are renamed to be {plate_name}_{well}.
#'
#' @name read_individual_plate
#' @rdname read_individual_plate
#' @export
read_individual_plate <- function(plate_file, path_target, time_vector=c(), sep='\t'){
  
  # Read the data for the file
  plate_df <- read.delim(file.path(path_target, plate_file), sep=sep)
  
  # Make sure there is a Time column
  if(!"Time" %in% colnames(plate_df)){
    stop(paste("No column named 'Time' found in", plate_file))
  }
  
  
  # Assign the values of the common time_vector to the Time column
  if(length(time_vector) > 0){
    if(nrow(plate_df) != length(time_vector)){
      stop(paste("The length of time_vector does not equal the number of found in", plate_name))
      
    }else{
      plate_df$Time <- time_vector
    }
  }
  
  # Update the well identifiers to include the plate name
  plate_name <- tools::file_path_sans_ext(plate_file, compression = TRUE)
  plate_df <- plate_df %>% 
    rename_with(~paste0(plate_name, "_", .x, recycle0 = TRUE), -Time)
  
  
  return(plate_df)
}


#' @title Read a collection growth plate data.
#'
#' @description This function receives a path to a folder that contains data for several plates. Each plate should have its own file.
#'
#' @param path_target The path to plate_file.
#' @param time_vector A numeric vector that indicates the Time at which each measurement was taken. If empty then the Time column from one of the plates is taken as reference.
#' @param sep The delimiter character that separates columns of each file in path_target

#' @return A list of two dataframes. An "od_data" dataframe with a common time column and each growth curve in a separate column named with it's "curve_id" ({plate_name}_{well}). 
#' Also a "metadata" dataframe indicating the well and the original plate from which each "curve_id" is gathered from.
#'
#' @name read_multiple_plates
#' @rdname read_multiple_plates
#' @export
read_multiple_plates <- function(path_target, time_vector = c(), sep='\t'){
  # List files in path_target
  plate_files <- list.files(path_target)
  message(paste("The following OD data was found:", paste0(plate_files, collapse = ", ")))
  
  # Send message if time_vector is empty
  if(length(time_vector) == 0){
    message(paste("No common time_vector given. Will use the Time column from", plate_files[1], "as reference."))
    time_vector <- read_individual_plate(path_target = path_target,
                                         plate_file = plate_files[1],
                                         time_vector = time_vector,
                                         sep = sep)$Time
  }
  
  # Read the data from the other lists
  # TODO: can be optimized so as to not re-read the first file
  plates_read_list <- lapply(plate_files, read_individual_plate, path_target=path_target, time_vector=time_vector, sep='\t')
  
  # Combine all dataframes into one big dataframe
  complete_od <- reduce(plates_read_list, full_join, by="Time")
  
  # Create a mapping dataframe
  curve_ids <- colnames(complete_od)[colnames(complete_od) != "Time"]
  mapper_df <- data.frame(curve_id = curve_ids,
                          well = str_split_i(curve_ids, "_", -1)) %>% 
    mutate(plate_name = str_remove(curve_id, paste0("_", well)))
  
  # Return OD and metadata df
  
  return(list("od_data" = complete_od,
              "metadata" = mapper_df))
  
  
}