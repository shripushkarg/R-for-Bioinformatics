library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
# library(purrr)


# ----------------------- Helper Functions to Implement ------------------------

#' Read the expression data "csv" file.
#'
#' Function to read microarray expression data stored in a csv file. The
#' function should return a sample x gene tibble, with an extra column named
#' "subject_id" that contains the geo accession ids for each subject.
#'
#' @param filename (str): the file to read.
#'
#' @return
#' @export
#'
#' @examples expr_mat <- read_expression_table('example_intensity_data_subset.csv')
read_expression_table <- function(filename) {
  expression_data <- read_csv('data/example_intensity_data_subset.csv')
    expression_mat_1 = expression_data %>% pivot_longer(cols= -1)
    expression_mat_1 = rename(expression_mat_1,subject_id=name)
    return (expression_mat_1)
}



#' Load Metadata from Specified CSV File
#'
#' This function reads the provided CSV file into a dataframe.
#'
#' @param filepath (character) The path to the CSV file.(data/proj_metadata.csv)
#'
#' @return A dataframe containing the loaded metadata.
#'

load_metadata <- function(filepath) {
  
  # TODO: Use appropriate function to read in the CSV file
  #metadata <-
  metadata <- read_csv('data/proj_metadata.csv')
  
  return(metadata)
    
    # Return the loaded metadata
    #return(metadata)
}




#' Replaces all '.' in a string with '_'
#'
#' @param str String to operate upon.
#'
#' @return reformatted string.
#' @export
#'
#' @examples
#' period_to_underscore("foo.bar")
#' "foo_bar"
period_to_underscore <- function(str) {

  return (str_replace_all(str, "\\.", "_"))
}


# rename variables:
# Age_at_diagnosis to Age
# SixSubtypesClassification to Subtype
# normalizationcombatbatch to Batch

#' Rename and select specified columns.
#'
#' Function to rename Age_at_diagnosis, SixSubtypesClassification, and
#' normalizationcombatbatch columns to Age, Subtype, and Batch, respectively. A
#' subset of the data should be returned, containing only the Sex, Age, TNM_Stage,
#' Tumor_Location, geo_accession, KRAS_Mutation, Subtype, and Batch columns.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) renamed and subsetted tibble
#' @export
#'
#' @examples rename_and_select(metadata)
#' 
#' 
rename_and_select <- function(data) {
  data_1 = rename(data, Age = Age_at_diagnosis, Subtype = SixSubtypesClassification, Batch = normalizationcombatbatch)
  data_1 <- dplyr::select(data_1,Sex,Age,TNM_Stage,Tumor_Location,geo_accession,KRAS_Mutation,Subtype,Batch)
  return (data_1)
}


#' Create new "Stage" column containing "stage " prefix.
#'
#' Creates a new column "Stage" with elements following a "stage x" format, where
#' `x` is the cancer stage data held in the existing TNM_Stage column. Stage
#' should have a factor data type.
#'
#' @param data  (tibble) metadata information for each sample
#'
#' @return (tibble) updated metadata with "Stage" column
#' @export
#'
#' @examples metadata <- stage_as_factor(metadata)
stage_as_factor <- function(data) {
    Stage <- paste("stage",data$TNM_Stage)
    data_modified <- add_column(data, Stage) %>% mutate_at(vars(Stage),factor)
  return (data_modified)
}


#' Calculate age of samples from a specified sex.
#'
#' @param data (tibble) metadata information for each sample
#' @param sex (str) which sex to calculate mean age. Possible values are "M"
#' and "F"
#'
#' @return (float) mean age of specified samples
#' @export
#'
#' @examples mean_age_by_sex(metadata, "F")
mean_age_by_sex <- function(data, sex) {
  data_1 <- dplyr::filter(data, Sex == sex)
  return (mean(data_1$Age))
}


#' Calculate average age of samples within each cancer stage. Stages should be
#' from the newly created "Stage" column.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) summarized tibble containing average age for all samples from
#' each stage. Name the newly created column containing the average, 'mean_avg'
#' @export
#'
#' @examples age_by_stage(data)
age_by_stage <- function(data) {
  avg_age_by_stage <- aggregate(Age ~ Stage, data = data, FUN = mean)
  result_data <- avg_age_by_stage %>%
    group_by(Stage) %>%
    summarise(mean_avg = mean(Age, na.rm = TRUE))
  return (result_data)
}

#' Create a cross tabulated table for Subtype and Stage using dplyr methods.
#'
#' @param data (tibble) metadata information for each sample
#'
#' @return (tibble) table where rows are the cancer stage of each sample, and the
#' columns are each cancer subtype. Elements represent the number of samples from
#' the corresponding stage and subtype. If no instances of a specific pair are
#' observed, a zero entry is expected.
#' @export
#'
#' @examples cross_tab <- dplyr_cross_tab(metadata)
subtype_stage_cross_tab <- function(data) {
    subtype = table(data)
    cross_table = tibble(Stage = rownames(subtype), C3 = subtype[1:4], C4 = subtype[5:8])
    return (cross_table)
}

#' Summarize average expression and probe variability over expression matrix.
#'
#' @param exprs An (n x p) expression matrix, where n is the number of samples,
#' and p is the number of probes.
#'
#' @return A summarized tibble containing `main_exp`, `variance`, and `probe`
#' columns documenting average expression, probe variability, and probe ids,
#' respectively.
summarize_expression <- function(exprs) {
  summary <-
    dplyr::summarise(exprs,
                     mean_exp = colMeans(exprs[,-1]),
                     variance = apply(exprs[,-1], 2, var),
                     probe = colnames(exprs[-1])
    )
  return(summary)
}
