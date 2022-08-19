
#' List files
#'
#' A function to serve as a wrapper around the base function list.files() to ensure returned results consist of ONLY file names
#'
#' @param path Character. Path to folder to be searched
#' @param pattern Character. Pattern to be used to select subsets of files
#' @param all.files Logical. Should all files be searched for? Default is FALSE
#' @param full.names Logical. Should the full pathnames of files be returned? Default is TRUE
#' @param recursive Logical. Should sub-folders be searched recursive for files matching the search pattern? Default is FALSE
#' @param ignore.case Logical. Should chase in file names be ignored? Default is FALSE
#' @param include.dirs Logical. Should folders be included int eh returned results? Default is FALSE
#' @param no.. Logical. Should both "." and ".." be excluded also from non-recursive listings?
#' @param incl_dirs Logical. Should subdirectory names be included in recursive listings? (They always are in non-recursive ones).
#'
#' @details From: https://stackoverflow.com/questions/22069095/r-get-list-of-files-but-not-of-directories
#' and answer by Dunois posted 2021-01-24
#'
#' @return
#' @export
#'
#' @examples
list_files <- function(path = ".",
                       pattern = NULL,
                       all.files = FALSE,
                       full.names = TRUE,
                       recursive = FALSE,
                       ignore.case = FALSE,
                       include.dirs = FALSE,
                       no.. = FALSE,
                       incl_dirs = FALSE)
  {

  #Set incl_dirs = TRUE to revert to default list.files() behavior.

  if(path == ".") { path = getwd() }

  #Include directories if recursive is set.
  if(incl_dirs & recursive) { include.dirs = TRUE }

  #Needs to have full.names = TRUE in order to get full path to pass to dir.exists().
  files <- list.files(path = path, pattern = pattern, all.files = all.files, full.names = TRUE,
                      recursive = recursive, ignore.case = ignore.case, include.dirs = include.dirs,
                      no.. = no..)

  if(!incl_dirs){
    files <- files[!dir.exists(files)]
  }

  if(!full.names){
    return(basename(files))
  } else{
    return(files)
  }

}


