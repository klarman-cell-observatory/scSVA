#' Run scSVA
#'
#' @param background Whether to run Shiny App in background (access to the console while running scSVA), default=FALSE
#' @return Open Shiny App
#' @examples
#' scSVA()
#' @export

scSVA <- function(background=FALSE)
{
  if(background){
    PATH  =  paste0("R -e \"shiny::runApp(\'",  system.file("scSVA", package = "scSVA"), "\', launch.browser = TRUE)\"")
    system(PATH, wait = FALSE)
  } else {
    shiny::runApp(appDir = system.file("scSVA", package = "scSVA"), launch.browser =TRUE)
  }
}

