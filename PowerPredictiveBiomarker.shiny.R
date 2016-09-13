#' Run for Shiny Applictaion for Power Calculation of Predictive Biomarker and generate a statistical plan to justify the sample size
#' @export
PowerPredictiveBiomarker.shiny <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "PowerPredictiveBiomarker")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `PowerPredictiveBiomarker`.", call. = FALSE)
  }

  shiny::runApp(paste(appDir,'/predictive_power_shiny_clean.R',sep=''), launch.browser =T,host = getOption( "127.0.0.1"))
}


