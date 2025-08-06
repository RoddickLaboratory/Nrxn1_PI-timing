library(here)
library(rmarkdown)
library(bookdown)

renv::restore()

if (!(file.exists(here("data", "longDataNoTail.csv")))) {
  cat(
    "No long data found, running PIAnalysis.m with MATLAB to generate it.\nThis may take a few minutes..."
  )
  system(
    "matlab -nodisplay -nosplash -r \"run('./code/PIAnalysis.m'); exit;\"",
    intern = TRUE
  )
}

if (!(file.exists(here("data", "PI_timing_poke_state.RDS")))) {
  cat(
    "No poke state data found, running Activity.R to generate it.\nThis may take a few minutes..."
  )
  source(here("code", "Activity.R"))
}

render(
  here("code", "PI_timing.Rmd"),
  output_format = pdf_document2(
    toc = FALSE,
    latex_engine = "xelatex",
    includes = includes(in_header = here("code", "PI_timing_preamble.tex"))
  ),
  output_file = here(
    "docs",
    stringr::str_c("PI_timing ", lubridate::today(), ".pdf")
  )
)
