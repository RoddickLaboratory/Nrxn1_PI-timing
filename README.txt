# Interval Timing is altered in male Nrxn1+/- mice: A Model of Autism Spectrum Disorder

Kyle M Roddick, Department of Psychology, Mount Allison University
Elias B Habib, Department of Epidemiology, Dalhousie University
Richard E Brown, Department of Psychology and Neuroscience, Dalhousie University
Fuat Balcı, Department of Biological Sciences, University of Manitoba

Dates of Data Collection: Oct 2019 - Dec 2023

Location: Department of Psychology and Neuroscience, Dalhousie University, Halifax, NS, Canada

## Instructions:

Run the code/render.R script to render the Rmarkdown file PI_timing.Rmd into a .pdf report. This output file will be placed in the docs/ folder. This will also restore the R environment and run the necessary scripts to generate the data files if they are missing.

## Files:

.gitignore - File listing the files and directories to ignore when using git for version control.
.renvignore - File listing the files and directories to ignore when using the renv package for R package management.
.Rprofile - R profile file that sets the default options for the R session.
code/ - Directory containing code for the analysis of the data.
  Activity.R - R script with code to analyze the activity of the mice. Outputs the PI_timing_poke_state.RDS file used by the Rmarkdown file.
  apa.csl - Citation style language file describing how the citations and bibliography should be formatted.
  getStartStop.m - MATLAB script with code to determine the start and stop times for the single trial analysis.
  PI_timing_preamble.tex - Tex file that provides the preamble for the rendering of the title, authors, affiliations, the header in the output file.
  PI_timing.Rmd - Rmarkdown file containing the code for the analysis and visualization of the data, renders the output .pdf file when render.R is run.
  PIAnalysis.m - MATLAB script with code to analyze the data. Outputs the longDataNoTail.csv data file used by the Rmarkdown file. Requires the Curve Fitting Toolbox and the Signal Processing Toolbox.
  references.bib - BibTeX file containing metadata on references cited.
  render.R - R script with code to render the Rmarkdown file into the output .pdf file. Will also restore the renv environment and run the Activity.R and PIAnalysis.m scripts if the data files they output are missing.
data/ - Directory containing the data files.
  genotypes.csv - Data file listing the genotypes of the mice. Columns are the mouse subject number (mouse), the genotype (geno), the date of birth (dob), the date of testing (dot), the age in days at the start of testing (age), and whether the mouse is too old to include in the analysis (old).
  longDataNoTail.csv - Data file containing the long format data for the analysis. It is the output of the PIAnalysis.m script and read into the PI_timing.Rmd file. 
  peak_data/ - Directory containing the data files for the peak interval test. Each data file is a .csv file with the subject ID, the name of the program used to test the mouse, and the trial number([subject ID] [program] Trial [trial].csv). Following metadata, each data file contains a log of the events that occurred during testing (both programmed events and inputs from the instruments in the test chamber), the time since the start of the test when the events occurred.
  PI_timing_poke_state.RDS - R data file containing the output of the Activity.R script. It is read into the PI_timing.Rmd file.
docs/ - Directory containing the output .pdf report created by the Rmarkdown file.
figures/ - Directory containing the output image files created by the Rmarkdown file.
LICENSE - The CC-BY-4.0 license.
Nrxn1_PI-timing.Rproj - The R project file containing the settings for the R project.
renv/ - Directory containing the R environment settings for the project.
  .gitignore - File listing the files and directories to ignore when using git for version control.
  activate.R - R script that activates the renv environment when run.
  settings.json - File containing the settings for the renv package.
renv.lock - File containing the R package dependencies for the project, used by the renv package to restore the environment.
