# Interval Timing is altered in male Nrxn1+/- mice: A Model of Autism Spectrum Disorder

Kyle M Roddick, Department of Psychology, Mount Allison University, ORCID: 0000-0003-2701-8166, corresponding author, kroddick@mta.ca
Elias B Habib, Department of Epidemiology, Dalhousie University, ORCID: 0000-0002-9911-8679
Richard E Brown, Department of Psychology and Neuroscience, Dalhousie University, ORCID: 0000-0002-4762-8175
Fuat Balcı, Department of Biological Sciences, University of Manitoba, ORCID: 0000-0003-3390-9352

Dates of Data Collection: Oct 2019 - Dec 2023

Location: Department of Psychology and Neuroscience, Dalhousie University, Halifax, NS, Canada

Funding was provided by the Natural Sciences and Engineering Research Council of Canada to REB (A7441), the Simons Foundation Autism Research Initiative to REB (608066), and NSERC Discovery Grant (RGPIN-2021-03334) to FB.

## Instructions:

Run the code/render.R script to render the Rmarkdown file PI_timing.Rmd into a .pdf report. This output file will be placed in the docs/ folder. This will also restore the R environment and run the necessary scripts to generate the data files if they are missing. If data/longDataNoTail.csv is missing, the MATLAB script code/PIAnalysis.m will be run (requires MATLAB with the Curve Fitting Toolbox, the Signal Processing Toolbox, and struct2csv add-ons) to regenerate it from the data/peak_data/ files. If data/PI_timing_poke_state.RDS is missing, the R script code/Activity.R will be run to regenerate it from the data/peak_data/ files.

## Files:

.gitignore - File listing the files and directories to ignore when using git for version control.
.renvignore - File listing the files and directories to ignore when using the renv package for R package management.
.Rprofile - R profile file that sets the default options for the R session.
code/ - Directory containing code for the analysis of the data.
  Activity.R - R script with code to analyze the activity of the mice. Outputs the PI_timing_poke_state.RDS file used by the Rmarkdown file.
  apa.csl - Citation style language file describing how the citations and bibliography should be formatted.
  getStartStop.m - MATLAB script with code to determine the start and stop times for the single trial analysis.
  NewExpWaldFittingforAutismResearchRevision.m - MATLAB script with code to fit Exponential Wald distributions to the single trial data.
  NewWaldFittingforAutismResearchRevision.m - MATLAB script with code to fit Wald distributions to the single trial data.
  PI_timing_preamble.tex - Tex file that provides the preamble for the rendering of the title, authors, affiliations, the header in the output file.
  PI_timing.Rmd - Rmarkdown file containing the code for the analysis and visualization of the data, renders the output .pdf file when render.R is run.
  PIAnalysis.m - MATLAB script with code to analyze the data. Outputs the longDataNoTail.csv data file used by the Rmarkdown file. Requires the Curve Fitting Toolbox and the Signal Processing Toolbox.
  references.bib - BibTeX file containing metadata on references cited.
  render.R - R script with code to render the Rmarkdown file into the output .pdf file. Will also restore the renv environment and run the Activity.R and PIAnalysis.m scripts if the data files they output are missing.
data/ - Directory containing the data files.
  genotypes.csv - Data file listing the genotypes of the mice. Columns are the mouse subject number (mouse), the genotype (geno).
    For the genotypes, wt indicates wildtype mice, tg indicates Nrxn1+/- mice, and ds5/- indicates Nrxn1ΔS5/- mice.
  longDataNoTail.csv - Data file containing the long format data for the analysis. It is the output of the PIAnalysis.m script and read into the PI_timing.Rmd file. 
  peak_data/ - Directory containing the data files for the peak interval test. Each data file is a .csv file with the subject ID, the name of the program used to test the mouse, and the trial number([subject ID] [program] Trial [trial].csv). Following a line of metadata, each data file contains a log of the events that occurred during testing (both programmed events and inputs from the instruments in the test chamber) in three columns. 
      - The first column is the time in milliseconds since the program started. 
      - The second column indicates the source type of the event; DIGITAL_IN indicates it is the result of input from the test chamber, while PROGRAM, PRGRAMMATIC, and ALWAYS indicate events initiated by the program running the chamber. 
      - The third column indicates the event. "ChoiceWallIR" and "TrayWallIR" are changes in the state of the infrared beams agains the choice wall and the tray wall of the test chamber respectivly. "Hole_NP" are changes of state of the infrared beam inside the _th nose poke hole in the choice wall. Other events are programmatic events initiated by the program running the test chamber.
      - The fourth column indicates the state of the event for DIGITAL_IN events. For the infrared beam events, a value of -1 indicates the beam was broken (i.e., the mouse interrupted the beam) while a value of 0 indicates the beam was unbroken. PROGRAMMATIC and ALWAYS events have a value of -1.
    Contains 1438 data files.
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

## Methodological Information:

The apparatus and procedures were the same as those described by Gür et al. (2019). Testing was conducted in a mouse nine-hole box (Cambridge Cognition Ltd., England) with a lick tube attached to a peristaltic pump, house lights, and speakers (mounted on both sides of the inner walls). The box was placed in a sound- and light-attenuating chamber. The test chamber contained a grid floor with a removable tray underneath. Six of nine holes (Holes 1–3, 7–9) were plugged, and the remaining three holes in the middle (Holes 4–6) were left open. Each of the three holes and the reinforcement tube/tray could be illuminated, and nose pokes to the open holes were detected via infrared beams. Computer software (Cambridge Cognition Ltd., England) was used to control the test box and record time-stamped nose pokes.

Reference for Apparatus:
Gür, E., Fertan, E., Kosel, F., Wong, A. A., Balcı, F., & Brown, R. E. (2019). Sex differences in the timing behavior performance of 3xTg-AD and wild-type mice in the peak interval procedure. Behavioural Brain Research, 360, 235–243. https://doi.org/10.1016/j.bbr.2018.11.047
