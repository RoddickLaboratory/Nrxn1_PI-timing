#### Activity Analysis ####

library(tidyverse)
library(here)

#### Read data file ####

# create list of files
fileList <- list.files(
  './data/peak_data',
  pattern = '\\.csv'
)

fileDat <- tibble(files = fileList)

fileDat <- fileDat %>%
  mutate(
    subj = str_extract(files, "^[:digit:]{4}"),
    session = str_extract(files, "[:digit:]{2}(?=\\.csv$)")
  )

# read data files and add subject and session numbers
dat <- lapply(fileDat$files, function(fileName) {
  tmpDat <- read_csv(
    here("data", "peak_data", fileName),
    skip = 1,
    col_types = 'dccdc',
    col_names = c('time', 'eventType', 'event', 'state', 'X')
  )[, c(-2, -5)]
  tmpDat$subj <- fileDat$subj[fileDat$files == fileName]
  tmpDat$session <- fileDat$session[fileDat$files == fileName]
  tmpDat
})

# combine files into one frame
dat <- do.call(rbind, dat)

# dat$session <- as.factor(dat$session)
dat <- dat %>%
  group_by(subj, session)

# read in genotypes
genoDat <- read_csv('data/genotypes.csv')
colnames(genoDat)[1] <- 'subj'
genoDat$subj <- as.factor(genoDat$subj)

#### Beam Breaks ####

beamDat <- dat %>%
  filter(
    event %in% c('ChoiceWallIR', 'TrayWallIR') & state == -1 & time != 0
  ) %>%
  mutate('breaks' = -cumsum(state))

sumBeam <- beamDat %>%
  filter(time <= 20 * 60 * 1000) %>%
  summarise(totalBreaks = max(breaks)) %>%
  left_join(genoDat)

#### Just Choice Wall ####

cwDat <- dat %>%
  filter(event == 'ChoiceWallIR' & state == -1 & time != 0) %>%
  mutate('breaks' = -cumsum(state))

sumCWBeam <- cwDat %>%
  filter(time <= 20 * 60 * 1000) %>%
  summarise(totalBreaks = max(breaks)) %>%
  left_join(genoDat)

#### Just Tray Wall ####

twDat <- dat %>%
  filter(event == 'TrayWallIR' & state == -1 & time != 0) %>%
  mutate('breaks' = -cumsum(state))

sumTWBeam <- twDat %>%
  filter(time <= 20 * 60 * 1000) %>%
  summarise(totalBreaks = max(breaks)) %>%
  left_join(genoDat)

#### Time at Choice Wall ####

cwTimeDat <- dat %>%
  filter(event == 'ChoiceWallIR') %>%
  mutate('time' = time / 1000) %>%
  mutate('timeIntervale' = time - lag(time, default = 0)) %>%
  filter(state == 0) %>%
  mutate('timeSum' = cumsum(timeIntervale))

sumCWTime <- cwTimeDat %>%
  summarise(totalTime = max(timeSum)) %>%
  left_join(genoDat)

#### Nose Pokes ####

pokeDat <- dat %>%
  filter(event == 'Hole5NP' & state == -1 & time != 0) %>%
  mutate('pokes' = -cumsum(state))

sumPokes <- pokeDat %>%
  filter(time <= 20 * 60 * 1000) %>%
  summarise(totalPokes = max(pokes)) %>%
  left_join(genoDat)

#### Add trial and ITI state ####

pokeStateDat <- dat %>%
  mutate('trialState' = 'ITI')

piTrialTimes <- dat %>%
  filter(event %in% c('PIEvent', 'PITimerEvent') & time <= (20 * 60 * 1000)) %>%
  mutate('timeEnd' = lead(time, default = (20 * 60 * 1000))) %>%
  filter(event == 'PIEvent')

fiTrialTIme <- dat %>%
  filter(event %in% c('FIEvent', 'Hole5Event') & time <= (20 * 60 * 1000)) %>%
  mutate('timeEnd' = lead(time, default = (20 * 60 * 1000))) %>%
  filter(event == 'FIEvent')

for (i in 1:nrow(piTrialTimes)) {
  tmpSubj <- piTrialTimes$subj[i]
  tmpSession <- piTrialTimes$session[i]
  tmpStart <- piTrialTimes$time[i]
  tmpEnd <- piTrialTimes$timeEnd[i]
  pokeStateDat$trialState[
    pokeStateDat$subj == tmpSubj &
      pokeStateDat$session == tmpSession &
      pokeStateDat$time >= tmpStart &
      pokeStateDat$time < tmpEnd
  ] <- 'PITrial'
}

for (i in 1:nrow(fiTrialTIme)) {
  tmpSubj <- fiTrialTIme$subj[i]
  tmpSession <- fiTrialTIme$session[i]
  tmpStart <- fiTrialTIme$time[i]
  tmpEnd <- fiTrialTIme$timeEnd[i]
  pokeStateDat$trialState[
    pokeStateDat$subj == tmpSubj &
      pokeStateDat$session == tmpSession &
      pokeStateDat$time >= tmpStart &
      pokeStateDat$time < tmpEnd
  ] <- 'FITrial'
}

pokeStateDat <- pokeStateDat %>%
  mutate('block' = ((as.numeric(session) - 1) %/% 5) + 1)

####  Save and Load Data ####

# save(pokeStateDat, file = 'pokeStateDat.RData')

# load('pokeStateDat.RData')

pokeStateDat %>%
  left_join(genoDat) %>%
  saveRDS(here("data", "PI_timing_poke_state.RDS"))
