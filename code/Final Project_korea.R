library(tidyverse)
library(lubridate)
library(incidence)
library(MASS)
library(survminer)
library(survival)
library(tidyverse)
library(lubridate)
library(icenReg)
library(igraph)
library(visNetwork)
library(mvtnorm)
library(viridis)
options(digits=3)
library(ggplot2)

# Load the data
data <- read_csv("Documents/JHU/term4/ID dynamics/final project/SK_cleaned.csv")

# Convert the date variable to a Date object
data$date_confirmation <- as.Date(data$date_confirmation, format = "%Y/%m/%d")
data$date_onset_symptoms <- as.Date(data$date_onset_symptoms, format = "%Y/%m/%d")

# Calculate daily counts
daily_counts <- data %>%
  count(date_confirmation) %>%
  rename(date = date_confirmation, count = n)

N <- 50000000 # population at risk
C <- sum(daily_counts$count)
total_i <- 0
for(i in (N-C+1):(N-1)){
  total_i <- 1/i + total_i
}

R01 <- ((N-1)/C)*total_i

# Calculate the exponential growth rate
growth_rate <- lm(log(count) ~ decimal_date(date), data = daily_counts)$coefficients[2]

#### calculating serial interval ####

#Load data with symptom onset 
spdata <- read_csv("Documents/JHU/term4/ID dynamics/final project/SK_revised.csv", 
                   col_types = list(presumed_infected_date = col_datetime(format = "%Y/%m/%d")))
spdata <- rename(spdata, date_onset_symptoms = S1)
spdata$date_onset_symptoms <- as.Date(spdata$date_onset_symptoms, format = "%Y/%m/%d")
spdata <- rename(spdata, end_source = CR)
spdata$end_source <- as.Date(spdata$end_source, format = "%Y/%m/%d")
spdata <- rename(spdata, start_source = CL)
spdata$start_source <- as.Date(spdata$start_source, format = "%Y/%m/%d")
spdata$S2 <- as.Date(spdata$S2, format = "%Y/%m/%d")
###delete？？
#table(spdata$`Related cases`) # There is one cell with "\n", needs to be changed to 'NA'
spdata$`Related cases`[which(spdata$`Related cases` == "\n")] <- NA
# Rename columns 2, 3 and 4 so no spaces
spdata <- rename(spdata, related_cases = starts_with("Related"),
                 cluster_links = "Cluster links",
                 relationship_notes = starts_with("Relation"))
# Remove all the cases that do not have info on date of symptom onset 
spdata <- filter(spdata, !is.na(date_onset_symptoms)) 
###I have informations of infector
# assume that, in each pair, the one who showed symptoms first was the infector
spnodes <- spdata$CaseID
# Split into separate columns
spdata <- separate(spdata,col = related_cases,into = paste("contactID", 1:7, sep = "_"),
                   fill = "right")
# Turn into numeric values
spdata <- mutate(spdata, 
                 contactID_1 = as.numeric(contactID_1),
                 contactID_2 = as.numeric(contactID_2),
                 contactID_3 = as.numeric(contactID_3),
                 contactID_4 = as.numeric(contactID_4),
                 contactID_5 = as.numeric(contactID_5),
                 contactID_6 = as.numeric(contactID_6),
                 contactID_7 = as.numeric(contactID_7))
# Select down to columns of interest 
spedges <- spdata[,1:8]
# Remove rows with NAs for at least one contact
spedges <- filter(spedges, !is.na(spedges$contactID_1)) 
singedges = data.frame(from=2,to=1) 
for (n in 1:nrow(spedges)) {
  for (k in 2:ncol(spedges)) { 
    if (!is.na(spedges[n,k])) {
      singedges=rbind(singedges, c(spedges[[n,k]],spedges[[n,1]])) 
    }  
  }
}
singedges=singedges[-1,]
# create undirected graph by removing duplicates
undir=data.frame(from = pmin(singedges[,1],singedges[,2]),  
                 to = pmax(singedges[,1], singedges[,2]))
undir = unique(undir)
undir = undir[-which(undir[,1]==undir[,2]),]
spdata_sympt <- spdata[,c(1,22)] # selecting CaseID and date_onset_symptoms
names(spdata_sympt) <- str_replace(names(spdata_sympt), "CaseID", "from")
undir_dates <- left_join(undir, spdata_sympt, by = "from")
names(undir_dates) <- str_replace(names(undir_dates), "date_onset_symptoms", "from_sympt_date")
names(spdata_sympt) <- str_replace(names(spdata_sympt), "from", "to")
undir_dates <- left_join(undir_dates, spdata_sympt, by = "to")
names(undir_dates) <- str_replace(names(undir_dates), "date_onset_symptoms", "to_sympt_date")
undir_dates <- mutate(undir_dates, earliest_sympt_onset = pmin(to_sympt_date, from_sympt_date, na.rm = T), 
                      raw_serial_interval = to_sympt_date - from_sympt_date,   #5 NAs because only 1 case in the pair has a date of symptom onset
                      abs_serial_interval = abs(raw_serial_interval))
pos <- filter(undir_dates, raw_serial_interval >= 0)
neg <- filter(undir_dates, raw_serial_interval < 0)
onlyone <- filter(undir_dates, is.na(raw_serial_interval))
names(neg)
names(neg)[1] <- "to"
names(neg)[2] <- "from"
names(neg)[3] <- "to_sympt_date"
names(neg)[4] <- "from_sympt_date"
names(neg)
undir_dates <- bind_rows(pos, neg, onlyone)
undir_dates$pto <- str_pad(undir_dates$to, width = 2, side = "left", pad = "0")
undir_dates$pfrom <- str_pad(undir_dates$from, width = 2, side = "left", pad = "0")
undir_dates <- mutate(undir_dates, pairID = factor(paste("case", pfrom, "-", "case", pto, sep = "")))
rm(pos, neg, onlyone)

# calculating incubation periods
spdata$minIncTimes <- spdata$date_onset_symptoms - spdata$end_source
spdata$maxIncTimes <- spdata$date_onset_symptoms - spdata$start_source
spdata$maxIncTimes = pmax(3, spdata$maxIncTimes)
spdata$minIncTimes = pmax(1, spdata$minIncTimes)
spdata$serial_interval <- spdata$date_onset_symptoms - spdata$S2

# making a data frame with a row for every suspected infector-infectee pair - and including 
# the serial interval for this pair, and the incubation period of both infector and infectee. 
#???
sing.data <- data.frame(infector = spdata$Infector_ID, infectee = spdata$ID, serial.interval = spdata$serial_interval, inc.infector.min = spdata$minIncTimes[match(spdata$infector_ID,spdata$ID)], inc.infector.max = spdata$maxIncTimes[match(spdata$Infector_ID,spdata$ID)])
sing.data = sing.data[!is.na(sing.data$serial.interval),] # filter out NAs

sing.data$serial.interval = as.numeric(sing.data$serial.interval)

library(fitdistrplus)
fit2 <- fitdist(sing.data$serial.interval,distr="pois")
# is removing 0s the right way? 

# Extract mean and coefficient of variation
mu <- fit2$estimate[1]
nu <- 1/sqrt(fit2$estimate[1])

# Calculate reproductive number using exponential growth rate - method
r <- growth_rate # example value
R <- (1 + r*mu*nu^2)^0.5

cat("Estimated R0: ",R)

