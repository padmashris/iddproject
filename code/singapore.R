# load libraries 
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
library(EpiEstim)
library(R0)
library(reshape2)

# set working directory
setwd("~/Desktop/IDD/Project_COVID")

# Load the data
data <- read_csv("~/Desktop/IDD/Project_COVID/singapore_cleaned.csv")

# Convert the date variable to a Date object
data$date_confirmation <- as.Date(data$date_confirmation, format = "%Y.%m.%d")
data$date_onset_symptoms <- as.Date(data$date_onset_symptoms, format = "%Y-%m-%d")

#### first method ####

# Calculate daily counts
daily_counts <- data %>%
  count(date_confirmation) %>%
  rename(date = date_confirmation, count = n)

N <- 5000000 # population at risk
C <- sum(daily_counts$count) # total cases

total_i <- 0
for(i in (N-C+1):(N-1)){
  total_i <- 1/i + total_i
}

R01 <- ((N-1)/C)*total_i

#### second method ####

# Calculate the exponential growth rate
growth_rate <- lm(log(count) ~ decimal_date(date), data = daily_counts)$coefficients[2]

# calculating serial interval 

#Load data with symptom onset 
spdata <- read_csv("~/Desktop/IDD/Project_COVID/COVID-19_Singapore_data_revised.csv", col_types = list(presumed_infected_date = col_datetime()))
#table(spdata$`Related cases`) # There is one cell with "\n", needs to be changed to 'NA'
spdata$`Related cases`[which(spdata$`Related cases` == "\n")] <- NA
# Rename columns 2, 3 and 4 so no spaces
spdata <- rename(spdata, related_cases = starts_with("Related"),
                 cluster_links = "Cluster links",
                 relationship_notes = starts_with("Relation"))

# Remove all the cases that do not have info on date of symptom onset 
spdata <- filter(spdata, !is.na(date_onset_symptoms)) 

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

# making a data frame with a row for every suspected infector-infectee pair - and including 
# the serial interval for this pair, and the incubation period of both infector and infectee. 

sing.data <- data.frame(infector = undir_dates$from, infectee = undir_dates$to, serial.interval = undir_dates$abs_serial_interval, inc.infector.min = spdata$minIncTimes[match(undir_dates$from,spdata$CaseID)], inc.infector.max = spdata$maxIncTimes[match(undir_dates$from,spdata$CaseID)], inc.infectee.min =  spdata$minIncTimes[match(undir_dates$to,spdata$CaseID)], inc.infectee.max = spdata$maxIncTimes[match(undir_dates$to,spdata$CaseID)])
sing.data = sing.data[!is.na(sing.data$serial.interval),] # filter out NAs

sing.data$serial.interval = as.numeric(sing.data$serial.interval)

fit.si <- generation.time("gamma",sing.data$serial.interval)

# Extract mean and coefficient of variation
mu <- fit.si$mean
nu <- fit.si$sd/fit.si$mean

# Calculate reproductive number using exponential growth rate - method
r <- growth_rate # example value
R <- (1 + r*mu*(nu^2))^(1/(nu^2))

cat("Estimated R0: ",R)

# third method

daily_counts <- daily_counts %>% mutate(day=as.numeric(difftime(daily_counts$date,daily_counts$date[1],
                                                                units="days")))
df <- daily_counts[,2:3]

# calculate R0 using White and Pagano #
R0.ML <- est.R0.ML(df$count,t=df$day,begin=1,end=34,GT=fit.si)
cat("Estimated R0: ",R0.ML$R)
R0.ML$R

# calculate R0 using sequential bayesian approach: 
R0.SB <- est.R0.SB(df$count,GT=fit.si)
cat("Estimated R0: ",R0.SB$R)

flu.R0.results <- as.data.frame(cbind(df[1:18,], c(NA,R0.SB$R)))
names(flu.R0.results) <- c("dayNum","N","R0.SB")
tmp <- melt(flu.R0.results,id="dayNum")
flu.R0.long <- dcast(tmp,dayNum~variable)
flu.R0.long <- as.data.frame(rbind(as.matrix(flu.R0.long[,1:3]),as.matrix(flu.R0.long[,c(1:2,4)])))
names(R0.long) <- c("dayNum","N","R0")
flu.R0.long$Method <- c(rep("SB",15),rep("ML",15))

flu.plot <- ggplot(flu.R0.long,aes(x=dayNum,y=N/2))+
  geom_bar(stat="identity")+labs(y="Number of cases",x="Day of outbreak",title="(a)")+
  geom_line(aes(y=R0*8,linetype=Method),size=1)+
  scale_y_continuous(sec.axis=sec_axis(trans=~./8,name=expression(hat(R)[0])))+
  geom_hline(yintercept=8,color="gray")
flu.plot
