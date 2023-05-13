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
data <- read.csv("~/Desktop/IDD/Project_COVID/SK_cleaned.csv")

data <- data %>% rename(date_onset_symptoms = "S1")

# Convert the date variable to a Date object
data$date_confirmation <- as.Date(data$date_confirmed, format = "%d/%m/%y")
data$date_onset_symptoms <- as.Date(data$date_onset_symptoms, format = "%d/%m/%y")

# Calculate daily counts
daily_counts <- data %>%
  count(date_confirmation) %>%
  rename(date = date_confirmation, count = n)

N <- 51840000 # population at risk - > need to change for SK 
C <- sum(daily_counts$count)
total_i <- 0
for(i in (N-C+1):(N-1)){
  total_i <- 1/i + total_i
}

R01 <- ((N-1)/C)*total_i
R01

# Calculate the exponential growth rate
growth_rate <- lm(log(count) ~ decimal_date(date), data = daily_counts)$coefficients[2]

#### calculating serial interval ####

skdata <- data
skdata <- skdata[which(!is.na(skdata$date_onset_symptoms)),] # removes 3 observations

skdata <- rename(skdata, end_source = CR)
skdata$end_source <- as.Date(skdata$end_source, format = "%d/%m/%y")
skdata <- rename(skdata, start_source = CL)
skdata$start_source <- as.Date(skdata$start_source, format = "%d/%m/%y")
skdata$S2 <- as.Date(skdata$S2, format = "%d/%m/%y")

# Rename columns 2, 3 and 4 so no spaces
skdata <- rename(skdata, related_cases = starts_with("Related"),
                 cluster_links = "Cluster.Links")

# Remove all the cases that do not have info on date of symptom onset 
skdata <- filter(skdata, !is.na(date_onset_symptoms)) 


# assume that, in each pair, the one who showed symptoms first was the infector
sknodes <- skdata$ID

# Split into separate columns
skdata <- separate(skdata,col = related_cases,into = paste("contactID", 1:5, sep = "_"),
                   fill = "right")
# Turn into numeric values
skdata <- mutate(skdata, 
                 contactID_1 = as.numeric(contactID_1),
                 contactID_2 = as.numeric(contactID_2),
                 contactID_3 = as.numeric(contactID_3),
                 contactID_4 = as.numeric(contactID_4),
                 contactID_5 = as.numeric(contactID_5))

# Select down to columns of interest 
skedges <- skdata[,c(1,9:13)]

# Remove rows with NAs for at least one contact
skedges <- filter(skedges, !is.na(skedges$contactID_1)) 
SK_edges = data.frame(from=2,to=1) 
for (n in 1:nrow(skedges)){
  for (k in 2:ncol(skedges)){ 
    if (!is.na(skedges[n,k])){
      SK_edges=rbind(SK_edges, c(skedges[[n,k]],skedges[[n,1]])) 
    }  
  }
}
SK_edges=SK_edges[-1,]

# create undirected graph by removing duplicates
undir=data.frame(from = pmin(SK_edges[,1],SK_edges[,2]),  
                 to = pmax(SK_edges[,1], SK_edges[,2]))
undir = unique(undir)
undir = undir[-which(undir[,1]==undir[,2]),]
skdata_sympt <- skdata[,c(1,4)] # selecting CaseID and date_onset_symptoms
names(skdata_sympt) <- str_replace(names(skdata_sympt), "ID", "from")
undir_dates <- left_join(undir, skdata_sympt, by = "from")
names(undir_dates) <- str_replace(names(undir_dates), "date_onset_symptoms", "from_sympt_date")
names(skdata_sympt) <- str_replace(names(skdata_sympt), "from", "to")
undir_dates <- left_join(undir_dates, skdata_sympt, by = "to")
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
skdata$minIncTimes <- skdata$date_onset_symptoms - skdata$end_source
skdata$maxIncTimes <- skdata$date_onset_symptoms - skdata$start_source
skdata$maxIncTimes = pmax(3, skdata$maxIncTimes)
skdata$minIncTimes = pmax(1, skdata$minIncTimes)
skdata$serial_interval <- abs(skdata$date_onset_symptoms - skdata$S2)

# making a data frame with a row for every suspected infector-infectee pair - and including 
# the serial interval for this pair, and the incubation period of both infector and infectee. 
#???
SK.data <- data.frame(infector = skdata$Infector_ID, infectee = skdata$ID, serial.interval = skdata$serial_interval,
                      inc.infector.min = skdata$minIncTimes[match(skdata$infector_ID,skdata$ID)], inc.infector.max = skdata$maxIncTimes[match(skdata$Infector_ID,skdata$ID)])
SK.data = sing.data[!is.na(sing.data$serial.interval),] # filter out NAs


sing.data$serial.interval = as.numeric(sing.data$serial.interval)

library(fitdistrplus)
fit2 <- fitdist(skdata$serial_interval,distr="pois")
# is removing 0s the right way? 

fit.sk <- generation.time("gamma",as.numeric(skdata$serial_interval))

# Extract mean and coefficient of variation
mu <- fit.sk$mean
nu <- fit.sk$sd/fit.sk$mean

# Calculate reproductive number using exponential growth rate - method
r <- growth_rate # example value
R <- (1 + r*mu*nu^2)^0.5

cat("Estimated R0: ",R)

##### third method #####

daily_counts <- daily_counts %>% mutate(day=as.numeric(difftime(daily_counts$date,daily_counts$date[1],
                                                                units="days")))
df <- daily_counts[,2:3]

# calculate R0 using White and Pagano #
R0.ML <- est.R0.ML(df$count,t=df$day,begin=1,end=34,GT=fit.sk)
cat("Estimated R0: ", R0.ML$R)

# calculate R0 using sequential bayesian approach: 
R0.SB <- est.R0.SB(df$count,GT=fit.sk)
cat("Estimated R0: ",mean(R0.SB$R))
