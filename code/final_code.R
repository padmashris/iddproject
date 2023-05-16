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

setwd("~/Desktop/IDD/Project_COVID")


# Singapore ---------------------------------------------------------------

# Load the data
data <- read_csv("~/Desktop/IDD/Project_COVID/singapore_cleaned.csv")

# Convert the date variable to a Date object
data$date_confirmation <- as.Date(data$date_confirmation, format = "%Y.%m.%d")
data$date_onset_symptoms <- as.Date(data$date_onset_symptoms, format = "%Y-%m-%d")

##### first method - final epidemic size #####

# Calculate daily counts
daily_counts <- data %>%
  count(date_confirmation) %>%
  rename(date = date_confirmation, count = n)

sg_plotinc <- daily_counts

N <- 5000000 # population at risk
C <- sum(daily_counts$count) # total cases

total_i <- 0
for(i in (N-C+1):(N-1)){
  total_i <- 1/i + total_i
}

R01.sg <- ((N-1)/C)*total_i

#### second method - exponential growth rate ####

# calculate the exponential growth rate using a regression model
growth_rate <- lm(log(count) ~ decimal_date(date), data = daily_counts)$coefficients[2]

# calculating serial interval 
# all serial interval calculations from line list data are based on Tindale et al., 2020

# load data with symptom onset 
data_sg <- read_csv("~/Desktop/IDD/Project_COVID/COVID-19_Singapore_data_revised.csv", col_types = list(presumed_infected_date = col_datetime()))
#table(spdata$`Related cases`)
data_sg$`Related cases`[which(data_sg$`Related cases` == "\n")] <- NA

data_sg <- rename(data_sg, related_cases = starts_with("Related"),
                  cluster_links = "Cluster links")
data_sg <- filter(data_sg, !is.na(date_onset_symptoms)) 

# we assume that in each pair, the person who exhibited symptoms first is considered the one who transmitted the infection.
nodes.sg <- data_sg$CaseID

# separating the contacts each case had into columns
data_sg <- separate(data_sg,col = related_cases,into = paste("contactID", 1:7, sep = "_"),
                    fill = "right")
data_sg <- mutate(data_sg, 
                  contactID_1 = as.numeric(contactID_1),
                  contactID_2 = as.numeric(contactID_2),
                  contactID_3 = as.numeric(contactID_3),
                  contactID_4 = as.numeric(contactID_4),
                  contactID_5 = as.numeric(contactID_5),
                  contactID_6 = as.numeric(contactID_6),
                  contactID_7 = as.numeric(contactID_7))

edges.sg <- data_sg[,1:8] # selects the contact columns in the dataset
edges.sg <- filter(edges.sg, !is.na(edges.sg$contactID_1)) 

SG_edges = data.frame(from=3,to=4) # create dummy edge in data frame

# fill the data frame with 'edges' as contacts for each case, and 'node' as each case
for (n in 1:nrow(edges.sg)) {
  for (k in 2:ncol(edges.sg)) { 
    if (!is.na(edges.sg[n,k])) {
      SG_edges=rbind(SG_edges, c(edges.sg[[n,k]],edges.sg[[n,1]])) 
    }  
  }
}

SG_edges=SG_edges[-1,] # removing the initial dummy edge

# remove duplicates and ensures edges aren't repeated
sg_undir=data.frame(from = pmin(SG_edges[,1],SG_edges[,2]),  
                    to = pmax(SG_edges[,1], SG_edges[,2]))
sg_undir = unique(sg_undir)
sg_undir = sg_undir[-which(sg_undir[,1]==sg_undir[,2]),]

data_sg_sympt <- data_sg[,c(1,22)] # selecting CaseID and date_onset_symptoms

# cleans up the data frame: rename columns, add info on symptom onset data for both 'from' and 'to' nodes
names(data_sg_sympt) <- str_replace(names(data_sg_sympt), "CaseID", "from")
sg_undir_dates <- left_join(sg_undir, data_sg_sympt, by = "from")
names(sg_undir_dates) <- str_replace(names(sg_undir_dates), "date_onset_symptoms", "from_sympt_date")
names(data_sg_sympt) <- str_replace(names(data_sg_sympt), "from", "to")
sg_undir_dates <- left_join(sg_undir_dates, data_sg_sympt, by = "to")
names(sg_undir_dates) <- str_replace(names(sg_undir_dates), "date_onset_symptoms", "to_sympt_date")
sg_undir_dates <- mutate(sg_undir_dates, earliest_sympt_onset = pmin(to_sympt_date, from_sympt_date, na.rm = T), 
                         raw_serial_interval = to_sympt_date - from_sympt_date,  
                         abs_serial_interval = abs(raw_serial_interval)) # serial interval calculation

# groups positive and negative serial interval values
pos <- filter(sg_undir_dates, raw_serial_interval >= 0)
neg <- filter(sg_undir_dates, raw_serial_interval < 0)
onlyone <- filter(sg_undir_dates, is.na(raw_serial_interval))

# renaming columns
names(neg)[1] <- "to"
names(neg)[2] <- "from"
names(neg)[3] <- "to_sympt_date"
names(neg)[4] <- "from_sympt_date"

sg_undir_dates <- bind_rows(pos, neg, onlyone)
sg_undir_dates$pto <- str_pad(sg_undir_dates$to, width = 2, side = "left", pad = "0")
sg_undir_dates$pfrom <- str_pad(sg_undir_dates$from, width = 2, side = "left", pad = "0")
sg_undir_dates <- mutate(sg_undir_dates, pairID = factor(paste("case", pfrom, "-", "case", pto, sep = "")))
rm(pos, neg, onlyone)

# calculating incubation periods
data_sg$minIncTimes <- data_sg$date_onset_symptoms - data_sg$end_source
data_sg$maxIncTimes <- data_sg$date_onset_symptoms - data_sg$start_source
data_sg$maxIncTimes = pmax(3, data_sg$maxIncTimes)
data_sg$minIncTimes = pmax(1, data_sg$minIncTimes)

# cleaning and storing the final dataset

data_sgclean <- data.frame(infector = sg_undir_dates$from, infectee = sg_undir_dates$to, serial.interval = sg_undir_dates$abs_serial_interval, inc.infector.min = data_sg$minIncTimes[match(sg_undir_dates$from,data_sg$CaseID)], inc.infector.max = data_sg$maxIncTimes[match(sg_undir_dates$from,data_sg$CaseID)], inc.infectee.min =  data_sg$minIncTimes[match(sg_undir_dates$to,data_sg$CaseID)], inc.infectee.max = data_sg$maxIncTimes[match(sg_undir_dates$to,data_sg$CaseID)])
data_sgclean = data_sgclean[!is.na(data_sgclean$serial.interval),] # filter out NAs
data_sgclean$serial.interval = as.numeric(data_sgclean$serial.interval)

#Pivot the to/from dates of symptom onset column to a long format, so that can make a legend based on this variable
undir_dotplot <- pivot_longer(sg_undir_dates, 
                              cols = contains("sympt_date"),
                              names_to = "pair_member",
                              values_to = "onset_date")

#Let's rename the values so it makes more sense in the legend
undir_dotplot$pair_member <- str_replace(undir_dotplot$pair_member, pattern = "from_sympt_date", replacement = "Presumed infector")
undir_dotplot$pair_member <- str_replace(undir_dotplot$pair_member, pattern = "to_sympt_date", replacement = "Presumed infectee")

#Make the Cleveland dotplot
s <- ggplot(undir_dotplot, aes(y = reorder(pairID, earliest_sympt_onset))) +
  geom_segment(aes(x = earliest_sympt_onset, xend = earliest_sympt_onset + abs_serial_interval, yend = pairID), 
               color = "#404788FF") +
  geom_point(aes(x = onset_date, color = pair_member, fill = pair_member, shape = pair_member)) +
  scale_x_date(date_breaks = "1 day") +
  scale_color_manual(name = "Pair member for \ndate of symptom onset", values = c("#D44842FF", "#FAC127FF")) +
  scale_fill_manual(name = "Pair member for \ndate of symptom onset", values = c("#D44842FF", "#FAC127FF")) +
  scale_shape_manual(name = "Pair member for \ndate of symptom onset", values = c(23, 21)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey80", linetype = "dashed"),
        panel.background = element_rect(fill = "white")) +
  labs(title = "Serial intervals of possible case pairs from Singapore, over time",
       x = "Date of symptom onset",
       y = "Case pairs")
s

# fitting a gamma distribution to the serial interval data
fit.si <- generation.time("gamma",data_sgclean$serial.interval)

# store mean and coefficient of variation
# method of calculation follows from Nishiura et al., 2009
mu <- fit.si$mean
nu <- fit.si$sd/fit.si$mean

# calculate reproductive number using exponential growth rate - method
r <- growth_rate
R.exp.sg <- (1 + r*mu*(nu^2))^(1/(nu^2))
R.exp.sg

##### MLE and SBE #####

# defining the MLE method function

# set an arbitrary range of R0 values to optimize from
range = c(0.01,50)

# check.incid() formats incidence data for further estimation
epid <- check.incid(df$count, t=df$day, time.step=1)

# converting dates to number of days since first case
daily_counts <- daily_counts %>% mutate(day=as.numeric(difftime(daily_counts$date,daily_counts$date[1],
                                                                units="days")))
df <- daily_counts[,2:3]

# estimating R0 using the maximum likelihood function as developed by White and Pagano, 2009 and implemented by Boelle and Obadia, 2011.

# inside the function below: 
## fit.epid() computes the Poisson log-likelihood of the epidemic and predicts the curve
# res.R <- optimize(f=fit.epid, log(range),GT=GT,epid=epid,import=NULL,maximum=TRUE) 
# pred <- fit.epid(res.R$maximum,epid,GT,import=import,pred=TRUE)
# returns the estimate of the reproduction number using MLE and its confidence interval

# calculate R0 using White and Pagano 
R0.ML.sg <- R0::est.R0.ML(df$count,t=df$day,begin=1,end=34,GT=fit.si)
R0.ML.sg

# calculate R0 using sequential bayesian approach: 
## we assume the inital R to be following a unfiorm distribution
R0.SB.sg <- est.R0.SB(df$count,GT=fit.si)
mean(R0.SB.sg$R)


# South Korea -------------------------------------------------------------

# load the data
data <- read.csv("~/Desktop/IDD/Project_COVID/SK_cleaned.csv")

data <- data %>% rename(date_onset_symptoms = "S1")

##### first method - final epidemic size #####
# Convert the date variable to a Date object
data$date_confirmation <- as.Date(data$date_confirmed, format = "%d/%m/%y")
data$date_onset_symptoms <- as.Date(data$date_onset_symptoms, format = "%d/%m/%y")

# calculate daily counts and format data
daily_counts <- data %>%
  count(date_confirmation) %>%
  rename(date = date_confirmation, count = n)

sk_plotinc <- daily_counts

N <- 51840000 # total population of south korea in 2020
C <- sum(daily_counts$count)
total_i <- 0
for(i in (N-C+1):(N-1)){
  total_i <- 1/i + total_i
}

R01.sk <- ((N-1)/C)*total_i

##### second method - exponential growth rate ##### 

growth_rate <- lm(log(count) ~ decimal_date(date), data = daily_counts)$coefficients[2]

# calculating serial interval 

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
sk_dates_und <- left_join(undir, skdata_sympt, by = "from")
names(sk_dates_und) <- str_replace(names(sk_dates_und), "date_onset_symptoms", "from_sympt_date")
names(skdata_sympt) <- str_replace(names(skdata_sympt), "from", "to")
sk_dates_und <- left_join(sk_dates_und, skdata_sympt, by = "to")
names(sk_dates_und) <- str_replace(names(sk_dates_und), "date_onset_symptoms", "to_sympt_date")
sk_dates_und <- mutate(sk_dates_und, earliest_sympt_onset = pmin(to_sympt_date, from_sympt_date, na.rm = T), 
                       raw_serial_interval = to_sympt_date - from_sympt_date,   #5 NAs because only 1 case in the pair has a date of symptom onset
                       abs_serial_interval = abs(raw_serial_interval))
pos <- filter(sk_dates_und, raw_serial_interval >= 0)
neg <- filter(sk_dates_und, raw_serial_interval < 0)
onlyone <- filter(sk_dates_und, is.na(raw_serial_interval))
names(neg)
names(neg)[1] <- "to"
names(neg)[2] <- "from"
names(neg)[3] <- "to_sympt_date"
names(neg)[4] <- "from_sympt_date"
names(neg)
sk_dates_und <- bind_rows(pos, neg, onlyone)
sk_dates_und$pto <- str_pad(sk_dates_und$to, width = 2, side = "left", pad = "0")
sk_dates_und$pfrom <- str_pad(sk_dates_und$from, width = 2, side = "left", pad = "0")
sk_dates_und <- mutate(sk_dates_und, pairID = factor(paste("case", pfrom, "-", "case", pto, sep = "")))
rm(pos, neg, onlyone)

# calculating incubation periods
skdata$minIncTimes <- skdata$date_onset_symptoms - skdata$end_source
skdata$maxIncTimes <- skdata$date_onset_symptoms - skdata$start_source
skdata$maxIncTimes = pmax(3, skdata$maxIncTimes)
skdata$minIncTimes = pmax(1, skdata$minIncTimes)
skdata$serial_interval <- abs(skdata$date_onset_symptoms - skdata$S2)

#Pivot the to/from dates of symptom onset column to a long format, so that can make a legend based on this variable
undir_dotplot <- pivot_longer(sk_dates_und, 
                              cols = contains("sympt_date"),
                              names_to = "pair_member",
                              values_to = "onset_date")

#Let's rename the values so it makes more sense in the legend
undir_dotplot$pair_member <- str_replace(undir_dotplot$pair_member, pattern = "from_sympt_date", replacement = "Presumed infector")
undir_dotplot$pair_member <- str_replace(undir_dotplot$pair_member, pattern = "to_sympt_date", replacement = "Presumed infectee")

#Make the Cleveland dotplot
kor <- ggplot(undir_dotplot, aes(y = reorder(pairID, earliest_sympt_onset))) +
  geom_segment(aes(x = earliest_sympt_onset, xend = earliest_sympt_onset + abs_serial_interval, yend = pairID), 
               color = "#404788FF") +
  geom_point(aes(x = onset_date, color = pair_member, fill = pair_member, shape = pair_member)) +
  scale_x_date(date_breaks = "1 day") +
  scale_color_manual(name = "Pair member for \ndate of symptom onset", values = c("#D44842FF", "#FAC127FF")) +
  scale_fill_manual(name = "Pair member for \ndate of symptom onset", values = c("#D44842FF", "#FAC127FF")) +
  scale_shape_manual(name = "Pair member for \ndate of symptom onset", values = c(23, 21)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey80", linetype = "dashed"),
        panel.background = element_rect(fill = "white")) +
  labs(title = "Serial intervals of possible case pairs from South Korea, over time",
       x = "Date of symptom onset",
       y = "Case pairs")
kor

fit.sk <- generation.time("gamma",as.numeric(skdata$serial_interval))

# Extract mean and coefficient of variation
mu <- fit.sk$mean
nu <- fit.sk$sd/fit.sk$mean

# Calculate reproductive number using exponential growth rate - method
r <- growth_rate 
R.sk <- (1 + r*mu*nu^2)^0.5

##### MLE and SBE #####

daily_counts <- daily_counts %>% mutate(day=as.numeric(difftime(daily_counts$date,daily_counts$date[1],
                                                                units="days")))
df <- daily_counts[,2:3]


# estimating R0 using the maximum likelihood function as developed by White and Pagano, 2009 and implemented by Boelle and Obadia, 2011.

# inside the function below: 
## fit.epid() computes the Poisson log-likelihood of the epidemic and predicts the curve
# res.R <- optimize(f=fit.epid, log(range),GT=GT,epid=epid,import=NULL,maximum=TRUE) 
# pred <- fit.epid(res.R$maximum,epid,GT,import=import,pred=TRUE)
# returns the estimate of the reproduction number using MLE and its confidence interval

R0.ML.sk <- R0::est.R0.ML(df$count,t=df$day,begin=1,end=34,GT=fit.sk)

# calculate R0 using sequential bayesian approach: 
R0.SB.sk <- est.R0.SB(df$count,GT=fit.sk)

# Tianjin -----------------------------------------------------------------

# Load the data
data <- read_csv("Tianjin135cases_revised.csv")

#### first method - final epidemic size ####

# Convert the date variable to a Date object
data$date_confirmation <- as.Date(data$confirm_date, format = "%d/%m/%Y")
data$date_onset_symptoms <- as.Date(data$symptom_onset, format = "%d/%m/%Y")

# Calculate daily counts
daily_counts <- data %>%
  count(date_confirmation) %>%
  rename(date = date_confirmation, count = n)

tj_plotinc <- daily_counts

N <- 15000000 # population at risk
C <- sum(daily_counts$count)
total_i <- 0
for(i in (N-C+1):(N-1)){
  total_i <- 1/i + total_i
}

R01.tj <- ((N-1)/C)*total_i

#### second method - exponential growth rate ####
# all serial interval calculations from line list data are based on Tindale et al., 2020

# Calculate the exponential growth rate
growth_rate <- lm(log(count) ~ decimal_date(date), data = daily_counts)$coefficients[2]

tj_data = read.csv("Tianjin135cases_revised.csv",na.strings = "", stringsAsFactors = F)
tj_data$symptom_onset=as.Date(tj_data$symptom_onset, format = "%d/%m/%Y")
tj_data$start_source=as.Date(tj_data$start_source, format = "%d/%m/%Y")
tj_data$end_source=as.Date(tj_data$end_source,format = "%d/%m/%Y" )
tj_data$confirm_date=as.Date(tj_data$confirm_date,format = "%d/%m/%Y" )

# remove any missing data of symptom onset
tj_data <- tj_data[which(!is.na(tj_data$symptom_onset)),]  #removes 10 observations

# convert to lowercase and remove whitespace
tj_data$Infection_source <- str_to_lower(tj_data$Infection_source)
tj_data$Infection_source <- str_trim(tj_data$Infection_source)

# store in duplicate column 
tj_data$Infection_source_dup <- tj_data$Infection_source
unique(tj_data$Infection_source_dup)

# changing 'person' to 'individual'
tj_data$Infection_source_dup <- str_replace_all(tj_data$Infection_source_dup, pattern = "person", replacement = "individual")

# renaming source to make the contact clearer
tj_data$Infection_source_dup <- str_replace(tj_data$Infection_source_dup, 
                                            pattern = "coworker of a individual from wuhan",
                                            replacement = "coworker")

# making the clusters more clear
tj_data <- mutate(tj_data, source_group = case_when(!is.na(str_match(Infection_source_dup, "wuhan|hubei")) ~ "Wuhan and Hubei", 
                                                    !is.na(str_match(Infection_source_dup, "mall|store|shopper|shopping")) ~ "Mall", 
                                                    !is.na(str_match(Infection_source_dup, "family|relative|wife|mother|son|sister|daughter|brother|husband|duaghtor|whife|hunsband")) ~ "Relative", 
                                                    !is.na(str_match(Infection_source_dup, "coworker|business|workplace|colleague|colleage")) ~ "Coworker",
                                                    !is.na(str_match(Infection_source_dup, "tj|patient")) ~ "Other relationship", 
                                                    !is.na(str_match(Infection_source_dup, "train|travel|trip|hebei|dalian")) ~ "Other travel", 
                                                    !is.na(str_match(Infection_source_dup, "unknown|unclear")) ~ "Unknown", 
                                                    is.na(Infection_source_dup) ~ "Unknown",
                                                    T ~ "other"))

tjnodes <- tj_data$case_id

# Transform everything to lower case to make sure there aren't any issues with matching due to case inconsistencies
tjnodes <- str_to_lower(tjnodes) 
tj_data$case_id <- str_to_lower(tj_data$case_id)

edges = data.frame(from=tjnodes[9],to=tjnodes[21],stringsAsFactors = F) # i read this one manually 

for (id in 1:nrow(tj_data)) {
  tonode=tj_data$case_id[id]
  fromnodes=str_extract_all(tj_data$Infection_source[id], "tj\\d+", simplify = T) #in lower case due to above early/late split on infection source
  if (length(fromnodes)>0) {
    for (k in 1:length(fromnodes)) {
      edges=rbind(edges, c(fromnodes[k], tonode))
    }
  }
}
head(edges)
edges=edges[-1,] #Remove the initial relationship we gave so it isn't duplicated
edges=edges[-which(is.na(edges[,1])),] # NAs arose from a few empty entries for Infection_source 

# Make a smaller dataset of original tpdata that contains only the CaseID and date of symptom onset
tj_data_sympt <- tj_data[,c(1,4)]

# renaming columns and joining the edges 
names(tj_data_sympt) <- str_replace(names(tj_data_sympt), "case_id", "from")
tj_undir_dates <- left_join(edges, tj_data_sympt, by = "from")
names(tj_undir_dates) <- str_replace(names(tj_undir_dates), "symptom_onset", "from_sympt_date")
names(tj_data_sympt) <- str_replace(names(tj_data_sympt), "from", "to")
tj_undir_dates <- left_join(tj_undir_dates, tj_data_sympt, by = "to")
names(tj_undir_dates) <- str_replace(names(tj_undir_dates), "symptom_onset", "to_sympt_date")

# calculating raw serial interval
tj_undir_dates <- mutate(tj_undir_dates, earliest_sympt_onset = pmin(to_sympt_date, from_sympt_date, na.rm = T), 
                         raw_serial_interval = to_sympt_date - from_sympt_date,   
                         abs_serial_interval = abs(raw_serial_interval))

# splitting dataset into positive and negative serial interval values and binding back values
pos <- filter(tj_undir_dates, raw_serial_interval >= 0)
neg <- filter(tj_undir_dates, raw_serial_interval < 0)
onlyone <- filter(tj_undir_dates, is.na(raw_serial_interval)) #3 NAs where date of symptom onset is not known; keep for making columns for now

names(neg)
names(neg)[1] <- "to"
names(neg)[2] <- "from"
names(neg)[3] <- "to_sympt_date"
names(neg)[4] <- "from_sympt_date"
tj_undir_dates <- bind_rows(pos, neg, onlyone)

# prepare to plot
# adding a column of caseID numbers for each edge to and from
# adding a zero to the left of the number so all numbers have three digits  
tj_undir_dates$pto <- str_replace(tj_undir_dates$to, pattern = "tj", replacement = "")
tj_undir_dates$pto <- str_pad(tj_undir_dates$pto, width = 3, side = "left", pad = "0")

tj_undir_dates$pfrom <- str_replace(tj_undir_dates$from, pattern = "tj", replacement = "")
tj_undir_dates$pfrom <- str_pad(tj_undir_dates$pfrom, width = 3, side = "left", pad = "0")

# making a new column with case pair ID
tj_undir_dates <- mutate(tj_undir_dates, pairID = factor(paste("tj", pfrom, "-", "tj", pto, sep = "")))

rm(pos, neg, onlyone)

###  dotplot of raw serial intervals for each case pair

tj_undir_dates <- filter(tj_undir_dates, !is.na(raw_serial_interval)) # filter out NAs

# turning data from wide to long format
undir_dotplot <- pivot_longer(tj_undir_dates, 
                              cols = contains("sympt_date"),
                              names_to = "pair_member",
                              values_to = "onset_date")

# rename variables
undir_dotplot$pair_member <- str_replace(undir_dotplot$pair_member, pattern = "from_sympt_date", replacement = "Presumed infector")
undir_dotplot$pair_member <- str_replace(undir_dotplot$pair_member, pattern = "to_sympt_date", replacement = "Presumed infectee")

# running the Cleveland dotplot

t <- ggplot(undir_dotplot, aes(y = reorder(pairID, earliest_sympt_onset))) +
  geom_segment(aes(x = earliest_sympt_onset, xend = earliest_sympt_onset + abs_serial_interval, yend = pairID), 
               color = "#404788FF") +
  geom_point(aes(x = onset_date, color = pair_member, fill = pair_member, shape = pair_member)) +
  scale_x_date(date_breaks = "1 day") +
  scale_color_manual(name = "Pair member for \ndate of symptom onset", values = c("#D44842FF", "#FAC127FF")) +
  scale_fill_manual(name = "Pair member for \ndate of symptom onset", values = c("#D44842FF", "#FAC127FF")) +
  scale_shape_manual(name = "Pair member for \ndate of symptom onset", values = c(23, 21)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey80", linetype = "dashed"),
        panel.background = element_rect(fill = "white")) +
  labs(title = "Serial intervals of possible case pairs from Tianjin, over time",
       x = "Date of symptom onset",
       y = "Case pairs")
t

# arranging by the date of earliest symptom onset for each case pair
tj_undir_dates <- arrange(tj_undir_dates, earliest_sympt_onset)

fit.ti <- generation.time("gamma",as.numeric(tj_undir_dates$raw_serial_interval))

# Extract mean and coefficient of variation
mu <- fit.ti$mean
nu <- fit.ti$sd/fit.ti$mean

# Calculate reproductive number using exponential growth rate - method
r <- growth_rate # -8.05
R.tj <- (1 + r*mu*(nu^2))^(1/(nu^2)) ## a non value - due to negative exponential growth rate!! 

##### third method #####

tj_counts <- daily_counts %>% mutate(day=as.numeric(difftime(daily_counts$date,daily_counts$date[1],
                                                             units="days")))

df <- tj_counts[,2:3] # selecting the columns we need (count and day)

# estimating R0 using the maximum likelihood function as developed by White and Pagano, 2009 and implemented by Boelle and Obadia, 2011.

# inside the function below: 
## fit.epid() computes the Poisson log-likelihood of the epidemic and predicts the curve
# res.R <- optimize(f=fit.epid, log(range),GT=GT,epid=epid,import=NULL,maximum=TRUE) 
# pred <- fit.epid(res.R$maximum,epid,GT,import=import,pred=TRUE)
# returns the estimate of the reproduction number using MLE and its confidence interval

# calculate R0 using White and Pagano #
R0.ML.tj <- R0::est.R0.ML(df$count,t=df$day,begin=1,end=32,GT=fit.si)

# calculate R0 using sequential bayesian approach: 
R0.SB.tj <- est.R0.SB(df$count,GT=fit.si)
mean(R0.SB.tj$R)

ggarrange(s,kor,p,
          ncol=1,nrow=3,labels=c("a","b","c"))
